// sir2d_mpi_omp.c — SIR 2D (reacción-difusión) híbrido MPI + OpenMP, halos periódicos
// Compilar: mpicc -O3 -march=native -fopenmp -o sir2d_mpi_omp sir2d_mpi_omp.c -lm
//
// Ejemplo:
//   mpirun -np 4 ./sir2d_mpi_omp --nx 1024 --ny 1024 --beta 0.6 --gamma 0.08 --D 0.2 \
//       --dt 0.1 --tmax 200 --snap_every 50 --out_prefix out/sir2d --i0frac 0.01 --sigma 3 --threads 8
//
// Notas:
//  - Requiere que nx % Px == 0 y ny % Py == 0 (Px*Py = #ranks).
//  - Guarda: out/sir2d_I_tXXXXX_rankRRRRR.csv y out/sir2d_R_tXXXXX_rankRRRRR.csv
//  - MLUPS = million lattice updates per second (global).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <sys/stat.h>
#include <errno.h>

typedef struct {
    int NX, NY;           // tamaño global
    double beta, gamma;
    double D;
    double dt, tmax;
    int snap_every;
    const char* out_prefix;
    double i0frac, sigma;
    int threads;
    const char* centers;  // lista de centros globales "x,y;x,y;..." (opcional)
} Params;

static void usage(const char* p){
    if (!p) p="sir2d_mpi_omp";
    fprintf(stderr,
    "Uso: %s [opciones]\n"
    "  --nx <int>        tamaño global X (ej. 1024)\n"
    "  --ny <int>        tamaño global Y (ej. 1024)\n"
    "  --beta <val>      contagio/día (0.6)\n"
    "  --gamma <val>     recuperación/día (0.08)\n"
    "  --D <val>         difusión (0.2)\n"
    "  --dt <val>        paso de tiempo (0.1)\n"
    "  --tmax <val>      tiempo total (60)\n"
    "  --snap_every <n>  snapshot cada n pasos (10; 0=solo t=0)\n"
    "  --out_prefix <p>  prefijo salida (frames/sir2d)\n"
    "  --i0frac <val>    pico infectado inicial (0.01)\n"
    "  --sigma <val>     sigma gauss inicial (3.0)\n"
    "  --threads <int>   hilos OpenMP por proceso\n", p);
}

// mkdir -p
static int mkdir_p(const char* dir) {
    if (!dir || !*dir) return 0;
    char tmp[1024];
    snprintf(tmp, sizeof(tmp), "%s", dir);
    for (int i=(int)strlen(tmp)-1; i>=0 && tmp[i]=='/'; --i) tmp[i]=0;
    if (!*tmp) return 0;
    char* p = tmp;
    if (*p=='/') ++p;
    for (; *p; ++p) {
        if (*p == '/') {
            *p = 0;
            if (mkdir(tmp, 0775) && errno != EEXIST) return -1;
            *p = '/';
        }
    }
    if (mkdir(tmp, 0775) && errno != EEXIST) return -1;
    return 0;
}

static inline int IDX(int y, int x, int lx) { return y*(lx+2) + x; }

// laplaciano 5-puntos con halos ya actualizados
static inline double lap5(const double* A, int x, int y, int lx, int ly){
    return A[IDX(y-1,x,  lx)] + A[IDX(y+1,x,  lx)] +
           A[IDX(y,  x-1,lx)] + A[IDX(y,  x+1,lx)] -
           4.0 * A[IDX(y,x,lx)];
}

// snapshot local (rank): guarda fracción A/N en CSV
static void write_csv_comp_rank(const char* prefix, int step, int rank,
                                const double* S, const double* I, const double* R,
                                int lx, int ly, char tag)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s_%c_t%05d_rank%05d.csv", prefix, tag, step, rank);

    // crear carpeta padre si hay '/'
    const char* slash = strrchr(path, '/');
    if (slash) {
        char dir[1024];
        size_t len = (size_t)(slash - path);
        if (len >= sizeof(dir)) len = sizeof(dir)-1;
        memcpy(dir, path, len); dir[len] = 0;
        (void)mkdir_p(dir);
    }

    FILE* f = fopen(path, "w");
    if (!f) { perror("open snapshot"); return; }

    for (int y=1; y<=ly; ++y) {
        for (int x=1; x<=lx; ++x) {
            int k = IDX(y,x,lx);
            double N = S[k] + I[k] + R[k];
            if (N <= 0) N = 1e-15;
            double A = (tag=='I') ? I[k] : (tag=='R') ? R[k] : S[k];
            double frac = A / N;
            fprintf(f, "%.8g%s", frac, (x==lx) ? "" : ",");
        }
        fputc('\n', f);
    }
    fclose(f);
}

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    Params P = { .NX=1024, .NY=1024, .beta=0.6, .gamma=0.08, .D=0.2,
                 .dt=0.1, .tmax=60.0, .snap_every=10, .out_prefix="frames/sir2d",
                 .i0frac=0.01, .sigma=3.0, .threads=0, .centers=NULL };

    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"--nx") && i+1<argc) P.NX = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--ny") && i+1<argc) P.NY = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--beta") && i+1<argc) P.beta = atof(argv[++i]);
        else if (!strcmp(argv[i],"--gamma") && i+1<argc) P.gamma = atof(argv[++i]);
        else if (!strcmp(argv[i],"--D") && i+1<argc) P.D = atof(argv[++i]);
        else if (!strcmp(argv[i],"--dt") && i+1<argc) P.dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"--tmax") && i+1<argc) P.tmax = atof(argv[++i]);
        else if (!strcmp(argv[i],"--snap_every") && i+1<argc) P.snap_every = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--out_prefix") && i+1<argc) P.out_prefix = argv[++i];
        else if (!strcmp(argv[i],"--i0frac") && i+1<argc) P.i0frac = atof(argv[++i]);
        else if (!strcmp(argv[i],"--sigma") && i+1<argc) P.sigma = atof(argv[++i]);
    else if (!strcmp(argv[i],"--threads") && i+1<argc) P.threads = atoi(argv[++i]);
    else if (!strcmp(argv[i],"--centers") && i+1<argc) P.centers = argv[++i];
        else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) { if(world_rank==0) usage(argv[0]); MPI_Finalize(); return 0; }
        else { if(world_rank==0) usage(argv[0]); MPI_Finalize(); return 1; }
    }

    if (P.threads > 0) omp_set_num_threads(P.threads);

    // 2D malla de procesos con periodicidad (X,Y) = (1,1)
    int dims[2] = {0,0};
    MPI_Dims_create(world_size, 2, dims); // dims[0]=Py (Y), dims[1]=Px (X)
    int periods[2] = {1,1};
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart);

    int rank, coords[2];
    MPI_Comm_rank(cart,&rank);
    MPI_Cart_coords(cart, rank, 2, coords);
    int Py = dims[0], Px = dims[1];

    // vecinos (N/S en Y, W/E en X)
    int nbr_north, nbr_south, nbr_west, nbr_east;
    MPI_Cart_shift(cart, 0, 1, &nbr_north, &nbr_south); // Y
    MPI_Cart_shift(cart, 1, 1, &nbr_west,  &nbr_east ); // X

    // tamaños locales
    if (P.NX % Px || P.NY % Py) {
        if (rank==0) fprintf(stderr,"[Error] nx %% Px != 0 o ny %% Py != 0. Elige tamaños divisibles (%d x %d)\n", Px, Py);
        MPI_Abort(cart, 2);
    }
    int lx = P.NX / Px;  // ancho local
    int ly = P.NY / Py;  // alto local

    // arrays con halos: (ly+2) x (lx+2)
    size_t nloc = (size_t)(lx+2)*(size_t)(ly+2);
    double *S = (double*)malloc(nloc*sizeof(double));
    double *I = (double*)malloc(nloc*sizeof(double));
    double *R = (double*)malloc(nloc*sizeof(double));
    double *Snew = (double*)malloc(nloc*sizeof(double));
    double *Inew = (double*)malloc(nloc*sizeof(double));
    double *Rnew = (double*)malloc(nloc*sizeof(double));
    if (!S||!I||!R||!Snew||!Inew||!Rnew) { fprintf(stderr,"Rank %d: memoria insuficiente\n", rank); MPI_Abort(cart,3); }

    // inicialización gaussiana: uno o varios centros globales
    // Si no se especifica --centers, usar el centro del dominio por defecto.
    // Formato de P.centers: "x,y;x,y;..." en coordenadas globales (0..NX-1, 0..NY-1)
    // Si hay múltiples centros, se reparte i0frac entre ellos (i0frac/num_centers).
    double default_cx = (P.NX-1)/2.0, default_cy = (P.NY-1)/2.0;

    // Parseo simple de la cadena centers
    int ncenters = 0;
    typedef struct { double x; double y; } CPair;
    CPair* CC = NULL;
    if (P.centers && *P.centers) {
        // contar ';' para estimar cantidad
        int count = 1;
        for (const char* s=P.centers; *s; ++s) if (*s==';') ++count;
        CC = (CPair*)malloc(sizeof(CPair)*count);
        if (CC) {
            const char* s = P.centers;
            while (s && *s) {
                double cx=default_cx, cy=default_cy;
                char* endptr;
                cx = strtod(s, &endptr);
                if (!endptr || *endptr!=',') break;
                s = endptr+1;
                cy = strtod(s, &endptr);
                ncenters++;
                CC[ncenters-1].x = cx; CC[ncenters-1].y = cy;
                if (!endptr) break;
                if (*endptr==';') s = endptr+1; else break;
            }
        }
    }
    if (ncenters==0) {
        // usar centro por defecto
        CC = (CPair*)malloc(sizeof(CPair));
        ncenters = 1;
        CC[0].x = default_cx; CC[0].y = default_cy;
    }
    // offset global de este bloque:
    int gx0 = coords[1]*lx; // X global de la columna local 0 (interior x=1)
    int gy0 = coords[0]*ly; // Y global de la fila local 0 (interior y=1)

    #pragma omp parallel for collapse(2) schedule(static)
    for (int y=0; y<ly+2; ++y){
        for (int x=0; x<lx+2; ++x){
            int k = IDX(y,x,lx);
            // por defecto halos 0; se rellenan por MPI
            S[k]=0.0; I[k]=0.0; R[k]=0.0;
        }
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for (int y=1; y<=ly; ++y){
        for (int x=1; x<=lx; ++x){
            int k = IDX(y,x,lx);
            double Xg = gx0 + (x-1);
            double Yg = gy0 + (y-1);
            // suma de múltiples gaussianas
            double i0 = 0.0;
            double frac_each = (ncenters>0) ? (P.i0frac / (double)ncenters) : P.i0frac;
            for (int c=0; c<ncenters; ++c) {
                double dx = Xg - CC[c].x;
                double dy = Yg - CC[c].y;
                double r2 = dx*dx + dy*dy;
                i0 += frac_each * exp(- r2 / (2.0*P.sigma*P.sigma));
            }
            double s0 = 1.0 - i0;
            S[k] = (s0>0)? s0 : 0.0;
            I[k] = (i0>0)? i0 : 0.0;
            R[k] = 0.0;
        }
    }

    long long steps = (long long) llround(P.tmax / P.dt);
    if (rank==0) {
        fprintf(stderr,"[Info] Procs=%d (Py=%d Px=%d)  NX=%d NY=%d  lx=%d ly=%d  beta=%.3g gamma=%.3g D=%.3g dt=%.3g steps=%lld threads=%d\n",
                world_size, Py, Px, P.NX, P.NY, lx, ly, P.beta, P.gamma, P.D, P.dt, steps,
                (P.threads>0?P.threads:omp_get_max_threads()));
        if (P.D>0){
            double dt_cfl = 1.0/(4.0*P.D);
            if (P.dt > dt_cfl) fprintf(stderr,"[Aviso] dt=%.3g > ~1/(4D)=%.3g (posible inestabilidad difusiva)\n", P.dt, dt_cfl);
        }
    }

    // tipos derivados para columnas y filas locales
    MPI_Datatype COL, ROW;
    MPI_Type_vector(ly, 1, (lx+2), MPI_DOUBLE, &COL); MPI_Type_commit(&COL);
    MPI_Type_contiguous(lx, MPI_DOUBLE, &ROW);      MPI_Type_commit(&ROW);

    // snapshot t=0
    if (P.snap_every >= 0) {
        write_csv_comp_rank(P.out_prefix, 0, rank, S,I,R, lx,ly, 'I');
        write_csv_comp_rank(P.out_prefix, 0, rank, S,I,R, lx,ly, 'R');
    }

    MPI_Barrier(cart);
    double t0 = MPI_Wtime();

    for (long long st=1; st<=steps; ++st){

        // === HALO EXCHANGE para S, I, R (4 direcciones) ===
        // Horizontal (x): columnas
        // Enviamos col x=1 a oeste; recibimos ghost derecha x=lx+1 desde este (este nos manda su x=1), etc.
        // S
        MPI_Sendrecv( &S[IDX(1,1,   lx)], 1, COL, nbr_west,  10,
                      &S[IDX(1,lx+1,lx)], 1, COL, nbr_east,  10, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv( &S[IDX(1,lx,  lx)], 1, COL, nbr_east,  11,
                      &S[IDX(1,0,   lx)], 1, COL, nbr_west,  11, cart, MPI_STATUS_IGNORE);
        // I
        MPI_Sendrecv( &I[IDX(1,1,   lx)], 1, COL, nbr_west,  20,
                      &I[IDX(1,lx+1,lx)], 1, COL, nbr_east,  20, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv( &I[IDX(1,lx,  lx)], 1, COL, nbr_east,  21,
                      &I[IDX(1,0,   lx)], 1, COL, nbr_west,  21, cart, MPI_STATUS_IGNORE);
        // R
        MPI_Sendrecv( &R[IDX(1,1,   lx)], 1, COL, nbr_west,  30,
                      &R[IDX(1,lx+1,lx)], 1, COL, nbr_east,  30, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv( &R[IDX(1,lx,  lx)], 1, COL, nbr_east,  31,
                      &R[IDX(1,0,   lx)], 1, COL, nbr_west,  31, cart, MPI_STATUS_IGNORE);

        // Vertical (y): filas
        // S
        MPI_Sendrecv( &S[IDX(1,    1,lx)], 1, ROW, nbr_north, 12,
                      &S[IDX(ly+1, 1,lx)], 1, ROW, nbr_south, 12, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv( &S[IDX(ly,   1,lx)], 1, ROW, nbr_south, 13,
                      &S[IDX(0,    1,lx)], 1, ROW, nbr_north, 13, cart, MPI_STATUS_IGNORE);
        // I
        MPI_Sendrecv( &I[IDX(1,    1,lx)], 1, ROW, nbr_north, 22,
                      &I[IDX(ly+1, 1,lx)], 1, ROW, nbr_south, 22, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv( &I[IDX(ly,   1,lx)], 1, ROW, nbr_south, 23,
                      &I[IDX(0,    1,lx)], 1, ROW, nbr_north, 23, cart, MPI_STATUS_IGNORE);
        // R
        MPI_Sendrecv( &R[IDX(1,    1,lx)], 1, ROW, nbr_north, 32,
                      &R[IDX(ly+1, 1,lx)], 1, ROW, nbr_south, 32, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv( &R[IDX(ly,   1,lx)], 1, ROW, nbr_south, 33,
                      &R[IDX(0,    1,lx)], 1, ROW, nbr_north, 33, cart, MPI_STATUS_IGNORE);

        // === Actualización explícita (OpenMP) en todo el dominio interior (1..ly, 1..lx) ===
        #pragma omp parallel for collapse(2) schedule(static)
        for (int y=1; y<=ly; ++y){
            for (int x=1; x<=lx; ++x){
                int k = IDX(y,x,lx);
                double s = S[k], i = I[k], r = R[k];
                double Nloc = s + i + r; if (Nloc <= 0) Nloc = 1e-15;

                double inf = P.beta * s * i / Nloc;
                double dS = -inf + P.D * lap5(S,x,y,lx,ly);
                double dI =  inf - P.gamma*i + P.D * lap5(I,x,y,lx,ly);
                double dR =  P.gamma*i + P.D * lap5(R,x,y,lx,ly);

                double sn = s + P.dt*dS;
                double in = i + P.dt*dI;
                double rn = r + P.dt*dR;
                if (sn < 0) sn = 0;
                if (in < 0) in = 0;
                if (rn < 0) rn = 0;

                Snew[k]=sn; Inew[k]=in; Rnew[k]=rn;
            }
        }

        // swap
        double* tmp;
        tmp=S; S=Snew; Snew=tmp;
        tmp=I; I=Inew; Inew=tmp;
        tmp=R; R=Rnew; Rnew=tmp;

        if (P.snap_every>0 && (st % P.snap_every)==0){
            write_csv_comp_rank(P.out_prefix, (int)st, rank, S,I,R, lx,ly, 'I');
            write_csv_comp_rank(P.out_prefix, (int)st, rank, S,I,R, lx,ly, 'R');
        }
    }

    double t1 = MPI_Wtime();
    double local_updates = (double)lx * (double)ly * (double)steps;
    double local_time = t1 - t0;
    double global_time, global_updates;

    MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart);
    MPI_Reduce(&local_updates, &global_updates, 1, MPI_DOUBLE, MPI_SUM, 0, cart);

    // conservación promedio local
    double local_sum = 0.0;
    for (int y=1;y<=ly;++y)
        for (int x=1;x<=lx;++x){
            int k=IDX(y,x,lx);
            local_sum += S[k]+I[k]+R[k];
        }
    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, cart);

    if (rank==0){
        double Ntot = (double)P.NX*(double)P.NY;
        double avg = global_sum / Ntot;
        double MLUPS = (global_updates / global_time) / 1e6;
        fprintf(stderr,"[Resumen] Conservación promedio por celda ~ %.6g (≈1)\n", avg);
        fprintf(stderr,"[Tiempo] %.3f s  |  Throughput: %.2f MLUPS (global)\n", global_time, MLUPS);
    }

    MPI_Type_free(&COL); MPI_Type_free(&ROW);
    free(S); free(I); free(R); free(Snew); free(Inew); free(Rnew);
    if (CC) free(CC);

    MPI_Comm_free(&cart);
    MPI_Finalize();
    return 0;
}
