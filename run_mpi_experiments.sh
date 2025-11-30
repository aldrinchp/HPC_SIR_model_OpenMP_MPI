#!/usr/bin/env bash
set -euo pipefail

# run_mpi_experiments.sh
# Experimento: N fijo (1024), variar procesos MPI y threads OpenMP
# Eje X: Procesos MPI
# Eje Y: Tiempo
# Líneas: diferentes números de threads por proceso

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/sir2d_mpi_omp"
OUTDIR="$HERE/results_mpi"
mkdir -p "$OUTDIR"

# N fijo (máximo del experimento anterior)
N_FIXED=8192

# Arrays de procesos y threads
# NOTA: Usando --oversubscribe para permitir más procesos que cores físicos
# P debe ser potencia de 4 para mallas 2D cuadradas: 1, 4, 16, 64
PROCS=(1 2 4 8 16 32 64)
THREADS=(64 32 16 4 2 1)

# Duración del experimento
TMAX=5.0
DT=0.1

RESULT_CSV="$OUTDIR/mpi_results.csv"
echo "# procs,threads,time_s,mlups" > "$RESULT_CSV"

if [ ! -x "$BIN" ]; then
  echo "Binario $BIN no encontrado o no ejecutable. Intentando compilar con mpicc..."
  if command -v mpicc >/dev/null 2>&1; then
    (cd "$HERE" && mpicc -O3 -march=native -fopenmp -o sir2d_mpi_omp sir2d_mpi_omp.c -lm)
  else
    echo "mpicc no está disponible en PATH. Instala MPI o compila manualmente." >&2
    exit 1
  fi
fi

echo "Ejecutando barrido MPI: N=$N_FIXED PROCS=${PROCS[*]} THREADS=${THREADS[*]} TMAX=${TMAX}"
for P in "${PROCS[@]}"; do
  for T in "${THREADS[@]}"; do
    echo "--> procs=$P  threads=$T"
    # Ejecutar con P procesos MPI y T threads por proceso
    STDERR="$OUTDIR/run_P${P}_T${T}.log"
    
    # Usamos --oversubscribe para permitir más procesos que cores físicos
    mpirun --oversubscribe -np "$P" "$BIN" --nx "$N_FIXED" --ny "$N_FIXED" --threads "$T" --tmax "$TMAX" --dt "$DT" --snap_every 0 2>"$STDERR" || {
      echo "   -> ERROR: P=$P no es compatible con N=$N_FIXED (ver log)"
      echo "$P,$T,0,0" >> "$RESULT_CSV"
      continue
    }

    # Extraer tiempo y MLUPS desde stderr
    TIME_S=$(grep "\[Tiempo\]" "$STDERR" | sed -E 's/.*\[Tiempo\] ([0-9.]+) s.*/\1/' || echo "0")
    MLUPS=$(grep "\[Tiempo\]" "$STDERR" | sed -E 's/.*Throughput: *([0-9.]+) MLUPS.*/\1/' || echo "0")

    echo "$P,$T,${TIME_S},${MLUPS}" >> "$RESULT_CSV"
    echo "   -> time=${TIME_S}s  mlups=${MLUPS}"
  done
done

echo "Experimentos MPI completados. Resultados en $RESULT_CSV"
echo "Graficar con: /media/aldrinchp/KAli-HDD/Subjects/8vo_AA_MA/.venv/bin/python plot_mpi_results.py $RESULT_CSV"
