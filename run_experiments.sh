#!/usr/bin/env bash
set -euo pipefail

# run_experiments.sh
# Ejecuta barridos de tamaño N (NX=NY=N) y de número de hilos, guardando tiempos en results/results.csv
# Ajusta los vectores NS y THREADS según tu máquina.

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/sir2d_mpi_omp"
OUTDIR="$HERE/results"
mkdir -p "$OUTDIR"

# Default arrays (modifica según conveniencia). N debe ser divisible por Px y Py (aquí usamos 1 proc).
NS=(1024 2048 4096 8192)
THREADS=(1 2 4 8)

# Duración del experimento (acortar para pruebas rápidas)
TMAX=5.0
DT=0.1

RESULT_CSV="$OUTDIR/results.csv"
echo "# N,threads,time_s,mlups" > "$RESULT_CSV"

if [ ! -x "$BIN" ]; then
  echo "Binario $BIN no encontrado o no ejecutable. Intentando compilar con mpicc..."
  if command -v mpicc >/dev/null 2>&1; then
    (cd "$HERE" && mpicc -O3 -march=native -fopenmp -o sir2d_mpi_omp sir2d_mpi_omp.c -lm)
  else
    echo "mpicc no está disponible en PATH. Instala MPI o compila manualmente." >&2
    exit 1
  fi
fi

echo "Ejecutando barrido: NS=${NS[*]} THREADS=${THREADS[*]} TMAX=${TMAX}"
for N in "${NS[@]}"; do
  for T in "${THREADS[@]}"; do
    echo "--> N=$N  threads=$T"
    # Ejecutar con 1 proceso (mpirun -np 1) y pasar --threads
    # Deshabilitamos snapshots para evitar I/O pesado usando --snap_every 0
    STDERR="$OUTDIR/run_N${N}_T${T}.log"
    mpirun -np 1 "$BIN" --nx "$N" --ny "$N" --threads "$T" --tmax "$TMAX" --dt "$DT" --snap_every 0 2>"$STDERR" || true

    # Extraer tiempo y MLUPS desde stderr
    TIME_S=$(grep "\[Tiempo\]" "$STDERR" | sed -E 's/.*\[Tiempo\] ([0-9.]+) s.*/\1/' || true)
    MLUPS=$(grep "\[Tiempo\]" "$STDERR" | sed -E 's/.*Throughput: *([0-9.]+) MLUPS.*/\1/' || true)
    # Si no se encontró la línea anterior, intentar con otras líneas
    if [ -z "$TIME_S" ]; then
      TIME_S=$(grep "\[Resumen\]" -n "$STDERR" || true)
    fi
    if [ -z "$MLUPS" ]; then MLUPS=0; fi

    echo "$N,$T,${TIME_S:-0},${MLUPS}" >> "$RESULT_CSV"
    echo "   -> time=${TIME_S:-0}s  mlups=${MLUPS}"
  done
done

echo "Experimentos completados. Resultados en $RESULT_CSV"
echo "Graficar con: python3 plot_results.py $RESULT_CSV"
