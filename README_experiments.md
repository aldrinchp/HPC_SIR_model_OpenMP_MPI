Experimentos y gráficas para `sir2d_mpi_omp`

1) Generar datos

- Edita `run_experiments.sh` para ajustar las listas `NS` y `THREADS` acorde a tu máquina.
- Ejecuta (desde la carpeta `Hibrido`):

```bash
./run_experiments.sh
```

Esto compilará el binario si no existe y ejecutará los barridos, guardando `results/results.csv` y logs en `results/`.

2) Graficar

Instala dependencias:

```bash
python3 -m pip install -r requirements.txt
```

Luego:

```bash
python3 plot_results.py results/results.csv
```

Las gráficas se guardarán en `results/plots/`:
- `time_vs_N.png` : tiempo vs N (líneas por número de hilos)
- `speedup_N<...>.png` : speedup vs threads para el N más grande medido

Notas:
- Para evitar I/O intensivo, por defecto `run_experiments.sh` ejecuta con `--snap_every 0`.
- Ajusta `TMAX` y `NS` para hacer pruebas más rápidas o más representativas.

2) Experimentos MPI (N fijo, variando procesos y hilos)

Script: `run_mpi_experiments.sh`

- Parámetros por defecto en el script:
	- `N_FIXED=8192` (puedes cambiarlo a 4096/2048 si quieres medir más rápido)
	- `PROCS=(1 2 4 8 16 32 64)` número de procesos MPI a probar
	- `THREADS=(64 32 16 4 2 1)` hilos OpenMP por proceso
	- `TMAX=5.0`, `DT=0.1`

Requisito de divisibilidad del dominio:
- El código exige `NX % Px == 0` y `NY % Py == 0`, donde `Px*Py = P` (número de procesos) y `dims = MPI_Dims_create(P,2, ...)`.
- Con `N` potencia de 2 (p.ej. 8192), cualquier `P` potencia de 2 suele tener descomposición válida (rectangular o cuadrada).
- Mallas cuadradas (p.ej. `P=4,16,64`) tienden a equilibrar mejor la comunicación.

Ejecución en máquina con suficientes recursos (sin oversubscribe):

```bash
# Recomendado: ajustar PROCS para no exceder los cores físicos
# Ejemplo en máquina con 32 cores
sed -i 's/^PROCS=.*/PROCS=(1 2 4 8 16 32)/' run_mpi_experiments.sh

# Ejecutar experimento MPI sin oversubscribe
./run_mpi_experiments.sh
```

Ejecución en máquina con pocos recursos (con oversubscribe):

```bash
# Mantener PROCS con más valores para construir una curva más densa
# El script ya incluye --oversubscribe en mpirun
./run_mpi_experiments.sh
```

Sugerencias de afinidad (opcional, cuando hay recursos suficientes):
- OpenMP: exportar `OMP_PROC_BIND=true` y `OMP_PLACES=cores`
- Open MPI: usar `--map-by ppr:1:core --bind-to core` para fijar cada proceso a cores distintos.

Memoria aproximada:
- N=8192 con P=1 puede usar ~3.2 GB por proceso (6 matrices double con halos).
- Con P=4, cada proceso usa ~0.8 GB aprox.

Generar gráficas de MPI:

```bash
# Usar el venv sugerido
/media/aldrinchp/KAli-HDD/Subjects/8vo_AA_MA/.venv/bin/python plot_mpi_results.py results_mpi/mpi_results.csv
```

Salida:
- `results_mpi/plots/time_vs_procs.png` : tiempo vs procesos (líneas por threads/proc)
- `results_mpi/plots/speedup_vs_procs.png` : speedup respecto a 1 proceso (por cada threads/proc)
