Experimentos y gráficas para `sir2d_mpi_omp`

Este documento describe cómo generar datos de rendimiento del código híbrido MPI+OpenMP, cómo graficar los resultados y recomendaciones básicas para ejecutar los experimentos.

## Requisitos

- Linux o similar
- GNU make / compilador (gcc/icc) para compilar el binario (los scripts compilan si hace falta)
- Open MPI o implementación MPI compatible (mpirun/mpiexec)
- Python 3 y las dependencias listadas en `requirements.txt`

## Estructura relevante del repositorio

- `run_experiments.sh` : script para barridos secuenciales / OpenMP (genera `results/results.csv`)
- `run_mpi_experiments.sh` : script para experimentos híbridos (genera `results_mpi/mpi_results.csv`)
- `plot_results.py` : genera gráficas desde `results/results.csv`
- `plot_mpi_results.py` : genera gráficas desde `results_mpi/mpi_results.csv`
- `results/`, `results_mpi/` : salidas de los experimentos

## 1) Generar datos (OpenMP / barridos)

- Edita `run_experiments.sh` para ajustar las listas `NS` y `THREADS` acorde a tu máquina.
- Ejecuta:

```bash
./run_experiments.sh
```

El script compilará el binario si no existe y ejecutará los barridos, guardando `results/results.csv` y logs en `results/`.

Notas importantes:

- Para evitar I/O intensivo, por defecto `run_experiments.sh` ejecuta con `--snap_every 0`.
- Ajusta `TMAX`, `DT` y `NS` para hacer pruebas más rápidas o más representativas según tus recursos.

## 2) Graficar

Instala dependencias (recomendado en un entorno virtual):

```bash
python3 -m pip install -r requirements.txt
```

Genera las gráficas:

```bash
python3 plot_results.py results/results.csv
```

Las gráficas resultantes se guardarán en `results/plots/`:
- `time_vs_N.png` : tiempo vs N (líneas por número de hilos)
- `speedup_N<...>.png` : speedup vs hilos para el N más grande medido

## 3) Experimentos MPI (N fijo, variando procesos y hilos)

Script: `run_mpi_experiments.sh`

Parámetros por defecto que suelen aparecer en el script (revisa el script para confirmar):

- `N_FIXED=8192` (puedes cambiar a 4096/2048 para mediciones más rápidas)
- `PROCS=(1 2 4 8 16 32 64)` : número de procesos MPI a probar
- `THREADS=(64 32 16 4 2 1)` : hilos OpenMP por proceso
- `TMAX=5.0`, `DT=0.1`

### Requisito de divisibilidad del dominio

- El código exige `NX % Px == 0` y `NY % Py == 0`, donde `Px*Py = P` (número de procesos) y `Px,Py` se obtienen mediante `MPI_Dims_create`.
- Con `N` potencia de 2 (p. ej. 8192), `P` potencias de 2 suelen dar factorizaciones válidas (rectangular/cuadrada).

### Ejecución (ejemplos)

Sin oversubscribe (recomendado cuando la máquina tiene suficientes cores físicos):

```bash
# Ajusta PROCS en el script si haces pruebas locales
# sed -i "s/^PROCS=.*/PROCS=(1 2 4 8 16 32)/" run_mpi_experiments.sh
./run_mpi_experiments.sh
```

Con oversubscribe (cuando hay pocos recursos físicos):

```bash
# El script puede incluir --oversubscribe en mpirun; revisa el script antes de usarlo
./run_mpi_experiments.sh
```

### Afinidad y binding (recomendado para resultados reproducibles)

- OpenMP: exportar `OMP_PROC_BIND=true` y `OMP_PLACES=cores`
- Open MPI: usar opciones como `--map-by ppr:1:core --bind-to core` para fijar procesos a cores

## 4) Generar gráficas de MPI

```bash
python3 plot_mpi_results.py results_mpi/mpi_results.csv
```

Salidas típicas (en `results_mpi/plots/`):
- `time_vs_procs.png` : tiempo vs procesos (líneas por threads/proc)
- `speedup_vs_procs.png` : speedup respecto a 1 proceso (por cada threads/proc)

## Notas finales

- Revisa los scripts `run_experiments.sh` y `run_mpi_experiments.sh` antes de ejecutar para ajustar parámetros por defecto a tu máquina.
- Si detectas discrepancias entre lo descrito en este README y el contenido de los scripts, indícalo y puedo actualizar el README o los scripts según convenga.

