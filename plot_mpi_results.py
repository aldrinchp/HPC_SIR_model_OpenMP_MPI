#!/usr/bin/env python3
"""
plot_mpi_results.py

Grafica resultados de experimentos MPI:
  - Eje X: Número de procesos MPI
  - Eje Y: Tiempo de ejecución (segundos)
  - Líneas: Diferentes números de threads por proceso

Uso: python3 plot_mpi_results.py [results_mpi/mpi_results.csv]
"""
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fn = sys.argv[1] if len(sys.argv) > 1 else 'results_mpi/mpi_results.csv'
df = pd.read_csv(fn, comment='#', header=None, names=['procs','threads','time_s','mlups'])

outdir = os.path.join(os.path.dirname(fn), 'plots')
os.makedirs(outdir, exist_ok=True)

# Cleanup: ensure numeric, remove failed runs (time=0)
df['procs'] = df['procs'].astype(int)
df['threads'] = df['threads'].astype(int)
df['time_s'] = df['time_s'].astype(float)
df = df[df['time_s'] > 0]  # Filtrar ejecuciones fallidas

if df.empty:
    print("No hay datos válidos para graficar. Verifica que las ejecuciones fueron exitosas.")
    sys.exit(1)

# Gráfica: tiempo vs procesos, líneas por threads
plt.figure(figsize=(10,7))
for t in sorted(df['threads'].unique(), reverse=True):  # Orden descendente para leyenda
    sub = df[df['threads']==t].sort_values('procs')
    if not sub.empty:
        plt.plot(sub['procs'], sub['time_s'], marker='o', label=f'{t} threads/proc', linewidth=2)

plt.xscale('log', base=2)
plt.xlabel('Número de procesos MPI', fontsize=12)
plt.ylabel('Tiempo de ejecución (s)', fontsize=12)
plt.title('Tiempo vs Procesos MPI (N=1024 fijo)', fontsize=14, fontweight='bold')
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.legend(fontsize=10)
plt.xticks(sorted(df['procs'].unique()), sorted(df['procs'].unique()))  # Etiquetas explícitas

f1 = os.path.join(outdir, 'time_vs_procs.png')
plt.savefig(f1, dpi=200, bbox_inches='tight')
print('Saved', f1)
plt.close()

# Gráfica adicional: Speedup vs procesos (respecto a 1 proc con mismo thread count)
plt.figure(figsize=(10,7))
for t in sorted(df['threads'].unique(), reverse=True):
    sub = df[df['threads']==t].sort_values('procs')
    if not sub.empty:
        time_1proc = sub[sub['procs']==1]['time_s'].values
        if time_1proc.size > 0:
            time_1 = time_1proc[0]
            sub_copy = sub.copy()
            sub_copy['speedup'] = time_1 / sub_copy['time_s']
            plt.plot(sub_copy['procs'], sub_copy['speedup'], marker='o', label=f'{t} threads/proc', linewidth=2)

plt.xscale('log', base=2)
plt.xlabel('Número de procesos MPI', fontsize=12)
plt.ylabel('Speedup (respecto a 1 proceso)', fontsize=12)
plt.title('Speedup vs Procesos MPI (N=1024 fijo)', fontsize=14, fontweight='bold')
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.legend(fontsize=10)
# Línea ideal (speedup lineal)
procs_range = sorted(df['procs'].unique())
plt.plot(procs_range, procs_range, 'k--', alpha=0.3, label='Speedup ideal (lineal)')
plt.xticks(procs_range, procs_range)

f2 = os.path.join(outdir, 'speedup_vs_procs.png')
plt.savefig(f2, dpi=200, bbox_inches='tight')
print('Saved', f2)
plt.close()

print('Gráficas generadas en', outdir)
