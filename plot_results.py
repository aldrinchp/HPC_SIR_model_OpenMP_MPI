#!/usr/bin/env python3
"""
plot_results.py

Lectura de CSV generado por run_experiments.sh y generación de dos gráficas:
  - tiempo vs N (líneas por número de hilos)
  - speedup (tiempo_1thread / tiempo_T) vs threads para un N fijo

Uso: python3 plot_results.py [results.csv]
"""
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fn = sys.argv[1] if len(sys.argv) > 1 else 'results/results.csv'
df = pd.read_csv(fn, comment='#', header=None, names=['N','threads','time_s','mlups'])

outdir = os.path.join(os.path.dirname(fn), 'plots')
os.makedirs(outdir, exist_ok=True)

# Cleanup: ensure numeric
df['N'] = df['N'].astype(int)
df['threads'] = df['threads'].astype(int)
df['time_s'] = df['time_s'].astype(float)

# 1) tiempo vs N para cada threads
plt.figure(figsize=(8,6))
for t in sorted(df['threads'].unique()):
    sub = df[df['threads']==t].sort_values('N')
    plt.plot(sub['N'], sub['time_s'], marker='o', label=f'{t} threads')

plt.xscale('log', base=2)
plt.yscale('log')
plt.xlabel('N (NX=NY)')
plt.ylabel('Tiempo (s)')
plt.title('Tiempo vs N por número de hilos')
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.legend()
f1 = os.path.join(outdir, 'time_vs_N.png')
plt.savefig(f1, dpi=200, bbox_inches='tight')
print('Saved', f1)
plt.close()

# 2) speedup vs threads para un N fijo (elegir mayor N disponible si no especificado)
# Buscar N fijo: usar el máximo N experimental (más interesante)
N_choices = sorted(df['N'].unique())
N_fixed = N_choices[-1]
subN = df[df['N']==N_fixed].sort_values('threads')
if subN.empty:
    print('No hay datos para N fijo. Abortando speedup plot.')
    sys.exit(1)

time_1 = subN[subN['threads']==1]['time_s'].values
if time_1.size == 0:
    # intentar primer hilo disponible
    time_1 = subN['time_s'].values[0]
else:
    time_1 = time_1[0]

subN = subN.copy()
subN['speedup'] = time_1 / subN['time_s']

plt.figure(figsize=(8,6))
plt.plot(subN['threads'], subN['speedup'], marker='o')
plt.xticks(subN['threads'])
plt.xlabel('Threads')
plt.ylabel('Speedup (time_1thread / time_T)')
plt.title(f'Speedup vs threads (N={N_fixed})')
plt.grid(True, ls='--', alpha=0.4)
f2 = os.path.join(outdir, f'speedup_N{N_fixed}.png')
plt.savefig(f2, dpi=200, bbox_inches='tight')
print('Saved', f2)
plt.close()

print('Plots generated in', outdir)
