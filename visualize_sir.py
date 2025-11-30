import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.animation import FuncAnimation, FFMpegWriter

def load_snapshot(prefix, timestep, component='I'):
    """Carga todos los archivos CSV de un timestep y los ensambla"""
    pattern = f"{prefix}_{component}_t{timestep:05d}_rank*.csv"
    files = sorted(glob.glob(pattern))
    
    if not files:
        return None
    
    data_pieces = []
    for f in files:
        data = np.loadtxt(f, delimiter=',')
        data_pieces.append(data)
    
    n_ranks = len(files)
    if n_ranks == 4:
        grid = np.block([
            [data_pieces[0], data_pieces[1]],
            [data_pieces[2], data_pieces[3]]
        ])
    elif n_ranks == 1:
        grid = data_pieces[0]
    else:
        grid = np.concatenate([np.concatenate(data_pieces, axis=1)])
    
    return grid

def create_animation_full(prefix, tmax, snap_every, dt, output='sir_full.mp4'):
    """Animación con S+I+R"""
    timesteps = list(range(0, int(tmax/dt)+1, snap_every))
    
    print(f"📊 Generando {len(timesteps)} frames (S+I+R)...")
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
    
    I0 = load_snapshot(prefix, 0, 'I')
    R0 = load_snapshot(prefix, 0, 'R')
    
    if I0 is None:
        print("❌ No se encontraron archivos.")
        return
    
    S0 = 1.0 - I0 - R0
    
    # Susceptibles
    im1 = ax1.imshow(S0, cmap='Greens', vmin=0, vmax=1, origin='lower')
    ax1.set_title('Susceptibles (S)', fontsize=14)
    plt.colorbar(im1, ax=ax1, label='S/N')
    
    # Infectados
    im2 = ax2.imshow(I0, cmap='Reds', vmin=0, vmax=0.15, origin='lower')
    ax2.set_title('Infectados (I)', fontsize=14)
    plt.colorbar(im2, ax=ax2, label='I/N')
    
    # Recuperados
    im3 = ax3.imshow(R0, cmap='Blues', vmin=0, vmax=1, origin='lower')
    ax3.set_title('Recuperados (R)', fontsize=14)
    plt.colorbar(im3, ax=ax3, label='R/N')
    
    time_text = fig.text(0.5, 0.96, '', ha='center', fontsize=14, weight='bold')
    
    def update(frame):
        t = timesteps[frame]
        I = load_snapshot(prefix, t, 'I')
        R = load_snapshot(prefix, t, 'R')
        
        if I is not None and R is not None:
            S = 1.0 - I - R
            im1.set_data(S)
            im2.set_data(I)
            im3.set_data(R)
            time_text.set_text(f'Tiempo: {t*dt:.1f} días (paso {t})')
        
        if frame % 3 == 0:
            print(f"  Frame {frame+1}/{len(timesteps)}...")
        
        return im1, im2, im3, time_text
    
    anim = FuncAnimation(fig, update, frames=len(timesteps), 
                         interval=200, blit=False, repeat=False)
    
    print(f"💾 Guardando...")
    writer = FFMpegWriter(fps=5, bitrate=3000)
    anim.save(output, writer=writer)
    plt.close(fig)
    
    print(f"✅ Video guardado: {output}")

def plot_single_frame(prefix, timestep, dt):
    """Frame individual S+I+R"""
    I = load_snapshot(prefix, timestep, 'I')
    R = load_snapshot(prefix, timestep, 'R')
    
    if I is None:
        print(f"❌ No encontrado t={timestep}")
        return
    
    S = 1.0 - I - R
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    im1 = axes[0].imshow(S, cmap='Greens', vmin=0, vmax=1, origin='lower')
    axes[0].set_title(f'Susceptibles - t={timestep*dt:.1f} días', fontsize=14)
    plt.colorbar(im1, ax=axes[0])
    
    im2 = axes[1].imshow(I, cmap='Reds', vmin=0, vmax=0.15, origin='lower')
    axes[1].set_title(f'Infectados - t={timestep*dt:.1f} días', fontsize=14)
    plt.colorbar(im2, ax=axes[1])
    
    im3 = axes[2].imshow(R, cmap='Blues', vmin=0, vmax=1, origin='lower')
    axes[2].set_title(f'Recuperados - t={timestep*dt:.1f} días', fontsize=14)
    plt.colorbar(im3, ax=axes[2])
    
    plt.tight_layout()
    filename = f'sir_frame_t{timestep:05d}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    print(f"✅ Guardado: {filename}")

if __name__ == "__main__":
    import sys
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    
    # Permitir especificar prefijo via argumento --prefix
    PREFIX = os.path.join(SCRIPT_DIR, "frames/sir2d")
    if len(sys.argv) > 1:
        for i, arg in enumerate(sys.argv):
            if arg == "--prefix" and i+1 < len(sys.argv):
                PREFIX = os.path.join(SCRIPT_DIR, "frames", sys.argv[i+1])
    
    TMAX = 200
    SNAP_EVERY = 50
    DT = 0.1
    
    print("🎬 Visualizador SIR 2D Completo (S+I+R)")
    print(f"📁 Buscando en: {PREFIX}\n")
    
    test_files = glob.glob(f"{PREFIX}_I_t00000_rank*.csv")
    if not test_files:
        print(f"❌ Sin archivos")
        exit(1)
    print(f"✅ {len(test_files)} ranks encontrados\n")
    
    print("1. Ver frame específico")
    print("2. Crear animación completa (S+I+R)")
    choice = input("Elige (1/2): ")
    
    if choice == "1":
        t = int(input(f"Timestep (0-{int(TMAX/DT)}, múltiplo de {SNAP_EVERY}): "))
        plot_single_frame(PREFIX, t, DT)
    else:
        create_animation_full(PREFIX, TMAX, SNAP_EVERY, DT,
                             output=os.path.join(SCRIPT_DIR, 'sir_complete.mp4'))