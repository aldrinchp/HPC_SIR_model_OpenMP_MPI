#!/usr/bin/env bash
set -euo pipefail

# run_simulation_video.sh
# Script para ejecutar simulación SIR y generar video automáticamente

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/sir2d_mpi_omp"
VENV_PYTHON="/media/aldrinchp/KAli-HDD/Subjects/8vo_AA_MA/.venv/bin/python"

echo "🎬 Generador de Simulación y Video SIR 2D"
echo "=========================================="
echo ""

# Verificar que el binario existe
if [ ! -x "$BIN" ]; then
  echo "⚠️  Binario no encontrado. Compilando..."
  mpicc -O3 -march=native -fopenmp -o sir2d_mpi_omp sir2d_mpi_omp.c -lm
  echo "✅ Compilación exitosa"
fi

# =======================
# CONFIGURACIÓN RÁPIDA
# =======================
echo "Selecciona un preset o configura manualmente:"
echo ""
echo "Presets disponibles:"
echo "  1) Demo rápida     - N=512,  t=100, 1 foco  (~30s)"
echo "  2) Demo estándar   - N=1024, t=200, 1 foco  (~2min)"
echo "  3) Tres focos      - N=1024, t=200, 3 focos (~2min)"
echo "  4) Alta resolución - N=2048, t=200, 1 foco  (~10min)"
echo "  5) Configuración manual"
echo ""
read -p "Elige (1-5): " PRESET

case $PRESET in
  1)
    NX=512; NY=512; TMAX=100; SNAP_EVERY=10; BETA=0.6; GAMMA=0.08
    D=0.2; I0FRAC=0.01; SIGMA=3; CENTERS=""; THREADS=8; NAME="demo_rapida"
    ;;
  2)
    NX=1024; NY=1024; TMAX=200; SNAP_EVERY=25; BETA=0.6; GAMMA=0.08
    D=0.2; I0FRAC=0.01; SIGMA=3; CENTERS=""; THREADS=8; NAME="demo_estandar"
    ;;
  3)
    NX=1024; NY=1024; TMAX=200; SNAP_EVERY=25; BETA=0.6; GAMMA=0.08
    D=0.2; I0FRAC=0.03; SIGMA=6; CENTERS="256,256;768,768;512,512"
    THREADS=8; NAME="tres_focos"
    ;;
  4)
    NX=2048; NY=2048; TMAX=200; SNAP_EVERY=25; BETA=0.6; GAMMA=0.08
    D=0.2; I0FRAC=0.01; SIGMA=4; CENTERS=""; THREADS=8; NAME="alta_resolucion"
    ;;
  5)
    echo ""
    echo "📝 Configuración manual:"
    read -p "  Tamaño de malla (NX=NY, potencia de 2, ej: 1024): " NX
    NY=$NX
    read -p "  Tiempo total (días, ej: 200): " TMAX
    read -p "  Snapshots cada N pasos (ej: 25 para ~80 frames): " SNAP_EVERY
    read -p "  Beta - tasa de contagio (ej: 0.6): " BETA
    read -p "  Gamma - tasa de recuperación (ej: 0.08): " GAMMA
    read -p "  D - coeficiente de difusión (ej: 0.2): " D
    read -p "  Fracción inicial infectada (ej: 0.01): " I0FRAC
    read -p "  Sigma gaussiana inicial (ej: 3): " SIGMA
    read -p "  Centros múltiples 'x,y;x,y;...' (dejar vacío para centro): " CENTERS
    read -p "  Número de hilos OpenMP (ej: 8): " THREADS
    read -p "  Nombre del experimento (sin espacios): " NAME
    ;;
  *)
    echo "❌ Opción inválida"
    exit 1
    ;;
esac

DT=0.1
PREFIX="frames/${NAME}"
OUTDIR="$HERE/frames"
VIDEODIR="$HERE/results/plots"
mkdir -p "$OUTDIR" "$VIDEODIR"

# Mostrar resumen
echo ""
echo "📋 Configuración seleccionada:"
echo "   Malla: ${NX}×${NY}"
echo "   Tiempo: 0 → ${TMAX} días (dt=${DT})"
echo "   Snapshots cada: ${SNAP_EVERY} pasos"
echo "   Beta: ${BETA}, Gamma: ${GAMMA}, D: ${D}"
echo "   i0frac: ${I0FRAC}, sigma: ${SIGMA}"
echo "   Centros: ${CENTERS:-"centro por defecto"}"
echo "   Hilos: ${THREADS}"
echo "   Prefijo: ${NAME}"
echo ""

# Estimar tiempo y frames
STEPS=$(echo "$TMAX / $DT" | bc)
FRAMES=$(echo "$STEPS / $SNAP_EVERY + 1" | bc)
echo "📊 Estimaciones:"
echo "   Pasos totales: ${STEPS}"
echo "   Frames en video: ~${FRAMES}"
echo ""

read -p "¿Continuar? (s/n): " CONFIRM
if [[ ! "$CONFIRM" =~ ^[sS]$ ]]; then
  echo "❌ Cancelado"
  exit 0
fi

# Limpiar frames antiguos del mismo experimento
echo ""
echo "🧹 Limpiando frames antiguos de '${NAME}'..."
rm -f "${OUTDIR}/${NAME}"*.csv 2>/dev/null || true

# =======================
# EJECUTAR SIMULACIÓN
# =======================
echo ""
echo "🚀 Ejecutando simulación..."
echo ""

CMD="mpirun -np 1 \"$BIN\" --nx $NX --ny $NY --beta $BETA --gamma $GAMMA --D $D --dt $DT --tmax $TMAX --snap_every $SNAP_EVERY --out_prefix \"$PREFIX\" --i0frac $I0FRAC --sigma $SIGMA --threads $THREADS"

if [ -n "$CENTERS" ]; then
  CMD="$CMD --centers \"$CENTERS\""
fi

echo "Comando: $CMD"
echo ""

eval $CMD

if [ $? -ne 0 ]; then
  echo "❌ Error en la simulación"
  exit 1
fi

echo ""
echo "✅ Simulación completada"
echo ""

# =======================
# GENERAR VIDEO
# =======================
echo "🎥 Generando video..."
echo ""

# Verificar que hay frames
FRAME_COUNT=$(ls -1 "${OUTDIR}/${NAME}"_I_t*.csv 2>/dev/null | wc -l)
if [ "$FRAME_COUNT" -eq 0 ]; then
  echo "❌ No se encontraron frames. Verifica la simulación."
  exit 1
fi

echo "   Frames encontrados: ${FRAME_COUNT}"
echo ""

# Ejecutar visualizador automáticamente (opción 2)
echo "2" | "$VENV_PYTHON" "$HERE/visualize_sir.py" --prefix "$NAME"

if [ $? -ne 0 ]; then
  echo "❌ Error generando video"
  exit 1
fi

# Mover y renombrar video
VIDEO_OUT="$VIDEODIR/${NAME}.mp4"
if [ -f "$HERE/sir_complete.mp4" ]; then
  mv "$HERE/sir_complete.mp4" "$VIDEO_OUT"
  echo ""
  echo "✅ Video generado exitosamente:"
  echo "   📹 $VIDEO_OUT"
else
  echo "⚠️  Video no encontrado en ubicación esperada"
fi

# =======================
# RESUMEN FINAL
# =======================
echo ""
echo "=========================================="
echo "🎉 Proceso completado"
echo "=========================================="
echo ""
echo "📁 Archivos generados:"
echo "   Frames CSV: ${OUTDIR}/${NAME}_*.csv (${FRAME_COUNT} archivos)"
echo "   Video:      ${VIDEO_OUT}"
echo ""
echo "💡 Para visualizar:"
echo "   mpv \"${VIDEO_OUT}\""
echo "   # o simplemente abre el archivo en tu reproductor"
echo ""
echo "🧹 Para limpiar frames y liberar espacio:"
echo "   rm ${OUTDIR}/${NAME}*.csv"
echo ""
