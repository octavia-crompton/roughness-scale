#!/usr/bin/env bash
# watch_registry.sh — Watch for changes to figure_registry_concise.txt files
# and auto-rebuild the corresponding LaTeX documents.
#
# Usage:
#   ./latex/watch_registry.sh            # foreground, Ctrl-C to stop
#   ./latex/watch_registry.sh &          # background
#
# Requires: fswatch (brew install fswatch) and pdflatex

set -euo pipefail

# MacTeX path (macOS)
[[ -d /Library/TeX/texbin ]] && export PATH="/Library/TeX/texbin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
FIGURES_DIR="$ROOT_DIR/figures"
BUILD_SCRIPT="$SCRIPT_DIR/build_figure_doc.py"

# ── check dependencies ─────────────────────────────────────────────────────────
if ! command -v fswatch &>/dev/null; then
    echo "fswatch not found. Install with:  brew install fswatch"
    exit 1
fi
if ! command -v pdflatex &>/dev/null; then
    echo "⚠  pdflatex not found — .tex files will be generated but not compiled."
fi

echo "Watching for registry changes in $FIGURES_DIR …"
echo "Press Ctrl-C to stop."
echo ""

# Initial build
python3 "$BUILD_SCRIPT"

# Watch all concise registry files; on change, rebuild that batch only
fswatch -0 --event Updated "$FIGURES_DIR"/*/figure_registry_concise.txt \
| while IFS= read -r -d '' changed_file; do
    batch="$(basename "$(dirname "$changed_file")")"
    echo ""
    echo "──────────────────────────────────────────────"
    echo "Registry changed: $batch  ($(date '+%H:%M:%S'))"
    echo "──────────────────────────────────────────────"
    python3 "$BUILD_SCRIPT" "$batch"
done
