#!/usr/bin/env bash
# Configure, build, and install LADlib.
#
# Usage: ./build.sh [-t BUILD_TYPE] [-j JOBS] [-H HCANA_PREFIX] [-c]
#   -t  CMake build type (Debug, Release, RelWithDebInfo; default: RelWithDebInfo)
#   -j  parallel build jobs (default: number of cores)
#   -H  hcana install prefix (default: $HCANA, else ../hcana/install if present)
#   -c  clean: remove the build and install directories first
#
# Debug builds go to build_dbg/ and install_dbg/, all others to build/ and install/.

set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BUILD_TYPE=RelWithDebInfo
JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
HCANA_PREFIX="${HCANA:-}"
CLEAN=0

while getopts "t:j:H:ch" opt; do
  case $opt in
    t) BUILD_TYPE=$OPTARG ;;
    j) JOBS=$OPTARG ;;
    H) HCANA_PREFIX=$OPTARG ;;
    c) CLEAN=1 ;;
    h) sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) exit 1 ;;
  esac
done

if [[ -z "$HCANA_PREFIX" && -d "$SRC_DIR/../hcana/install" ]]; then
  HCANA_PREFIX="$(cd "$SRC_DIR/../hcana/install" && pwd)"
fi

if [[ "$BUILD_TYPE" == Debug ]]; then
  BUILD_DIR="$SRC_DIR/build_dbg"
  INSTALL_DIR="$SRC_DIR/install_dbg"
else
  BUILD_DIR="$SRC_DIR/build"
  INSTALL_DIR="$SRC_DIR/install"
fi

if (( CLEAN )); then
  rm -rf "$BUILD_DIR" "$INSTALL_DIR"
fi

CMAKE_ARGS=(
  -B "$BUILD_DIR" -S "$SRC_DIR"
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR"
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
)
[[ -n "$HCANA_PREFIX" ]] && CMAKE_ARGS+=(-DCMAKE_PREFIX_PATH="$HCANA_PREFIX")

echo "LADlib: $BUILD_TYPE build, $JOBS jobs"
echo "  build:   $BUILD_DIR"
echo "  install: $INSTALL_DIR"
echo "  hcana:   ${HCANA_PREFIX:-<found by CMake>}"

cmake "${CMAKE_ARGS[@]}"
cmake --build "$BUILD_DIR" -j "$JOBS"
cmake --install "$BUILD_DIR"
