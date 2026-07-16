#!/usr/bin/env bash
set -eu

ROOT_DIR="$1"
FOX_NOGUI="$ROOT_DIR/Fox/src/Fox-nogui"
TMP_DIR="$ROOT_DIR/test/regression/.tmp"

# shellcheck source=/dev/null
. "$ROOT_DIR/test/scripts/test_runner.sh"

rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

if [ ! -x "$FOX_NOGUI" ]; then
  echo "Missing executable: $FOX_NOGUI"
  echo "Build it first with: make -f gnu.mak -C $ROOT_DIR/Fox Fox-nogui"
  exit 1
fi

run_test "regression::output-placeholder-expansion" \
  bash -c 'cd "$3" && "$1" "$2/Fox/example/Cimetidine-powder.xml" --nogui --silent -n 1 -o result-#cost.xml && [ ! -e "$3/result-#cost.xml" ] && [ "$(find "$3" -maxdepth 1 -type f -name "result-*" | wc -l | tr -d " ")" = "1" ] && ! find "$3" -maxdepth 1 -type f -name "result-*" | grep -q "#"' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR"
