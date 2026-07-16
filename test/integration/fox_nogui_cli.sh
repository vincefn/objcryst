#!/usr/bin/env bash
set -eu

ROOT_DIR="$1"
FOX_NOGUI="$ROOT_DIR/Fox/src/Fox-nogui"
TMP_DIR="$ROOT_DIR/test/integration/.tmp"
TIMEOUT_BIN=""
if command -v timeout >/dev/null 2>&1; then
  TIMEOUT_BIN="timeout"
elif command -v gtimeout >/dev/null 2>&1; then
  TIMEOUT_BIN="gtimeout"
fi

# shellcheck source=/dev/null
. "$ROOT_DIR/test/scripts/test_runner.sh"

rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

if [ ! -x "$FOX_NOGUI" ]; then
  echo "Missing executable: $FOX_NOGUI"
  echo "Build it first with: make -f gnu.mak -C $ROOT_DIR/Fox Fox-nogui"
  exit 1
fi

run_test "integration::cimetidine-cli" \
  bash -c 'if [ -n "$4" ]; then "$4" --signal=TERM --kill-after=20s "${FOX_TEST_TIMEOUT:-180s}" "$1" "$2/Fox/example/Cimetidine-powder.xml" --nogui --silent -n 0 -o "$3/cimetidine-out.xml"; else "$1" "$2/Fox/example/Cimetidine-powder.xml" --nogui --silent -n 0 -o "$3/cimetidine-out.xml"; fi && [ -s "$3/cimetidine-out.xml" ] && grep -q "Cimetidine" "$3/cimetidine-out.xml"' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR" "$TIMEOUT_BIN"

# Tutorial dataset: cimetidine tutorial XML with local relative files.
run_test "integration::tutorial-cimetidine" \
  bash -c 'cd "$2/Fox/example/tutorial-cimetidine" && if [ -n "$4" ]; then "$4" --signal=TERM --kill-after=20s "${FOX_TEST_TIMEOUT:-180s}" "$1" test.xml --nogui --silent -n 0 -o "$3/tutorial-cimetidine-out.xml"; else "$1" test.xml --nogui --silent -n 0 -o "$3/tutorial-cimetidine-out.xml"; fi && [ -s "$3/tutorial-cimetidine-out.xml" ]' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR" "$TIMEOUT_BIN"

# Tutorial dataset: ylid indexing input.
run_test "integration::tutorial-ylid-indexing" \
  bash -c 'if [ -n "$3" ]; then "$3" --signal=TERM --kill-after=20s "${FOX_TEST_TIMEOUT:-180s}" "$1" --index "$2/Fox/example/tutorial-ylid/ylid1.hkl"; else "$1" --index "$2/Fox/example/tutorial-ylid/ylid1.hkl"; fi' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TIMEOUT_BIN"
