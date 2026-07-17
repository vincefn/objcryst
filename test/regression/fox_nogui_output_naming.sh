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

TEST_TIMEOUT="${FOX_TEST_TIMEOUT:-30s}"

# Fox-nogui's "-o" output filename argument supports "#cost" / "#Rwp"
# placeholders (see Fox/src/Fox.cpp: filenameInsertCost/filenameInsertRwp),
# which get substituted with the actual final fit quality (e.g.
# "-Cost-<loglikelihood>" / "-Rwp-<value>") right before the result XML is
# saved. This lets scripted/batch runs name their output after the result
# without clashing on a fixed filename.
#
# This test checks that the substitution actually happens end-to-end:
#   1. run Fox-nogui with "-o result-#cost.xml" (a single, short trial)
#   2. the literal, unexpanded filename "result-#cost.xml" must NOT exist
#   3. exactly one "result-*" file must have been created
#   4. that file's name must not contain a literal "#" (i.e. the "#cost"
#      placeholder was properly replaced by a real value, not left as-is)
run_test "regression::output-placeholder-expansion" \
  bash -c 'cd "$3" && "$1" "$2/Fox/example/Cimetidine-powder.xml" --nogui --silent -n 1 -o result-#cost.xml && [ ! -e "$3/result-#cost.xml" ] && [ "$(find "$3" -maxdepth 1 -type f -name "result-*" | wc -l | tr -d " ")" = "1" ] && ! find "$3" -maxdepth 1 -type f -name "result-*" | grep -q "#"' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR"
