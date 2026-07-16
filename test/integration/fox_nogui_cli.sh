#!/usr/bin/env bash
set -eu

ROOT_DIR="$1"
FOX_NOGUI="$ROOT_DIR/Fox/src/Fox-nogui"
TMP_DIR="$ROOT_DIR/test/integration/.tmp"

# shellcheck source=/dev/null
. "$ROOT_DIR/test/scripts/test_runner.sh"

rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

if [ ! -x "$FOX_NOGUI" ]; then
  echo "Missing executable: $FOX_NOGUI"
  echo "Build it first with: make -f gnu.mak -C $ROOT_DIR/Fox Fox-nogui"
  exit 1
fi

TEST_TIMEOUT="${FOX_TEST_TIMEOUT:-180s}"

if [ "${CI:-}" = "true" ] || [ "${FOX_SKIP_CIMETIDINE_CLI_TEST:-}" = "1" ]; then
  print_result_line "integration::cimetidine-cli" "SKIP" ""
else
  run_test "integration::cimetidine-cli" \
    bash -c '"$1" "$2/Fox/example/Cimetidine-powder.xml" --nogui --silent -n 100 -o "$3/cimetidine-out.xml" && [ -s "$3/cimetidine-out.xml" ] && grep -q "Cimetidine" "$3/cimetidine-out.xml"' \
    -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR"
fi

# Tutorial-equivalent cimetidine workflow:
# crystal + xray powder pattern + MonteCarlo already encoded in example XML.
run_test "integration::tutorial-cimetidine" \
  bash -c '"$1" "$2/Fox/example/Cimetidine-powder.xml" --nogui --silent --randomize -n 1000 --nbrun 1 -o "$3/tutorial-cimetidine-out.xml" && [ -s "$3/tutorial-cimetidine-out.xml" ] && grep -q "Cimetidine" "$3/tutorial-cimetidine-out.xml"' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR"

# Tutorial-equivalent PbSO4 workflow with two powder patterns (xray + neutron)
# and one MonteCarlo optimization object.
run_test "integration::tutorial-pbso4-joint" \
  bash -c '"$1" "$2/Fox/example/pbso4-joint.xml" --nogui --silent --randomize -n 1000 --nbrun 1 -o "$3/tutorial-pbso4-out.xml" && [ -s "$3/tutorial-pbso4-out.xml" ] && grep -q "PbSO4" "$3/tutorial-pbso4-out.xml"' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR"

# Tutorial-equivalent YLID single-crystal workflow:
# crystal + molecule + DiffractionDataSingleCrystal + MonteCarlo from example XML.
run_test "integration::tutorial-ylid-singlecrystal" \
  bash -c 'gzip -cd "$2/Fox/example/ylid.xml.gz" | sed "s/<Option Name=\"Automatic Least Squares Refinement\" Choice=\"1\" ChoiceName=\"At the end of each run\"\\/>/<Option Name=\"Automatic Least Squares Refinement\" Choice=\"0\" ChoiceName=\"Never\"\\/>/" > "$3/ylid-work.xml" && "$1" "$3/ylid-work.xml" --nogui --silent --randomize -n 1000 --nbrun 1 -o "$3/tutorial-ylid-out.xml" && [ -s "$3/tutorial-ylid-out.xml" ] && grep -q "DiffractionDataSingleCrystal" "$3/tutorial-ylid-out.xml"' \
  -- "$FOX_NOGUI" "$ROOT_DIR" "$TMP_DIR"
