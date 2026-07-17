#!/usr/bin/env bash
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# shellcheck source=/dev/null
. "$ROOT_DIR/test/scripts/test_runner.sh"

TEST_TIMEOUT="${FOX_TEST_TIMEOUT:-30s}"

# Cimetidine tutorial (Fox/doc/source/tutorial-cimetidine.rst), reproduced
# entirely with API calls (crystal creation, molecule import from a
# Fenske-Hall Z-matrix, powder pattern import, automatic Bayesian
# background, crystalline phase, and a short global optimization run),
# rather than loading a pre-built XML file. See fox_nogui_cli.sh for the
# complementary XML-driven variant of the same tutorial.
run_test "integration::tutorial-cimetidine-api-build" "$SCRIPT_DIR/tutorial_cimetidine" "$ROOT_DIR"
