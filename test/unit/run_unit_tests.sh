#!/usr/bin/env bash
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# shellcheck source=/dev/null
. "$ROOT_DIR/test/scripts/test_runner.sh"

run_test "unit::unit_cell_smoke" "$SCRIPT_DIR/unit_cell_smoke"
run_test "unit::crystallography_workflow" "$SCRIPT_DIR/crystallography_workflow"
run_test "unit::spacegroup" "$SCRIPT_DIR/api_surface_tests" "spacegroup"
run_test "unit::scattering-power-atom" "$SCRIPT_DIR/api_surface_tests" "scattering-power-atom"
run_test "unit::crystal-atom" "$SCRIPT_DIR/api_surface_tests" "crystal-atom"
run_test "unit::molecule" "$SCRIPT_DIR/api_surface_tests" "molecule"
run_test "unit::scatteringdata-singlecrystal" "$SCRIPT_DIR/api_surface_tests" "scatteringdata-singlecrystal"
run_test "unit::powderpattern-background" "$SCRIPT_DIR/api_surface_tests" "powderpattern-background"
run_test "unit::powderpattern-diffraction" "$SCRIPT_DIR/api_surface_tests" "powderpattern-diffraction"
run_test "unit::powderpattern-import" "$SCRIPT_DIR/api_surface_tests" "powderpattern-import"
run_test "unit::cif-import" "$SCRIPT_DIR/api_surface_tests" "cif-import"
run_test "unit::refinablepar" "$SCRIPT_DIR/api_surface_tests" "refinablepar"
run_test "unit::refinableobj" "$SCRIPT_DIR/api_surface_tests" "refinableobj"
run_test "unit::optimizationobj" "$SCRIPT_DIR/api_surface_tests" "optimizationobj"
run_test "unit::lsqnumobj" "$SCRIPT_DIR/api_surface_tests" "lsqnumobj"
run_test "unit::cellexplorer" "$SCRIPT_DIR/api_surface_tests" "cellexplorer"
run_test "unit::singlecrystal-groundtruth-xray" "$SCRIPT_DIR/api_surface_tests" "singlecrystal-groundtruth-xray"
run_test "unit::singlecrystal-groundtruth-neutron" "$SCRIPT_DIR/api_surface_tests" "singlecrystal-groundtruth-neutron"
run_test "unit::powder-groundtruth-xray-pv-gaussian" "$SCRIPT_DIR/api_surface_tests" "powder-groundtruth-xray-pv-gaussian"
run_test "unit::powder-groundtruth-xray-pv-lorentzian" "$SCRIPT_DIR/api_surface_tests" "powder-groundtruth-xray-pv-lorentzian"
if [ "${CI:-}" = "true" ]; then
  TEST_LIVE_OUTPUT=1 OBJCRYST_TEST_TRACE=1 run_test "unit::reflectionprofile-pv-anisotropic-direct" "$SCRIPT_DIR/api_surface_tests" "reflectionprofile-pv-anisotropic-direct"
  TEST_LIVE_OUTPUT=1 OBJCRYST_TEST_TRACE=1 run_test "unit::powder-groundtruth-xray-anisotropic" "$SCRIPT_DIR/api_surface_tests" "powder-groundtruth-xray-anisotropic"
else
  run_test "unit::reflectionprofile-pv-anisotropic-direct" "$SCRIPT_DIR/api_surface_tests" "reflectionprofile-pv-anisotropic-direct"
  run_test "unit::powder-groundtruth-xray-anisotropic" "$SCRIPT_DIR/api_surface_tests" "powder-groundtruth-xray-anisotropic"
fi
run_test "unit::powder-groundtruth-neutron-pv-gaussian" "$SCRIPT_DIR/api_surface_tests" "powder-groundtruth-neutron-pv-gaussian"
run_test "unit::powder-groundtruth-neutron-pv-lorentzian" "$SCRIPT_DIR/api_surface_tests" "powder-groundtruth-neutron-pv-lorentzian"
