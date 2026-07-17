#!/usr/bin/env bash
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# shellcheck source=/dev/null
. "$ROOT_DIR/test/scripts/test_runner.sh"

TEST_TIMEOUT="${FOX_TEST_TIMEOUT:-30s}"

run_test "unit::unit_cell_smoke" "$SCRIPT_DIR/unit_cell_smoke"
run_test "unit::crystallography_workflow" "$SCRIPT_DIR/crystallography_workflow"

# --- api_spacegroup ---
run_test "unit::spacegroup" "$SCRIPT_DIR/api_spacegroup" "spacegroup"
run_test "unit::spacegroup-alternate-settings" "$SCRIPT_DIR/api_spacegroup" "spacegroup-alternate-settings"
run_test "unit::spacegroup-reflection-properties" "$SCRIPT_DIR/api_spacegroup" "spacegroup-reflection-properties"
run_test "unit::spacegroup-symmetry-operations" "$SCRIPT_DIR/api_spacegroup" "spacegroup-symmetry-operations"
run_test "unit::spacegroup-asymmetric-unit" "$SCRIPT_DIR/api_spacegroup" "spacegroup-asymmetric-unit"

# --- api_crystal ---
run_test "unit::scattering-power-atom" "$SCRIPT_DIR/api_crystal" "scattering-power-atom"
run_test "unit::crystal-atom" "$SCRIPT_DIR/api_crystal" "crystal-atom"
run_test "unit::crystal-scatterer-management" "$SCRIPT_DIR/api_crystal" "crystal-scatterer-management"
run_test "unit::unitcell-geometry" "$SCRIPT_DIR/api_crystal" "unitcell-geometry"

# --- api_molecule ---
run_test "unit::molecule-atoms-bonds" "$SCRIPT_DIR/api_molecule" "molecule-atoms-bonds"
run_test "unit::molecule-angles-dihedrals" "$SCRIPT_DIR/api_molecule" "molecule-angles-dihedrals"
run_test "unit::molecule-formula-loglikelihood" "$SCRIPT_DIR/api_molecule" "molecule-formula-loglikelihood"

# --- api_scattering ---
run_test "unit::scatteringdata-singlecrystal" "$SCRIPT_DIR/api_scattering" "scatteringdata-singlecrystal"
run_test "unit::scatteringdata-radiation-types" "$SCRIPT_DIR/api_scattering" "scatteringdata-radiation-types"
run_test "unit::diffractiondata-observed" "$SCRIPT_DIR/api_scattering" "diffractiondata-observed"
run_test "unit::singlecrystal-groundtruth-xray" "$SCRIPT_DIR/api_scattering" "singlecrystal-groundtruth-xray"
run_test "unit::singlecrystal-groundtruth-neutron" "$SCRIPT_DIR/api_scattering" "singlecrystal-groundtruth-neutron"

# --- api_cif ---
run_test "unit::cif-import" "$SCRIPT_DIR/api_cif" "cif-import"
run_test "unit::cif-data-fields" "$SCRIPT_DIR/api_cif" "cif-data-fields"
run_test "unit::cif-coordinate-conversion" "$SCRIPT_DIR/api_cif" "cif-coordinate-conversion"

# --- api_optimization ---
run_test "unit::refinablepar" "$SCRIPT_DIR/api_optimization" "refinablepar"
run_test "unit::refinableobj" "$SCRIPT_DIR/api_optimization" "refinableobj"
run_test "unit::optimizationobj" "$SCRIPT_DIR/api_optimization" "optimizationobj"
run_test "unit::optimizationobj-limits-options" "$SCRIPT_DIR/api_optimization" "optimizationobj-limits-options"
run_test "unit::lsqnumobj" "$SCRIPT_DIR/api_optimization" "lsqnumobj"
run_test "unit::lsqnumobj-residual-statistics" "$SCRIPT_DIR/api_optimization" "lsqnumobj-residual-statistics"

# --- api_indexing ---
run_test "unit::peaklist-simulate-volume" "$SCRIPT_DIR/api_indexing" "peaklist-simulate-volume"
run_test "unit::peaklist-add-remove" "$SCRIPT_DIR/api_indexing" "peaklist-add-remove"
run_test "unit::cellexplorer" "$SCRIPT_DIR/api_indexing" "cellexplorer"
run_test "unit::cellexplorer-configuration" "$SCRIPT_DIR/api_indexing" "cellexplorer-configuration"
run_test "unit::cellexplorer-dicvol-tetragonal" "$SCRIPT_DIR/api_indexing" "cellexplorer-dicvol-tetragonal"
run_test "unit::cellexplorer-dicvol-monoclinic" "$SCRIPT_DIR/api_indexing" "cellexplorer-dicvol-monoclinic"

# --- api_powderpattern ---
run_test "unit::powderpattern-background" "$SCRIPT_DIR/api_powderpattern" "powderpattern-background"
run_test "unit::powderpattern-diffraction" "$SCRIPT_DIR/api_powderpattern" "powderpattern-diffraction"
run_test "unit::powderpattern-diffraction-mur" "$SCRIPT_DIR/api_powderpattern" "powderpattern-diffraction-mur"
run_test "unit::powderpattern-import" "$SCRIPT_DIR/api_powderpattern" "powderpattern-import"
run_test "unit::scatteringcorr-subclasses" "$SCRIPT_DIR/api_powderpattern" "scatteringcorr-subclasses"
run_test "unit::reflectionprofile-pseudo-voigt" "$SCRIPT_DIR/api_powderpattern" "reflectionprofile-pseudo-voigt"
run_test "unit::reflectionprofile-double-exponential-pv" "$SCRIPT_DIR/api_powderpattern" "reflectionprofile-double-exponential-pv"
run_test "unit::powder-groundtruth-xray-pv-gaussian" "$SCRIPT_DIR/api_powderpattern" "powder-groundtruth-xray-pv-gaussian"
run_test "unit::powder-groundtruth-xray-pv-lorentzian" "$SCRIPT_DIR/api_powderpattern" "powder-groundtruth-xray-pv-lorentzian"
if [ "${CI:-}" = "true" ]; then
  TEST_LIVE_OUTPUT=1 OBJCRYST_TEST_TRACE=1 run_test "unit::reflectionprofile-pv-anisotropic-direct" "$SCRIPT_DIR/api_powderpattern" "reflectionprofile-pv-anisotropic-direct"
  TEST_LIVE_OUTPUT=1 OBJCRYST_TEST_TRACE=1 run_test "unit::powder-groundtruth-xray-anisotropic" "$SCRIPT_DIR/api_powderpattern" "powder-groundtruth-xray-anisotropic"
else
  run_test "unit::reflectionprofile-pv-anisotropic-direct" "$SCRIPT_DIR/api_powderpattern" "reflectionprofile-pv-anisotropic-direct"
  run_test "unit::powder-groundtruth-xray-anisotropic" "$SCRIPT_DIR/api_powderpattern" "powder-groundtruth-xray-anisotropic"
fi
run_test "unit::powder-groundtruth-neutron-pv-gaussian" "$SCRIPT_DIR/api_powderpattern" "powder-groundtruth-neutron-pv-gaussian"
run_test "unit::powder-groundtruth-neutron-pv-lorentzian" "$SCRIPT_DIR/api_powderpattern" "powder-groundtruth-neutron-pv-lorentzian"
