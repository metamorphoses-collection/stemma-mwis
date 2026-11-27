#!/bin/zsh
# Run the final audit of the report: every machine check, keyed by the result it tests (code/audit.py).
# From the paper directory:
#   zsh code/run_checks.sh            all checks on the saved simulation outputs (about 1.5 min)
#   FULL=1 zsh code/run_checks.sh     first regenerate the simulation outputs that the birth-death checks read (hours)
# audit.py runs the existing check scripts unchanged (verify_claims.py, kappa_quartet.py, ..., cutting_consistency_sim.py;
# logs in data/check_*.log) and the new checks of audit_new_checks.py, and writes data/final_audit_YYYYMMDD.json and
# data/final_audit_table.tex. It exits non-zero if any check fails or any script exits non-zero.
# Inputs outside this directory: the Lai units (data/lai_ombre/units_v2.json, in
# the MC repository) and the Heinrichi CASC data (~/casc_heinrichi, not redistributable) for the data part of the alpha check.
cd "${0:A:h}/.."
regen() { "$@" || { print "FAILED: regeneration step: $*"; exit 1; }; }
if [[ -n $FULL ]]; then
  regen python3 code/bedier_on_birthdeath.py data/bedier_birthdeath.json 300
  regen env E=0.03 SHAPE=0.38 python3 code/bedier_on_birthdeath.py data/bedier_birthdeath_spencer.json 300
  regen env E=0.03 SHAPE=0.38 L=300 ALPHAS=0.6 python3 code/bedier_on_birthdeath.py data/bedier_birthdeath_spencer_L300.json 300
  regen python3 code/bedier_calibrated.py data/bedier_calibrated.json 300
  regen env E=0.03 SHAPE=0.38 python3 code/bedier_calibrated.py data/bedier_calibrated_spencer.json 300
fi
python3 code/audit.py
