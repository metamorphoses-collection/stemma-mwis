#!/bin/zsh
# Run every machine check of the audit appendix and report each exit status. From the paper directory:
#   zsh code/run_checks.sh            quick checks only (C1-C25 on saved simulation outputs)
#   FULL=1 zsh code/run_checks.sh     also regenerate the simulation outputs that C21, C22 and C24 read (hours)
# Inputs outside this directory: the Lai units (~/MC/projects/priscus_excerpta/data/historic/lai_ombre/units.json,
# in the repository) and the Heinrichi CASC data (~/casc_heinrichi, not redistributable) for C20's data part.
cd "${0:A:h}/.."
typeset -A rc
run() { local tag=$1; shift; print -n "$tag ... "; "$@" > data/check_$tag.log 2>&1; rc[$tag]=$?; print ${rc[$tag]}; }
run C1-C14 python3 code/verify_claims.py data
run C15 python3 code/kappa_quartet.py data/kappa_quartet.json
run C16 python3 code/mechanisms_check.py data/mechanisms_check.json
run C17 python3 code/drift_quartet.py data/drift_quartet.json
run C18 python3 code/gap_check.py data/gap_check.json
run C19 python3 code/fragmentary_check.py data/fragmentary_check.json
run C20 python3 code/alpha_estimate.py data/alpha_estimate.json
if [[ -n $FULL ]]; then
  python3 code/bedier_on_birthdeath.py data/bedier_birthdeath.json 300
  E=0.03 SHAPE=0.38 python3 code/bedier_on_birthdeath.py data/bedier_birthdeath_spencer.json 300
  E=0.03 SHAPE=0.38 L=300 ALPHAS=0.6 python3 code/bedier_on_birthdeath.py data/bedier_birthdeath_spencer_L300.json 300
  python3 code/bedier_calibrated.py data/bedier_calibrated.json 300
  E=0.03 SHAPE=0.38 python3 code/bedier_calibrated.py data/bedier_calibrated_spencer.json 300
fi
run C21 python3 code/bedier_birthdeath_check.py data/bedier_birthdeath.json data/bedier_birthdeath_spencer.json data/bedier_birthdeath_spencer_L300.json
run C22 python3 code/gamma_confidence.py data/gamma_confidence.json
run C23 python3 code/polytomy_null.py data/polytomy_null.json
run C24 python3 code/bedier_calibrated_check.py data/bedier_calibrated.json data/bedier_calibrated_spencer.json
run C25 python3 code/cutting_consistency.py data/cutting_consistency.json
fails=0; for k v in ${(kv)rc}; do (( v != 0 )) && { print "FAILED: $k"; fails=1; }; done
(( fails == 0 )) && print "ALL CHECKS PASS"; exit $fails
