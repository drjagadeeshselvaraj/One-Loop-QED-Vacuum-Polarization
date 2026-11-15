#!/bin/sh
# One-loop and published-result checks (minutes). verify_two_loop_Z.py reads code from verify_extra.py, so keep this order.
for f in verify_one_loop verify_extra verify_closed_form verify_two_loop_Z check_new_claims real_data_checks timelike_cutoff sigma_R_pv sigma_R_pv_wide vertex_direct_check vertex_direct_gauss sharp_control twoloop_checks; do
  echo "== $f"; python3 $f.py > /tmp/$f.out 2>&1 && echo ok || echo FAILED
done
# Slow, optional (about 10 min / about 70 s per run); reference outputs are included:
#   python3 vertex_R_pv.py
#   python3 twoloop_run.py 0.5 64 40 single      # Q Lambda Rfac [single|double|gauss]
