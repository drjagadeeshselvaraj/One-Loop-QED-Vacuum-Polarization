# Replication files: Regulator Independence of the Renormalized One-Loop QED Vacuum Polarization

Independent re-implementation, from the paper's formulas, of its numerical checks. The final paper (v21, clean) is in paper/. The authors' original
scripts (verify_regulator_independence.py, verify_two_loop.py) were not available.
Author: Jagadeesh Selvaraj, NextGen Technologies Private Limited, Chennai, India (ORCID 0009-0007-2314-4057).
Install: `pip install mpmath numpy`. Each script writes to stdout; the *_output.txt files are the reference runs.

| Script | Covers | Result |
|---|---|---|
| verify_one_loop.py | Tables 2, 3, 5, 6, 8, 9, 10; Sec. 5.9 (incl. the -5/18 integral); Gaussian rows of Table 7; Sec. 5.10 (a), (b), (d) | all reproduce |
| verify_extra.py | Sec. 5.4 direct 4D integral (agreement 3e-10 and 8e-11, better than the paper's 5e-9); sharp-cutoff J_L rows (Tables 4/7); Sec. 5.10(c); Sec. 6 rows 1, 2, 4, 5, 6 via the Kallen-Sabry spectral function, Eq. (55) of Laporta-Jentschura | all reproduce |
| verify_two_loop_Z.py | Sec. 6 row 7: two-loop leptonic Delta alpha(MZ^2) = 0.37983, 0.23600, 0.16034 (sum 0.77617) e-4 via a timelike principal-value dispersion integral | reproduces. Run after verify_extra.py (reads its Eq. (55) code) |
| verify_closed_form.py | Sec. 6 rows 1-3: closed forms Eq. (26) and Eq. (52) of Laporta-Jentschura [18] against the dispersion integrals (relative differences below 3e-29, matching the paper's 7e-29 bound); imaginary part of Eq. (52) continued to s+i0 vs Eq. (55) at s=5, 10, 100 (diff ~1e-30) | all reproduce |
| real_data_checks.py | Sec. 5.11: muon g-2 electron loop, muonic-hydrogen Uehling shift, OPAL running of alpha | all reproduce |
| timelike_cutoff.py | Sec. 5.12: sharp proper-time cutoff at timelike momentum (real error = continuation of Eq. (13), imaginary error 0) | confirmed |
| sigma_R_pv.py | Sec. 6.1: renormalized one-loop self-energy with a Pauli-Villars photon converges to the Lambda->inf limit as Lambda^-2, subtractions vanish | confirmed |
| sigma_R_pv_wide.py | Sec. 6.1: same self-energy test at spacelike and above-threshold p^2 (-1e4 to 1e4 m^2): error exactly O(Lambda^-2), ~ -(alpha/2pi)p^2/Lambda^2 (identity part) | confirmed |
| vertex_R_pv.py | Sec. 6.1: on-shell vertex F1 with a Pauli-Villars photon reproduces Bethe's slope; error vs Lambda->inf is O(Lambda^-2 ln Lambda). Slow (about 10 min); reference output is a truncated run | confirmed |
| twoloop_lib.py, twoloop_checks.py, twoloop_run.py | Sec. 6.2: diagram-level two-loop vacuum polarization with two Pauli-Villars photon profiles, Taylor-subtracted outer loop. Usage: python3 twoloop_run.py Q Lambda Rfac [single|double|gauss]; each run about 70 s. Results in twoloop_results_output.txt | reproduces Kallen-Sabry to 3e-5..4e-4 |
| vertex_direct_check.py, vertex_direct_gauss.py | Sec. 6.2: one-loop vertex vs direct 4D loop integration (off-shell, non-forward): 7e-5 (PV), 2e-6 (Gaussian) | confirmed |
| sharp_control.py | Sec. 6.2: sharp outer cutoff gives a longitudinal R^2 term: one-loop -0.499 R^2, two-loop R^2 (a ln Lambda_in + b) | confirmed |
| two_loop_test_spec.md | Specification for the full two-loop test (not run) | n/a |
| check_new_claims.py | Remark 6 boundary term, (A3) bound for dimensional regularization, three-regulator Pauli-Villars weights | confirmed |

Notes
- Sec. 6 is now fully covered: rows 1-3 in verify_closed_form.py, rows 4-6 in verify_extra.py, row 7 in verify_two_loop_Z.py.
  Continuation for row 3: q^2 + i0 gives ln((v-1)/(v+1)) = ln|.| + i*pi; the polylogarithms have no cut there.
- Table 10: the sum agrees with Steinhauser to 3e-8 (with masses tuned to [14]), not 3e-7.
- Sec. 5.4 and J_L: the radial integral is cut at R=2000, 5000 and the R^-2 truncation is removed by Richardson extrapolation.
