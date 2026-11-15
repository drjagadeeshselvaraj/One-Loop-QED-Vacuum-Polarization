# Specification: two-loop regulator test (option B), for an independent code

Goal. Show that the on-shell-renormalized two-loop photon self-energy is the same in two regularizations,
including the overlapping vertex subdivergences, and equals the Kallen-Sabry result already reproduced in this package.

Diagrams (QED, one fermion of mass m): the electron loop with one internal photon: (a),(b) the two self-energy insertions, (c) the vertex correction. 
In (c) the two vertex subgraphs share propagators; this is the overlapping-divergence region that cannot be tested at one loop.

Regulators to compare
 R1: dimensional regularization, D = 4 - 2 eps, all loops.
 R2: Pauli-Villars on the internal photon (Lambda) AND on the outer electron loop (masses and weights with sum c_i = 0, sum c_i m_i^2 = 0, so that the tensor amplitude is transverse; see check_new_claims.py), Lambda >> m, extrapolated Lambda -> infinity.

Renormalization (identical in both): on-shell mass (pole), on-shell charge (Pi(0)=0), one-loop counterterm insertions into (a),(b),(c) in the forest-formula sense (subdivergences subtracted first).

Quantities to compute: Pi_hat_2(Q^2) = Pi_2(Q^2) - Pi_2(0), in the normalization of Eq. (52) of Laporta-Jentschura [18], m = 1.

Targets (reference values from verify_closed_form.py, relative accuracy 1e-29, spacelike Q^2):
  Q^2 = 0.5    0.11704351288357
  Q^2 = 4      0.62821652018414
  Q^2 = 100    2.013783842206
  Q^2 = 1e4    3.2935482856561

Pass criteria
 1. R1 and R2 each reproduce the targets (R2 at several Lambda, with the expected Lambda^-2 (ln Lambda) approach; extrapolate).
 2. Per-diagram: after subdivergence subtraction, the sum (a)+(b)+(c) is transverse, i.e. the q^mu q^nu part vanishes, for R2 at finite Lambda up to O(Lambda^-2).
 3. Control: replace R2's transverse (PV) outer regulator by a sharp cutoff on the loop momentum. A Lambda^2 term should appear in the longitudinal part (two-loop analogue of Proposition 4). Failing to see it means the test is not sensitive.

Suggested tools: FeynCalc + Package-X or LoopTools for R1; pySecDec or numerical Feynman-parameter integration for R2. Cross-check at least one diagram in a second package.

Status. Partly done: Section 6.2 of the paper reports a diagram-level computation with two Pauli-Villars photon profiles (twoloop_*.py), agreeing with the targets to 3e-5..4e-4. The sharp-cutoff control (criterion 3) and a Gaussian proper-time photon regulator were added. Still open: the dimensional-regularization run and a cross-check in an independent package. Partial one-loop ingredients are in sigma_R_pv.py, sigma_R_pv_wide.py and vertex_R_pv.py.
