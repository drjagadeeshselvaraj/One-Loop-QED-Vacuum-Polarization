#!/usr/bin/env python3
"""Replication script for the one-loop results of
'Regulator Independence of the Renormalized One-Loop QED Vacuum Polarization'.
Re-implemented from the formulas in the paper (the authors' original scripts were not available).
Covers Tables 2, 3, 5, 6, 8, 9, 10, Sec. 5.9 asymptotics, Sec. 5.10 (a),(b),(d), Gaussian part of Table 7.
Not covered: Sec. 5.4 (direct 4D integral), Table 7 sharp-cutoff rows, Sec. 6 (two loops; needs Laporta-Jentschura formulas).
Usage: python3 verify_one_loop.py > verify_output.txt     (units GeV, alpha(0)=1/137.036)
"""
import numpy as np
from mpmath import mp, mpf, quad, log, exp, pi, sqrt, e1, euler, fabs, inf, zeta
mp.dps = 30
al = 1/mpf('137.036'); me = mpf('0.000511')
def f(x): return x*(1-x)
def PiR(Q2, m=me):
    return 2*al/pi*quad(lambda x: f(x)*log(1+Q2*f(x)/m**2), [0, .5, 1])
def PiLam(Q2, L, m=me):
    return 2*al/pi*quad(lambda x: f(x)*(e1(m**2/L**2)-e1((m**2+Q2*f(x))/L**2)), [0, .5, 1])
def PiPV(Q2, M, m=me):
    return 2*al/pi*quad(lambda x: f(x)*(log(1+Q2*f(x)/m**2)-log(1+Q2*f(x)/M**2)), [0, .5, 1])
def PiPVw(Q2, Ms, cs, m=me):   # scalar-level PV with weights c_i (sum c_i = -1)
    return 2*al/pi*quad(lambda x: f(x)*(log(1+Q2*f(x)/m**2)+sum(c*log(1+Q2*f(x)/M**2) for M, c in zip(Ms, cs))), [0, .5, 1])
lead = lambda Q2, L: al*Q2/(15*pi*L**2)
def H(t): print('\n== '+t+' ==')

H('Table 2: regulator-removal error at Q^2=0.01 GeV^2')
Q2 = mpf('0.01'); R = PiR(Q2); print('Pi_R =', mp.nstr(R, 10))
for L in [1, 10, 100, 1000]:
    L = mpf(L); a, b = PiLam(Q2, L), PiPV(Q2, L)
    print(f'L={int(L):5d}  PT={mp.nstr(a,10)}  PV={mp.nstr(b,10)}  |dPT|={mp.nstr(fabs(a-R),4)}  |dPV|={mp.nstr(fabs(b-R),4)}  lead={mp.nstr(lead(Q2,L),4)}')

H('Table 3: dispersion integral (Eq. 8) vs Feynman parameters (Eq. 6)')
def ImPi(s): return al/3*sqrt(1-4*me**2/s)*(1+2*me**2/s)
def Pdisp(Q2): return Q2/pi*quad(lambda s: ImPi(s)/(s*(s+Q2)), [4*me**2, 8*me**2, 100*me**2, 1e4*me**2, inf])
for q2 in ['0.01', '1', '100', '10000']:
    q2 = mpf(q2); d, r = Pdisp(q2), PiR(q2)
    print(f'Q2={q2}  disp={mp.nstr(d,12)}  FP={mp.nstr(r,12)}  diff={mp.nstr(fabs(d-r),2)}')

H('Table 5: unsubtracted values, all scales 1 GeV')
def MS(Q2): return 2*al/pi*quad(lambda x: f(x)*log((me**2+Q2*f(x))), [0, .5, 1])
def PT(Q2): return -2*al/pi*quad(lambda x: f(x)*e1((me**2+Q2*f(x))), [0, .5, 1])
def PV(Q2): return 2*al/pi*quad(lambda x: f(x)*(log(me**2+Q2*f(x))-log(1+Q2*f(x))), [0, .5, 1])
for n, F in [('MS-bar', MS), ('proper time', PT), ('Pauli-Villars', PV)]:
    print(f'{n:14s} raw(0)={mp.nstr(F(0),9)}  raw(Q2)={mp.nstr(F(Q2),9)}  diff={mp.nstr(F(Q2)-F(0),8)}')
print('alpha*gamma/(3pi) =', mp.nstr(al*euler/(3*pi), 4))

H('Table 6: smooth profiles, Delta Pi vs -alpha Q^2 c_w/(15 pi L^2)')
def dPi_w(Q2, L, omw, upper):
    mp.dps = 15
    inner = lambda x: quad(lambda t: omw(t)/t*(exp(-t*me**2/L**2)-exp(-t*(me**2+Q2*f(x))/L**2)), upper)
    v = -2*al/pi*quad(lambda x: f(x)*inner(x), [0, .5, 1]); mp.dps = 30; return v
profiles = [('exp(-t^2)', lambda t: exp(-t**2), [0, 1, 10, 40], sqrt(pi)/2),
            ('1/(1+t^2)', lambda t: 1/(1+t**2), [0, 1, 10, 100, 1e4, 1e6], pi/2),
            ('sharp', lambda t: mpf(1) if t < 1 else mpf(0), [0, 1], mpf(1))]
for n, om, up, c in profiles:
    for L in [10, 100]:
        L = mpf(L); m_ = dPi_w(Q2, L, om, up); p = -lead(Q2, L)*c
        print(f'{n:10s} c_w={mp.nstr(c,5)} L={int(L):4d} measured={mp.nstr(m_,6)} predicted={mp.nstr(p,6)} ratio={mp.nstr(m_/p,6)}')

H('Table 8 / Cor. 3.1: extrapolation and improved PV weights')
for L1 in [1, 2, 5, 10]:
    L1 = mpf(L1); L2 = 2*L1
    ext = (L2**2*PiLam(Q2, L2)-L1**2*PiLam(Q2, L1))/(L2**2-L1**2)
    pred = al*Q2/(2*pi*L1**2*L2**2)*(me**2/15+Q2/140)
    print(f'L1={int(L1):3d} |err Pi(L1)|={mp.nstr(fabs(PiLam(Q2,L1)-R),4)}  |err ext|={mp.nstr(fabs(ext-R),4)}  Eq.17={mp.nstr(pred,4)}')
e1_ = PiPV(Q2, mpf(10))-R; e2_ = PiPVw(Q2, [mpf(10), mpf(20)], [mpf(1)/3, -mpf(4)/3])-R
print('PV single M=10:', mp.nstr(e1_, 6), '  weights 1/3,-4/3 at 10,20:', mp.nstr(e2_, 5), '  ratio', mp.nstr(e1_/e2_, 3))

H('Table 9: k and Lambda/Q for target accuracy (electron loop, c_w=1)')
for Q in ['0.1', '1', '10']:
    Q = mpf(Q); PR = PiR(Q**2); k = al/(15*pi*PR)
    print(f'Q={Q} Pi_R={mp.nstr(PR,6)} k={mp.nstr(k,4)}  L/Q for 1e-2,1e-4,1e-6: ' + ', '.join(mp.nstr(sqrt(k/t), 3) for t in [1e-2, 1e-4, 1e-6]))

H('Sec. 5.9: asymptotics (Eq. 24)')
print('int x(1-x)ln[x(1-x)] =', mp.nstr(quad(lambda x: f(x)*log(f(x)), [0, .5, 1]), 12), ' (-5/18 =', mp.nstr(mpf(-5)/18, 12), ')')
for r in [1e2, 1e4, 1e6, 1e8]:
    Q2_ = mpf(r)*me**2; res = PiR(Q2_)-al/(3*pi)*(log(Q2_/me**2)-mpf(5)/3)
    print(f'Q2/m2={r:.0e}  residual*Q2/m2={mp.nstr(res*r,6)}   2alpha/pi={mp.nstr(2*al/pi,6)}')

H('Table 10: leptonic Delta alpha(MZ^2), units 1e-4')
MZ = mpf('91.1876')
def timelike(s, m): return 2*al/pi*quad(lambda x: f(x)*log(fabs(1-s*f(x)/m**2)), [0, .5, 1])
def spacelike(s, m): return PiR(s, m)
for lab, ms, F in [('spacelike', (0.000511, 0.105658, 1.77686), spacelike), ('timelike', (0.000511, 0.105658, 1.77686), timelike),
                   ('timelike (rounded masses)', (0.0005109989, 0.10566, 1.777), timelike)]:
    v = [F(MZ**2, mpf(m))*1e4 for m in ms]; print(f'{lab:28s}', ', '.join(mp.nstr(t, 8) for t in v), ' sum', mp.nstr(sum(v), 8))
print('Steinhauser [14] (as quoted in paper): 174.34653, 91.78419, 48.05934, sum 314.19007')

H('Table 7 (Gaussian part): L_G/Lambda^2, m=Q=1')
def LG(L, m=1, Q=1):
    g = lambda x: (lambda D: quad(lambda s: s*(D+s/2)/(s+D)**2*exp(-s/L**2), [0, 1, L**2, 40*L**2]))(m**2+Q**2*f(x))
    return quad(g, [0, .5, 1])
for L in [10, 100, 1000]:
    L = mpf(L); print(f'L={int(L)}  L_G/L^2={mp.nstr(LG(L)/L**2,7)}   0.5-(m^2+Q^2/6)/(2L^2)={mp.nstr(0.5-(1+mpf(1)/6)/(2*L**2),7)}')

H('Sec. 5.10: tensor algebra (numpy, explicit Dirac matrices)')
rng = np.random.default_rng(1)
I2 = np.eye(2); Z = np.zeros((2, 2)); sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1])
g0 = np.block([[I2, Z], [Z, -I2]]).astype(complex)
gs = [np.block([[Z, s], [-s, Z]]).astype(complex) for s in (sx, sy, sz)]; G = [g0]+gs
metric = np.diag([1., -1, -1, -1])
slash = lambda p: sum(metric[i, i]*p[i]*G[i] for i in range(4))
m = 0.7; w1 = w2 = 0
for _ in range(1000):
    k = rng.normal(size=4); q = rng.normal(size=4)
    S = lambda p: np.linalg.inv(slash(p)-m*np.eye(4))
    lhs = S(k+q)@slash(q)@S(k); rhs = S(k)-S(k+q); w1 = max(w1, np.abs(lhs-rhs).max()/np.abs(rhs).max())
    kk = k+q
    for mu in range(4):
        for nu in range(4):
            tr = np.trace(G[mu]@(slash(k)+m*np.eye(4))@G[nu]@(slash(kk)+m*np.eye(4)))
            dot = k@metric@kk
            ref = 4*(k[mu]*kk[nu]+k[nu]*kk[mu]-metric[mu, nu]*(dot-m**2))   # upper-index components
            tr_up = tr
            w2 = max(w2, abs(tr_up-ref)/(abs(ref)+1e-12))
print('(a) Ward identity Eq.3, max rel residual:', f'{w1:.2e}')
print('(b) trace identity Eq.29, max rel residual:', f'{w2:.2e}')
eps = mpf('1e-9'); D = 4-2*eps
PiDR = lambda Q2_: -2*al/pi*quad(lambda x: f(x)*(mp.gamma(eps)*(4*pi)**eps*(me**2+Q2_*f(x))**(-eps)), [0, .5, 1])
print('(d) DR at eps=1e-9: Pi(Q2)-Pi(0) =', mp.nstr(PiDR(Q2)-PiDR(0), 12), ' vs Eq.6:', mp.nstr(R, 12))
