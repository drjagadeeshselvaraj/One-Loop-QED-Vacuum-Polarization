#!/usr/bin/env python3
"""Additional replication checks: Sec. 5.4 (direct 4D integral), sharp-cutoff Table 7/4 rows, Sec. 5.10(c),
and Sec. 6 (two loops) via the Kallen-Sabry spectral function, Eq. (55) of Laporta-Jentschura (arXiv:2403.07127).
Units m=1 unless stated. Usage: python3 verify_extra.py > verify_extra_output.txt"""
from mpmath import mp, mpf, quad, sqrt, pi, log, gamma, polylog, zeta, inf, mpc, nstr, fabs
mp.dps = 30
def H(t): print('\n== '+t+' ==')
al = 1/mpf('137.036'); e2 = 4*pi*al
f = lambda x: x*(1-x)

# ---------- angular integrals over S^3 (sin^2 weight), exact ----------
def J(a, b):
    J0 = pi*(a-sqrt(a*a-b*b))/b**2
    J1 = (pi/2-a*J0)/b
    J2 = (a*a*J0-a*pi/2)/b**2
    return J0, J1, J2

H('Sec. 5.4: direct 4D Euclidean integration, PV combination (weights 1,-1.2,0.2 at m=1,3,7)')
mp.dps = 60
ms = [mpf(1), mpf(3), mpf(7)]; cs = [mpf(1), mpf('-1.2'), mpf('0.2')]
print('sum c =', sum(cs), ' sum c m^2 =', sum(c*m*m for c, m in zip(cs, ms)))
def radial_trace(k, q, m):      # int dOmega3 of  N/(D1 D2) * k^3, N = 8 k.k' + 16 m^2 (g_{mu nu} contraction)
    a = k*k+q*q+m*m; b = 2*k*q
    J0, J1, J2 = J(a, b)
    return 4*pi*k**3/(k*k+m*m)*((8*k*k+16*m*m)*J0+4*b*J1)
def Pi4d_R(Q, R):
    pts = [0, 1, 10, 100, 1e3, 2e3, 5e3]; pts = [p for p in pts if p < R]+[R]
    tot = sum(c*quad(lambda k: radial_trace(k, Q, m), pts) for c, m in zip(cs, ms))
    return e2/(3*Q*Q*(2*pi)**4)*tot
def Pi4d(Q):
    # truncation error ~ R^-2 (+ R^-4): Richardson from R = 2000, 5000 with a + b/R^2
    r1, r2 = Pi4d_R(Q, mpf(2000)), Pi4d_R(Q, mpf(5000))
    return (r2*25e6-r1*4e6)/(25e6-4e6), r1, r2
def Pi_fp(Q):
    return 2*al/pi*sum(c*quad(lambda x: f(x)*log(m*m+Q*Q*f(x)), [0, .5, 1]) for c, m in zip(cs, ms))
for Q in [1, 2]:
    Q = mpf(Q); (d, r1, r2), r = Pi4d(Q), Pi_fp(Q)
    print(f'Q={Q}: raw(R=2000)={nstr(r1,9)} raw(R=5000)={nstr(r2,9)} rel.diff raw={nstr(fabs(r2-r)/fabs(r),2)}; extrapolated={nstr(d,10)}  Feynman-parameter={nstr(r,10)}  rel.diff={nstr(fabs(d-r)/fabs(r),3)}')

H('Longitudinal part J_L = -q_mu q_nu Pi^{mu nu}/e^2, |q|=0.5')
def JL_integrand(k, q, m):
    a = k*k+q*q+m*m; b = 2*k*q; J0, J1, J2 = J(a, b)
    # T_E = 4 q^2 [2k^2 c^2 - k^2 + q k c - m^2]
    ang = 4*q*q*(2*k*k*J2-(k*k+m*m)*J0+q*k*J1)
    return 4*pi*0+k**3/(k*k+m*m)*ang*1   # S^3 measure: 4 pi sin^2 -> factor 4pi applied below
def JL(q, cuts, masses, weights):
    tot = 0
    for c, m in zip(weights, masses):
        up = cuts if cuts is not None else inf
        tot += c*quad(lambda k: JL_integrand(k, q, m), [0, 1, 10, 100, 1e3, 1e4, up] if cuts is None else [0, 1, min(10, cuts), cuts] if cuts <= 10 else [0, 1, 10, cuts])
    return 4*pi*tot/(2*pi)**4
q = mpf('0.5')
def JLR(R):
    pts = [0, 1, 10, 100, 1e3, 2e3, 5e3]; pts = [p for p in pts if p < R]+[R]
    return 4*pi*sum(c*quad(lambda k: JL_integrand(k, q, m), pts) for c, m in zip(cs, ms))/(2*pi)**4
j1, j2 = JLR(mpf(2000)), JLR(mpf(5000))
print('PV combination, R=2000, 5000:', nstr(j1, 4), nstr(j2, 4), ' R^-2 extrapolated:', nstr((j2*25e6-j1*4e6)/21e6, 3), ' (Ward identity: 0)')
for L in [10, 100]:
    L = mpf(L); v = JL(q, L, [mpf(1)], [mpf(1)])
    pred = (1/(8*pi**2))*L**4/(L**2+1)*q*q   # magnitude from Eq. (18)
    print(f'sharp cutoff L={int(L)}: J_L = {nstr(v,8)}   |Eq.(18) prediction| = {nstr(pred,8)}   ratio |J_L|/pred = {nstr(fabs(v)/pred,8)}')

mp.dps = 30
H('Sec. 5.10(c): cancellation of Delta^{D/2-1} terms, non-integer D')
for D in [mpf('3.3'), mpf('3.77'), mpf('2.6')]:
    Dl = mpf('1.7'); q2 = mpf('-0.9'); fx = mpf('0.21')
    pref = 1/(4*pi)**(D/2)
    I0 = gamma(2-D/2)*Dl**(D/2-2)*pref            # i factored out
    I2 = -(D/2)*gamma(1-D/2)*Dl**(D/2-1)*pref
    A = 4*((2/D-1)*I2+(Dl+2*fx*q2)*I0); B = -8*fx*I0
    print(f'D={D}: A+q^2 B... |A - 8 f q^2 I0| = {nstr(fabs(A-8*fx*q2*I0),3)}  (A=-q^2B: {nstr(fabs(A+q2*B),3)})')

H('Sec. 6: two-loop via Kallen-Sabry spectral function, Eq. (55) of [18] (m=1)')
def Phi(n, x): return polylog(n, x)+2*polylog(n, -x)
def Im2(s):   # returns Im Pi^(2) / (alpha/pi)^2 = pi * {...}   (alpha^2/pi = (alpha/pi)^2 * pi)
    v = sqrt(1-4/s); x = (v-1)/(v+1); L = log((1-v)/(1+v)); P1 = -log(1-x); P2 = Phi(2, x); P1 = Phi(1, x)
    br = (v**4-2*v**2-3)/6*(P1*L-2*P2)+(3*v-v**3)/6*(2*P1+3*L)+(7*v**4-22*v**2-33)/48*L+(5*v-3*v**3)/8
    return pi*br
def Im1(s):
    v = sqrt(1-4/s); return pi*v*(3-v*v)/6
def F_spacelike(Im, Q2):
    return Q2/pi*quad(lambda s: Im(s)/(s*(s+Q2)), [4, 5, 10, 100, 1e4, 1e6, 1e9, inf])
print('F1(0.5) =', nstr(F_spacelike(Im1, mpf('0.5')), 11), ' (paper: 0.03166956943)')
for Q2 in ['0.5', '4', '100', '10000']:
    print(f'F2({Q2}) =', nstr(F_spacelike(Im2, mpf(Q2)), 11))
print('(paper rows: 0.11704351288, 0.62821652018, 2.0137838422, 3.2935482857)')
Q2 = mpf('1e-4'); print('row 4: F2/Q2 at Q2=1e-4 =', nstr(F_spacelike(Im2, Q2)/Q2, 8), ' 41/162 =', nstr(mpf(41)/162, 10))
Q2 = mpf('1e7'); F2 = F_spacelike(Im2, Q2); asy = log(Q2)/4+zeta(3)-mpf(5)/24
print('row 5: F2(1e7) =', nstr(F2, 13), ' Eq.(27) without mass term =', nstr(asy, 13), ' with -3 ln Q2/Q2:', nstr(asy-3*log(Q2)/Q2, 13))
print('row 6: (F2-asy)*Q2/ln Q2 =', nstr((F2-asy)*Q2/log(Q2), 8), ' (expected -3)')
