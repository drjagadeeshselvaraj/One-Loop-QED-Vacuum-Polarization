#!/usr/bin/env python3
"""Sec. 6 rows 1-3: closed forms of Laporta-Jentschura [18] (arXiv:2403.07127), m=1.
Eq. (26) one-loop, Eq. (52) two-loop (both in units alpha/pi and (alpha/pi)^2), against
Eq. (6)/(8) (dispersion) and Eq. (55) (imaginary part). Usage: python3 verify_closed_form.py > verify_closed_form_output.txt"""
from mpmath import mp, mpf, mpc, sqrt, log, polylog, zeta, pi, quad, inf, nstr, fabs
mp.dps = 30
def Phi(n, x): return polylog(n, x)+2*polylog(n, -x)
def P1(q2, im=False):   # Eq. (26) divided by alpha/pi
    v = sqrt(1-4/mpc(q2)); L = log((v-1)/(v+1))
    if im: L = log(abs((v-1)/(v+1)))+1j*pi
    return (8-3*v*v)/9-v*(v*v-3)*L/6
def P2(q2, im=False):   # Eq. (52) divided by (alpha/pi)^2
    v = sqrt(1-4/mpc(q2)); x = (v-1)/(v+1)
    L = log(x)
    if im: L = log(abs(x))+1j*pi          # q^2 + i0 above threshold: ln x = ln|x| + i pi
    xr = x.real if im else x
    P = [None, Phi(1, xr), Phi(2, xr), Phi(3, xr)]
    return ((v**4-2*v**2-3)/12*(P[1]*L**2-4*P[2]*L+6*P[3]+3*zeta(3))
            +(3*v-v**3)/12*(4*P[1]*L-4*P[2]+3*L**2)
            +(5*v-3*v**3)/8*L+(7*v**4-22*v**2-33)/96*L**2-13*(v*v-1)/24+mpf(5)/24)
# --- dispersion reference (Eqs. 8, 28, 55) ---
def Im1(s): v = sqrt(1-4/s); return pi*v*(3-v*v)/6
def Im2(s):
    v = sqrt(1-4/s); x = (v-1)/(v+1); L = log((1-v)/(1+v)); P1_ = Phi(1, x); P2_ = Phi(2, x)
    br = (v**4-2*v**2-3)/6*(P1_*L-2*P2_)+(3*v-v**3)/6*(2*P1_+3*L)+(7*v**4-22*v**2-33)/48*L+(5*v-3*v**3)/8
    return pi*br
def F(Im, Q2): return Q2/pi*quad(lambda s: Im(s)/(s*(s+Q2)), [4, 5, 10, 100, 1e4, 1e6, 1e9, inf])

print('== Row 1: Eq. (26) of [18] (spacelike q^2=-Q^2, sign flipped) vs Eq. (6)/(8) ==')
for Q2 in ['0.5', '4', '100']:
    Q2 = mpf(Q2); a = -P1(-Q2).real; b = F(Im1, Q2)
    print(f'Q2={Q2}: -Eq.(26) = {nstr(a,12)}   dispersion F1 = {nstr(b,12)}   diff = {nstr(fabs(a-b),2)}')
print('(paper row 1: 0.03166956943, 0.1932612712, 0.9984904260)')
print('\n== Row 2: Eq. (52) of [18] (spacelike, sign flipped) vs dispersion integral of Eq. (55) ==')
for Q2 in ['0.5', '4', '100', '10000']:
    Q2 = mpf(Q2); a = -P2(-Q2).real; b = F(Im2, Q2)
    print(f'Q2={Q2}: -Eq.(52) = {nstr(a,14)}   dispersion F2 = {nstr(b,14)}   rel.diff = {nstr(fabs(a-b)/fabs(b),2)}')
print('(paper row 2: 0.11704351288, 0.62821652018, 2.0137838422, 3.2935482857; relative difference below 7e-29)')
print('\n== Row 3: Im of Eq. (52) continued to q^2 = s + i0 vs Eq. (55) ==')
for s in ['5', '10', '100']:
    s = mpf(s); a = P2(s, im=True).imag; b = Im2(s)
    print(f's={s}: Im Eq.(52) = {nstr(a,14)}   pi*{{Eq.(55)}} = {nstr(b,14)}   diff = {nstr(fabs(a-b),2)}')
print('(paper row 3: 4.9045691889, 2.3345709382, 0.88928064128)')
print('Re part check at s=5 (Eq. 52 real part, for reference):', nstr(P2(mpf(5), im=True).real, 12))
