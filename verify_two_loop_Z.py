#!/usr/bin/env python3
"""Sec. 6 row 7: two-loop leptonic Delta alpha(MZ^2) from the Kallen-Sabry spectral function, Eq. (55) of [18],
via a principal-value (timelike) dispersion integral. Units 1e-4. Reference (Steinhauser [14], as quoted in the paper):
0.37983, 0.23600, 0.16034; sum 0.77617."""
import re
from mpmath import mp, mpf, quad, pi, log, sqrt, nstr
mp.dps = 20
src = open('verify_extra.py').read()
ns = {}
exec("from mpmath import *\nmp.dps=20\n"+src[src.index("def Phi(n, x)"):src.index("def F_spacelike")], ns)
Im2 = ns['Im2']
al = 1/mpf('137.036'); MZ = mpf('91.1876')
def F2_timelike(s0):
    g = lambda s: Im2(s)/s
    a = mpf(4); g0 = g(s0)
    pts = lambda lo, hi: [lo*(hi/lo)**(i/12) for i in range(13)]
    I1 = quad(lambda s: (g(s)-g0)/(s-s0), [a, 5, 50, s0/100, s0/2, s0*0.9, s0*0.99, s0, s0*1.01, s0*1.1, s0*2])
    I2 = g0*log(s0/(s0-a))
    I3 = quad(lambda s: g(s)/(s-s0), pts(2*s0, 1e6*s0)+[mp.inf])
    return -s0/pi*(I1+I2+I3)   # Pi_R(Q^2 = -s0): continuation of Eq. (8) to timelike momentum
tot = 0
for name, m in [('electron', '0.5109989e-3'), ('muon', '0.10566'), ('tau', '1.777')]:
    s0 = MZ**2/mpf(m)**2
    v = (al/pi)**2*F2_timelike(s0)*1e4; tot += v
    print(f'{name:9s} two-loop Delta alpha = {nstr(v,6)} e-4')
print('sum', nstr(tot, 6), 'e-4   (paper/[14]: 0.77617)')
