#!/usr/bin/env python3
"""Confronting the paper's one-loop formulas (Eq. 6, 8, 9) with published high-precision / measured values.
Usage: python3 real_data_checks.py > real_data_output.txt"""
from mpmath import mp, mpf, quad, log, pi, sqrt, inf, nstr
mp.dps = 25
al = 1/mpf('137.035999177'); me = mpf('0.51099895069'); mmu = mpf('105.6583755'); mtau = mpf('1776.86')
f = lambda x: x*(1-x)
def F1(Q2, m):                                  # Pi_R/(alpha/pi), Eq. (6); Q2, m in same units
    return 2*quad(lambda x: f(x)*log(1+Q2*f(x)/m**2), [0, .5, 1])
def Dalpha_lep(Q2_GeV2):                         # e, mu, tau loops, spacelike
    Q2 = mpf(Q2_GeV2)*1e6                        # -> MeV^2
    return al/pi*sum(F1(Q2, m) for m in (me, mmu, mtau))

print("== A. Muon g-2: electron vacuum-polarization loop, A2(m_mu/m_e) ==")
r = mpf('206.7682838')
A = quad(lambda x: (1-x)*F1(x**2*r**2/(1-x), 1), [0, .5, .9, .99, 1])   # a = (alpha/pi)^2 * A, from Eq. (6) in units m_e=1
print("computed from Eq. (6):", nstr(A, 11), "   literature (Laporta-Remiddi / CODATA eq. 42): 1.0942583111(84)")

print("\n== B. Muonic hydrogen: Uehling (electron-loop) shift of the 2P-2S Lamb shift ==")
x = mpf('0.73738368')                            # m_e/(alpha m_R), quoted in Jentschura (arXiv:1011.5275)
Escale = al*me*1e6/x                              # alpha^2 m_R  in eV  (= alpha m_e / x)
K = lambda t: (1+1/(2*t*t))*sqrt(t*t-1)/t**2      # spectral weight from Im Pi, Eq. (9)
I2P = lambda t: 1/(4*(1+2*x*t)**4)
I2S = lambda t: mpf(1)/2*(1/(1+2*x*t)**2-2/(1+2*x*t)**3+mpf(3)/2/(1+2*x*t)**4)
dE = Escale*(2*al/(3*pi))*quad(lambda t: K(t)*(I2S(t)-I2P(t)), [1, 2, 10, 100, inf])
print("2P-2S Uehling shift from Eq. (9):", nstr(dE*1e3, 8), "meV    literature (arXiv:1011.5275, Eq. 2.24): 205.0074 meV")

print("\n== C. OPAL small-angle Bhabha: change of alpha between t=-1.81 and -6.07 GeV^2 ==")
lep = (Dalpha_lep(6.07)-Dalpha_lep(1.81))*1e5
meas, st, sy, th = 440, 58, 43, 30
err = sqrt(st**2+sy**2+th**2)
print("leptonic-only prediction from Eq. (6):", nstr(lep, 5), "e-5   (OPAL quotes 202e-5)")
print("OPAL measured: (440 +- 58 stat +- 43 syst +- 30 theory) e-5;  total error", nstr(err, 4), "e-5;  SM full prediction quoted: 460 +- 8 e-5")
print("pull of leptonic-only vs data:", nstr((meas-lep)/err, 3), "sigma   (OPAL: 'about 3 sigma')")
print("pull of SM (leptonic+hadronic, 460) vs data:", nstr((meas-460)/err, 3), "sigma")
print("hadronic part implied by the data:", nstr(meas-lep, 4), "+-", nstr(err, 3), "e-5")
