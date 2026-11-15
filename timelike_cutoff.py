# Sharp proper-time cutoff at timelike momentum: Pi_L(q2=s+i0) = (2a/pi) Int x(1-x)[E1(m2/L2) - E1(A/L2)] dx,
# A = m2 - s x(1-x) - i0 (analytic continuation of the spacelike result, which is how the regulated function is defined).
from mpmath import mp, mpf, mpc, e1, quad, log, pi, im, re
mp.dps=30
a=1/mpf('137.036'); m=mpf('0.5109989e-3'); MZ=mpf('91.1876')
def PiR(s):  # renormalized, timelike, s+i0
    f=lambda x: x*(1-x)*log((m**2-s*x*(1-x)-mpc(0,mpf('1e-40')))/m**2)
    return 2*a/pi*quad(f,[0,0.5,1])
def PiL(s,L):
    f=lambda x: x*(1-x)*(e1(m**2/L**2)-e1((m**2-s*x*(1-x)-mpc(0,mpf('1e-40')))/L**2))
    return 2*a/pi*quad(f,[0,0.5,1])
print("sharp cutoff, timelike s; predicted real error = +a s/(15 pi L^2) (analytic continuation of -a Q^2 c/(15 pi L^2), c=1)")
for s in [MZ**2, 4*MZ**2]:
  for L in [mpf(10)*abs(s)**0.5, mpf(100)*abs(s)**0.5, mpf(1000)*abs(s)**0.5]:
    d=PiL(s,L)-PiR(s); pred=a*s/(15*pi*L**2)
    print(f"s={float(s):9.1f} L/sqrt(s)={float(L/abs(s)**0.5):6.0f}  Re err={float(re(d)):.4e}  pred={float(pred):.4e}  ratio={float(re(d)/pred):.5f}  Im err={float(im(d)):.2e} (rel to ImPi {float(im(d)/im(PiR(s))):.1e})")
print("\nspacelike control, Q2=MZ^2:")
for L in [mpf(10)*MZ, mpf(100)*MZ]:
    s=-MZ**2; d=PiL(s,L)-PiR(s) if False else None
