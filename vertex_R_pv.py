# One-loop on-shell vertex form factor F1(q^2), Feynman gauge, photon mass mu, Pauli-Villars photon Lambda.
# Per photon mass M (Peskin Eq. 6.56 structure), with z the Feynman parameter of the photon line, x=(1-z)u, y=(1-z)(1-u):
#   f_M(q2) = (alpha/2pi) Int dz du (1-z) [ -ln D_M + N/D_M ],  D_M = (1-z)^2 (m^2 - u(1-u) q2) + z M^2,
#   N = (1-x)(1-y) q2 + (1-4z+z^2) m^2.     PV: F1 = f_mu - f_Lambda.   Renormalized: F1_R(q2) = F1(q2) - F1(0).
# Lambda->inf limit: F1_R^inf = [f_mu(q2)-f_mu(0)] (the -ln(z Lambda^2) piece is q2-independent).
# Validation: Bethe slope  dF1/dq2 (0) = (alpha/(3 pi m^2)) [ln(m/mu) - 3/8].   Test: rate of convergence in Lambda.
from mpmath import mp, mpf, quad, log
import time
mp.dps=20
al=1/mpf('137.036'); m=mpf(1); mu=mpf('0.001')
def f(M,q2,zb):
    def inner(z):
        def g(u):
            x=(1-z)*u; y=(1-z)*(1-u)
            D=(1-z)**2*(m**2-u*(1-u)*q2)+z*M**2
            N=(1-x)*(1-y)*q2+(1-4*z+z*z)*m**2
            return -log(D)+N/D
        return (1-z)*quad(g,[0,0.5,1])
    return al/(2*mp.pi)*quad(inner,zb)
def zb_for(M): 
    pts=[mpf(0)]+[mpf(10)**k for k in range(-14,0) if mpf(10)**k < 1]+[mpf(1)]
    return pts
t0=time.time()
# validation of the formula via slope
h=mpf('1e-4'); zbm=zb_for(mu)
F=lambda q2: f(mu,q2,zbm)-f(mu,0,zbm)
slope=(F(h)-F(-h))/(2*h); bethe=al/(3*mp.pi)*(log(m/mu)-mpf(3)/8)
print(f"slope check (mu={float(mu)}): numeric {mp.nstr(slope,8)}  Bethe {mp.nstr(bethe,8)}  rel diff {mp.nstr((slope-bethe)/bethe,3)}   [{time.time()-t0:.0f}s]",flush=True)
qs=[mpf(-4),mpf(-1),mpf('0.5'),mpf(3)]
inf={q:f(mu,q,zbm)-f(mu,0,zbm) for q in qs}
print("\nq2/m2   Lambda/m   err=F1R^L-F1R^inf   err*L^2/m^2",flush=True)
for L in [mpf(10)**k for k in (2,3,4)]:
    zbl=zb_for(L)
    for q in qs:
        FL=lambda qq:(f(mu,qq,zbm)-f(mu,0,zbm))-(f(L,qq,zbl)-f(L,0,zbl))
        e=FL(q)-inf[q]
        print(f"{float(q):5.1f} {float(L):9.0e} {mp.nstr(e,5):>14} {mp.nstr(e*L**2,5):>12}",flush=True)
