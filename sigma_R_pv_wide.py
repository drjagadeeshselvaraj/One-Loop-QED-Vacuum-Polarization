# Self-energy with a Pauli-Villars photon at spacelike and above-threshold p^2.
# Sigma_Lambda(pslash) = (alpha/2pi) Int dx (2m - x pslash) ln D,  D=(1-x)m^2 + x L^2 - x(1-x)p^2 ; ln(x L^2) is polynomial -> drops after renormalization.
# Error vs the Lambda->inf limit is the renormalized Lambda-part of  Int (2m - x P) ln(D/(x L^2))  (independent of the IR photon mass, real for p^2 < L^2).
# Write Sigma = m A(p2) + pslash B(p2).  Renormalized (subtract at pslash=m):  Sigma_R = m S(p2) + pslash beta(p2),
#   K = 2 m^2 (A'+B')(m^2) + B(m^2);  S = A(p2)-A(m^2)-B(m^2)+K ;  beta = B(p2)-K.
from mpmath import mp, mpf, quad, log
mp.dps=25
al=1/mpf('137.036'); m=mpf(1)
def parts(p2,L):
    L2=L*L
    bp=[mpf(0)]+[mpf(10)**k for k in range(-12,0)]+[mpf(1)]
    bp=[b for b in bp if b<=1]
    def Dn(x): return ((1-x)*m**2+x*L2-x*(1-x)*p2)/(x*L2)
    A =al/(2*mp.pi)*quad(lambda x: 2*log(Dn(x)),bp)
    B =al/(2*mp.pi)*quad(lambda x: -x*log(Dn(x)),bp)
    Ap=al/(2*mp.pi)*quad(lambda x: 2*(-x*(1-x))/((1-x)*m**2+x*L2-x*(1-x)*p2),bp)
    Bp=al/(2*mp.pi)*quad(lambda x: -x*(-x*(1-x))/((1-x)*m**2+x*L2-x*(1-x)*p2),bp)
    return A,B,Ap,Bp
def err(p2,L):
    A,B,Ap,Bp=parts(p2,L); A0,B0,Ap0,Bp0=parts(m**2,L)
    K=2*m**2*(Ap0+Bp0)+B0
    return A-A0-B0+K, B-K
print("p2/m2  Lambda/m   S_err*L^2/m^2   beta_err*L^2/m^2      (S: identity part /m, beta: pslash part)")
for p2 in [-10000,-100,-1,2,4,100,10000]:
    for L in [mpf(10)**k for k in (3,4,5,6)]:
        if L*L<=100*abs(p2): continue
        S,b=err(mpf(p2),L)
        print(f"{p2:7d} {float(L):8.0e} {mp.nstr(S*L**2,6):>14} {mp.nstr(b*L**2,6):>16}",flush=True)
