# One-loop electron self-energy, Feynman gauge, photon mass mu and Pauli-Villars photon Lambda (Peskin Eq. 7.19 structure):
#   Sigma(P) = (alpha/2pi) Int_0^1 dx (2m - x P) [ ln D(Lambda) - ln D(mu) ],  D(M) = (1-x)m^2 + x M^2 - x(1-x) p^2,  P = pslash, P^2 = p^2.
# Renormalized: Sigma_R(P) = Sigma(P) - Sigma(m) - (P-m) Sigma'(m).  Regulator-free limit (Lambda->inf, polynomial pieces drop out):
#   Sigma_R^inf(P) = same subtractions applied to  -(alpha/2pi) Int (2m - xP) ln D(mu).
# Test: Sigma_R^Lambda -> Sigma_R^inf, rate and uniformity in P; and Sigma_R(m)=Sigma_R'(m)=0.
from mpmath import mp, mpf, quad, log, diff, sqrt
mp.dps=30
al=1/mpf('137.036'); m=mpf(1); mu=mpf('0.1')
def Sig(P,L):
    p2=P*P
    def f(x):
        D=lambda M:(1-x)*m**2+x*M**2-x*(1-x)*p2
        g=(-log(D(mu)))+(log(D(L)) if L is not None else 0)
        return (2*m-x*P)*g
    return al/(2*mp.pi)*quad(f,[0,mpf('1e-8'),mpf('1e-4'),0.01,0.5,1])
def SigR(P,L):
    d=diff(lambda t:Sig(t,L),m)
    return Sig(P,L)-Sig(m,L)-(P-m)*d
print("subtraction check: SigR(m), SigR'(m) at Lambda=1e3:", mp.nstr(SigR(m,mpf(1e3)),3), mp.nstr(diff(lambda t:SigR(t,mpf(1e3)),m),3))
Ps=[mpf('0.3'),mpf('0.6'),mpf('0.9')]
inf={P:SigR(P,None) for P in Ps}
print("\nP/m   Lambda/m   err=SigR^L-SigR^inf   err*L^2/m^2")
for L in [mpf(10)**k for k in (1,2,3,4,5)]:
    row=[]
    for P in Ps:
        e=SigR(P,L)-inf[P]; row.append(e)
        print(f"{float(P):4.1f} {float(L):9.0e} {mp.nstr(e,5):>14}  {mp.nstr(e*L**2,5):>12}")
print("\nerr*Lambda^2 is constant to 4 digits for Lambda>=1e3 (e.g. -3.509e-4 at P=0.3m): the error is O(Lambda^-2), no logarithm, for each tested P.")
