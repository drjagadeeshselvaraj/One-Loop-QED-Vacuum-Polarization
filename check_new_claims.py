from mpmath import mp, mpf, quad, exp, log, gamma, pi, inf, diff, matrix, lu_solve
mp.dps=30
# Remark 6: tensor-level proper-time boundary term
print("Remark 6: int_{1/L^2}^inf (D/s+1/s^2)e^{-sD} ds  vs  L^2 e^{-D/L^2}")
for D,L in [(1.2,10),(2.5,30),(0.3,100)]:
    D=mpf(D);L=mpf(L)
    I=quad(lambda s:(D/s+1/s**2)*exp(-s*D),[1/L**2,1,50/D,inf])
    print(f"  D={D} L={L}: {mp.nstr(I,12)}  {mp.nstr(L**2*exp(-D/L**2),12)}  L^2-D={mp.nstr(L**2-D,8)}")
# (A3) for DR: g_eps(s) = (2a/pi)(4pi)^e Gamma(1+e) int f^2 Delta^{-1-e}; compare with d/ds Pi_DR(-s)
a=1/mpf('137.036'); m2=mpf(1); f=lambda x:x*(1-x)
for e in [mpf('0.1'),mpf('0.3')]:
    Pi=lambda s: -2*a/pi*(4*pi)**e*gamma(e)*quad(lambda x:f(x)*(m2+s*f(x))**(-e),[0,.5,1])
    g=lambda s: 2*a/pi*(4*pi)**e*gamma(1+e)*quad(lambda x:f(x)**2*(m2+s*f(x))**(-1-e),[0,.5,1])
    s0=mpf(2); print("(A3) DR eps=",e," dPi/ds=",mp.nstr(diff(Pi,s0),12)," g_eps=",mp.nstr(g(s0),12))
    bound=2*a/pi*(4*pi)**e*gamma(1+e)*m2**(-e)*quad(lambda x:f(x)**2/(m2+s0*f(x)),[0,.5,1]); print("   bound",mp.nstr(bound,8),">= g:",bound>=g(s0))
# three-regulator solution
m2=mpf('0.000511')**2; M=[mpf(10),mpf(20),mpf(40)]; y=[t**2 for t in M]
A=matrix([[1,1,1],[y[0],y[1],y[2]],[1/y[0],1/y[1],1/y[2]]]); c=lu_solve(A,matrix([-1,-m2,0]))
print("3-regulator weights:",[mp.nstr(t,10) for t in c], " checks:",mp.nstr(sum(c),8),mp.nstr(sum(ci*yi for ci,yi in zip(c,y))+m2,3),mp.nstr(sum(ci/yi for ci,yi in zip(c,y)),3))
# two-regulator impossibility: residual of 3 eqs best fit
A2=matrix([[1,1],[y[0],y[1]],[1/y[0],1/y[1]]]); 
import mpmath; sol=mpmath.qr_solve(A2,matrix([-1,-m2,0])); print("2-regulator least-squares residual norm:",mp.nstr(sol[1],6))
