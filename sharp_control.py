# Control: sharp cutoff |k|<R on the OUTER loop (no Taylor subtraction), inner photon PV regulated at large Lambda_in.
# Longitudinal part  Pi_L = q_mu q_nu Pi_{mu nu}/q^2 = -(1/(4 pi^2)) Int_{|k|<R} d^4k  f_L(k),  f_L = tr[g4 (...) g4 (...)]  (q along axis 4).
# One loop expectation (paper Prop. 4 type): Pi_L ~ c1 R^2 ; test two-loop growth.
import sys,time
import twoloop_lib as T
from twoloop_lib import *
Q=2.0; Lam=float(sys.argv[1]) if len(sys.argv)>1 else 2000.0; mu=0.0
T.set_reg('single',Lam); nodes=inner_nodes(14,22); sos=sigma_os(mu,Lam)
def fL(k):
    q=np.array([0,0,0,Q]); kq=k+q
    Sk=S(k); Skq=S(kq); dSk=-Sk@Sigma_tilde(k,mu,Lam,sos)@Sk; dSkq=-Skq@Sigma_tilde(kq,mu,Lam,sos)@Skq
    g4=G[3]
    f1=np.trace(g4@Sk@g4@Skq).real
    f2=(np.trace(g4@dSk@g4@Skq)+np.trace(g4@Sk@g4@dSkq)).real
    V=vertex_contracted(k,kq,mu,Lam,nodes,None)
    f2+=np.trace(g4@Sk@V[3]@Skq).real
    return f1,f2
Rs=[1,2,4,8,16,32]
edges=np.concatenate([[np.log(1e-3)],np.log(Rs)])
th,thw=gl(14,0,np.pi)
cum1=0;cum2=0; t0=time.time(); res=[]
for a,b in zip(edges[:-1],edges[1:]):
    x,w=gl(8,a,b)
    for s_,w_ in zip(x,w):
        r=np.exp(s_); a1=a2=0
        for t_,tw in zip(th,thw):
            k=np.array([r*np.sin(t_),0,0,r*np.cos(t_)]); f1,f2=fL(k); meas=4*np.pi*np.sin(t_)**2*tw
            a1+=meas*f1; a2+=meas*f2
        cum1+=w_*r**4*a1; cum2+=w_*r**4*a2
    R=np.exp(b); res.append((R,-cum1/(4*np.pi**2),-cum2/(4*np.pi**2)))
print(f"Lambda_in={Lam}  Q={Q}  ({time.time()-t0:.0f}s)\n R    Pi_L^(1)    Pi_L^(1)/R^2    Pi_L^(2)    Pi_L^(2)/R^2")
for R,a,b in res: print(f"{R:4.0f} {a:11.5f} {a/R**2:12.5f} {b:11.5f} {b/R**2:12.6f}")
