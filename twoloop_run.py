import sys,time,os
ENV=lambda k,d: type(d)(os.environ.get(k,d))
from twoloop_lib import *
Q=float(sys.argv[1]) if len(sys.argv)>1 else 2.0       # q^2 = Q^2 ; Q=2 -> q^2=4 m^2
Lam=float(sys.argv[2]) if len(sys.argv)>2 else 8.0
Rfac=float(sys.argv[3]) if len(sys.argv)>3 else 60.0
kind=sys.argv[4] if len(sys.argv)>4 else 'single'
set_reg(kind,Lam)
mu=0.0
nodes=inner_nodes_log(ENV('NU',14),ENV('NPAN',10),ENV('NPER',6)) if kind=='gauss' else inner_nodes(ENV('NU',14),ENV('NV',22))
sos=sigma_os(mu,Lam)
h=ENV('H',0.05)
def mkts(hh): return np.array([Q,0.0,hh,-hh,2*hh,-2*hh])
ts=mkts(h)
def traces(k):
    global ts
    ts=mkts(h)
    """returns f1(t), f2(t) arrays over ts for outer momentum k (4,)"""
    f1=np.zeros(len(ts)); f2=np.zeros(len(ts))
    Sk=S(k); Stk=Sigma_tilde(k,mu,Lam,sos)
    dSk=-Sk@Stk@Sk
    for i,t in enumerate(ts):
        q=np.array([0,0,0,t]); kq=k+q
        Skq=S(kq); dSkq=-Skq@Sigma_tilde(kq,mu,Lam,sos)@Skq
        f1[i]=sum(np.trace(G[n]@Sk@G[n]@Skq) for n in range(4)).real
        f2[i]=sum(np.trace(G[n]@dSk@G[n]@Skq)+np.trace(G[n]@Sk@G[n]@dSkq) for n in range(4)).real
        V=vertex_contracted(k,kq,mu,Lam,nodes,None)       # Lambda_nu(p'=k, p=k+q)
        f2[i]+=sum(np.trace(G[n]@Sk@V[n]@Skq) for n in range(4)).real
    return f1,f2
def subtract(f,hh):
    fQ,f0,fh,fmh,f2h,fm2h=f
    d1=(-f2h+8*fh-8*fmh+fm2h)/(12*hh); d2=(-f2h+16*fh-30*f0+16*fmh-fm2h)/(12*hh*hh)
    return fQ-f0-Q*d1-0.5*Q*Q*d2
# outer grid: s=ln r composite Gauss, theta Gauss
smin,smax=-5.0,np.log(Rfac*max(Lam,1.0))
npan=ENV('OPAN',14); sn,sw=[],[]
edges=np.linspace(smin,smax,npan+1)
for a,b in zip(edges[:-1],edges[1:]):
    x,w=gl(ENV('ONOD',8),a,b); sn+=list(x); sw+=list(w)
sn=np.array(sn); sw=np.array(sw)
th,thw=gl(ENV('NTH',14),0,np.pi)
tot1=0; tot2=0; t0=time.time(); tail=[]
for s_,w_ in zip(sn,sw):
    r=np.exp(s_); a1=0; a2=0
    for t_,tw in zip(th,thw):
        k=np.array([r*np.sin(t_),0,0,r*np.cos(t_)])
        f1,f2=traces(k)
        meas=4*np.pi*np.sin(t_)**2*tw
        hh=h; a1+=meas*subtract(f1,hh); a2+=meas*subtract(f2,hh)
    tot1+=w_*r**4*a1; tot2+=w_*r**4*a2; tail.append((r,w_*r**4*a2))
pref=-1/(4*np.pi**2)/(3*Q*Q)
print(f"reg={kind} Q={Q} Lam={Lam} Rmax={np.exp(smax):.0f}  time {time.time()-t0:.0f}s")
print("one-loop  Pi1 =",pref*tot1,"   target 0.1932613 (Q=2)")
print("two-loop  Pi2 =",pref*tot2,"   target 0.6282165 (Q=2)")
print("tail check: contribution of last panel to Pi2:", pref*sum(x for r,x in tail[-8:]), " (log-tail region)")
