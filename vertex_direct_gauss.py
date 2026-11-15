# Independent check of the one-loop vertex: direct 4D numerical integration over the loop momentum l (no Feynman parameters,
# no gamma-matrix identities) at an off-shell, non-forward point, versus twoloop_lib.vertex_contracted.
# Lambda_nu(p',p) = -(g/(4 pi^2)) Int d^4l  gamma_a S(p'-l) gamma_nu S(p-l) gamma_a * Dreg(l),  Dreg = sum c_i/(l^2+M_i^2)   [units g=1]
import numpy as np, time
import twoloop_lib as T
from twoloop_lib import G,S,I4,gl
T.set_reg('gauss',8.0)
p=np.array([0.6,0.2,0,0]); pp=np.array([0.1,0.8,0,0])          # both in the (e1,e2) plane
# 4D coordinates adapted to the plane: l=(l1,l2,rho cos phi, rho sin phi); l1=r cos t1, l2=r sin t1 cos t2, rho=r sin t1 sin t2
# measure d^4l = r^3 sin^2 t1 sin t2 dr dt1 dt2 dphi
Rmax=8.0*60
edges=np.linspace(np.log(1e-3),np.log(Rmax),9); sn=[];sw=[]
for a,b in zip(edges[:-1],edges[1:]):
    x,w=gl(10,a,b); sn+=list(x); sw+=list(w)
sn=np.array(sn); sw=np.array(sw); r=np.exp(sn)
t1,w1=gl(28,0,np.pi); t2,w2=gl(28,0,np.pi); ph=(np.arange(8)+0.5)*2*np.pi/8; wph=np.full(8,2*np.pi/8)
Dreg=lambda l2: np.exp(-l2/64.0)/l2
tot=np.zeros((4,4,4),complex); t0=time.time()
for ri,rw in zip(r,sw):
    acc=np.zeros((4,4,4),complex)
    for a1,wa in zip(t1,w1):
        for a2,wb in zip(t2,w2):
            l=np.stack([np.full(8,ri*np.cos(a1)),np.full(8,ri*np.sin(a1)*np.cos(a2)),ri*np.sin(a1)*np.sin(a2)*np.cos(ph),ri*np.sin(a1)*np.sin(a2)*np.sin(ph)],axis=1)
            Sa=S(pp[None,:]-l); Sb=S(p[None,:]-l)
            meas=ri**3*np.sin(a1)**2*np.sin(a2)*wa*wb*wph*Dreg(ri*ri)*ri*rw   # extra ri from dr = r ds
            for nu in range(4):
                M=sum(G[al]@Sa@G[nu]@Sb@G[al] for al in range(4))              # (8,4,4)
                acc[nu]+=np.sum(meas[:,None,None]*M,axis=0)
    tot+=acc
direct=-(1/(4*np.pi**2))*tot
nodes=T.inner_nodes(16,32)
fp=T.vertex_contracted(pp,p,0.0,0.0,nodes,None)
print("time %.0fs"%(time.time()-t0))
for nu in range(4):
    print("nu",nu,"max|direct - FeynmanParam| =",np.max(np.abs(direct[nu]-fp[nu])).round(7)," |Lambda_nu|max=",np.max(np.abs(fp[nu])).round(5))
