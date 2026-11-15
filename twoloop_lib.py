# Euclidean QED, m=1, g=alpha/pi.  Gauge: Feynman.  Photon: Pauli-Villars (mass mu minus mass Lam) on the INNER photon only.
import numpy as np
from scipy.special import exp1
m=1.0
REG=None
def set_reg(kind,Lam):
    global REG
    if kind=='single': REG=[(1.0,0.0),(-1.0,Lam)]
    elif kind=='double': REG=[(1.0,0.0),(-4/3,Lam),(1/3,2*Lam)]
    elif kind=='gauss': REG=('gauss',1.0/Lam**2)   # photon propagator e^{-l^2/Lam^2}/l^2  (proper-time cutoff s>=1/Lam^2 on the photon only)
    else: raise ValueError
    return REG
s0=np.eye(2); sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)
Z2=np.zeros((2,2),complex)
G=[np.block([[Z2,-1j*s],[1j*s,Z2]]) for s in (sx,sy,sz)]+[np.block([[Z2,s0],[s0,Z2]])]
G=np.array(G); I4=np.eye(4,dtype=complex)
def slash(v): return np.einsum('...m,mab->...ab',v,G)
def S(p): # (-i pslash + m)/(p^2+m^2)
    p2=np.sum(p*p,axis=-1)
    return (-1j*slash(p)+m*I4)/(p2+m*m)[...,None,None]
def gl(n,a=0.,b=1.):
    x,w=np.polynomial.legendre.leggauss(n); return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w
# ---- self energy, PV photon.  Sigma(p)=g/4 * Int dx [2i(1-x) pslash + 4m] ln(D_Lam/D_mu),  D_M = x(1-x)p^2 + x m^2 + (1-x) M^2
def sigma_parts(p2,mu,Lam,n=64):
    # returns (sA,sB) with Sigma = m*sA + i*pslash*sB ; p2 may be negative (continuation) as long as D>0
    x,w=gl(n)
    X=x[:,None]; P=np.atleast_1d(p2)[None,:]
    Dm=lambda M: X*(1-X)*P+X*m*m+(1-X)*M*M
    if isinstance(REG,tuple): L=exp1(REG[1]*(X*(1-X)*P+X*m*m)/(1-X))
    else: L=-sum(c*np.log(Dm(M)) for c,M in REG)
    sA=0.25*np.sum(w[:,None]*4*L,axis=0); sB=0.25*np.sum(w[:,None]*2*(1-X)*L,axis=0)
    return sA,sB
def sigma_os(mu,Lam): # on-shell value (pslash -> i m  => i*pslash = -m): m*(sA - sB) at p^2=-m^2
    sA,sB=sigma_parts(-m*m,mu,Lam,n=96); return m*(sA[0]-sB[0])
def Sigma_tilde(k,mu,Lam,sos): # matrix, k (...,4);  mass counterterm only
    p2=np.sum(k*k,axis=-1); sA,sB=sigma_parts(p2.ravel(),mu,Lam); sA=sA.reshape(p2.shape); sB=sB.reshape(p2.shape)
    return (m*sA-sos)[...,None,None]*I4+1j*sB[...,None,None]*slash(k)
# ---- vertex (one-loop), p' outgoing, p incoming:  Lambda_nu = g/4 Int dxdy { -N0_nu (1/D_mu - 1/D_Lam) + 2 gamma_nu ln(D_Lam/D_mu) }
#  N0_nu = 2 bslash g_nu aslash - 4 i m (a_nu + b_nu) - 2 m^2 g_nu ;  a=p'-c, b=p-c, c=x p + y p' ;  D_M = x p^2 + y p'^2 - c^2 + (x+y) m^2 + z M^2
def inner_nodes(nu=16,nv=24):
    u,wu=gl(nu); v,wv=gl(nv); z=v**4; wz=4*v**3*wv
    Z_,U_=np.meshgrid(z,u,indexing='ij'); W=(wz[:,None]*(1-Z_))*wu[None,:]
    return Z_.ravel(),U_.ravel(),W.ravel()
def inner_nodes_log(nu=14,npan=10,nper=6,wmax=36.0):
    u,wu=gl(nu); ws=[];wws=[]
    for a,b in zip(np.linspace(0,wmax,npan+1)[:-1],np.linspace(0,wmax,npan+1)[1:]):
        x,w=gl(nper,a,b); ws+=list(x); wws+=list(w)
    ws=np.array(ws); wws=np.array(wws); z=np.exp(-ws); wz=z*wws
    Z_,U_=np.meshgrid(z,u,indexing='ij'); W=(wz[:,None]*(1-Z_))*wu[None,:]
    return Z_.ravel(),U_.ravel(),W.ravel()
def vertex_contracted(pp,p,mu,Lam,nodes,stream):
    """returns Lambda_nu matrices for nu=0..3 at momenta (pp,p) (single points).  shape (4,4,4)"""
    z,u,w=nodes; x=(1-z)*u; y=(1-z)*(1-u)
    c=x[:,None]*p+y[:,None]*pp; a=pp-c; b=p-c
    D0=x*np.dot(p,p)+y*np.dot(pp,pp)-np.sum(c*c,axis=1)+(x+y)*m*m
    Dmu=D0+z*mu*mu; DL=D0+z*Lam*Lam
    sa=slash(a); sb=slash(b)
    if isinstance(REG,tuple):
        A=REG[1]; W1=np.exp(-A*D0/z)/D0; W2=exp1(A*D0/z)
    else:
        W1=sum(c/(D0+z*M*M) for c,M in REG); W2=-sum(c*np.log(D0+z*M*M) for c,M in REG)
    out=np.zeros((4,4,4),complex)
    for nu in range(4):
        N0=2*sb@G[nu]@sa-4j*m*(a[:,nu]+b[:,nu])[:,None,None]*I4-2*m*m*G[nu]
        t1=np.sum((w*(-W1))[:,None,None]*N0,axis=0)
        t2=G[nu]*np.sum(w*2*W2)
        out[nu]=0.25*(t1+t2)
    return out
