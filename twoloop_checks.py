from twoloop_lib import *
# 1. gamma algebra
for a in range(4):
    for b in range(4):
        assert np.allclose(G[a]@G[b]+G[b]@G[a],2*(a==b)*I4)
print("Clifford algebra ok")
# 2. identities used: g_a (3 gammas) g_a = -2 reversed ; g_a X g_a for 2 gammas
rng=np.random.default_rng(1); a,b=rng.normal(size=4),rng.normal(size=4)
nu=2
lhs=sum(G[al]@(slash(a)@G[nu]@slash(b))@G[al] for al in range(4)); rhs=-2*slash(b)@G[nu]@slash(a)
print("3-gamma identity:",np.allclose(lhs,rhs))
lhs=sum(G[al]@(slash(a)@G[nu])@G[al] for al in range(4)); print("2-gamma identity:",np.allclose(lhs,4*a[nu]*I4))
# 3. Ward identity: Lambda_nu(p,p) = -i dSigma/dp_nu  (PV photon, mu=0.3, Lam=8)
mu,Lam=0.3,8.0; nodes=inner_nodes(16,32)
p=np.array([0.4,-0.7,0.5,0.9])
V=vertex_contracted(p,p,mu,Lam,nodes,None)
def Sig(p):
    p2=np.sum(p*p); sA,sB=sigma_parts(p2,mu,Lam,n=128); return m*sA[0]*I4+1j*sB[0]*slash(p)
h=1e-4
for nu in range(4):
    e=np.zeros(4); e[nu]=h
    d=(Sig(p+e)-Sig(p-e))/(2*h)
    print("nu",nu,"WI max|Lambda_nu(p,p) + i dSigma/dp_nu| =",np.max(np.abs(V[nu]+1j*d)).round(8)," |Lambda|=",np.max(np.abs(V[nu])).round(5))
