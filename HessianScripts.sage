def colorHessGraph(p):
	"""Plots the Hessian graph over Fp, with the color coding described in Appendix E of https://arxiv.org/abs/2407.17042."""
	import sage.graphs.graph_plot
	F0=GF(p)
	F=GF(p^6)
	T.<x>=PolynomialRing(F)
	E=EllipticCurve(F0,[0,-1728])
	uu=F.multiplicative_generator()
	sexRoot=uu^((p^6-1)/(6*(p-1)))
	if mod(p,3)==1:
		exps=[0,3,4,1,2,5] #used to sort the twists to match the color coding
	else:
		exps=[0,3]
	tws=[EllipticCurve(F0,[0,-1728*sexRoot^(6*ell)]) for ell in exps]
	tws = tws[1:] + tws[:1]
	twsMaps=[EE.change_ring(F).isomorphism_to(E.change_ring(F)) for EE in tws]
	S=[[twsMaps[i](P) for P in tws[i].points() if P!=tws[i](0)] for i in [0..len(tws)-1]]
	G = DiGraph(loops=1, multiedges=0)
	prevertices=[s  for SS in S for s in SS]
	sqrt3=(x^2+3).any_root()
	k=-6912
	E=E.change_ring(F)
	T1=E(0,2^3*3*sqrt3,1)
	for P in prevertices:
		if ((P[2] !=0) and (P != T1) and (P != -T1)):
			newP=(-(P[0]^3 +k) / (3 * P[0]^2), - P[1] * (P[0]^3 - 2 * k) / (3 * sqrt3 * P[0]^3),1)
			newP=E(newP)
			G.add_edge(F0(P[0]^3),F0(newP[0]^3))
		else:
			G.add_edge(F0(P[0]^3),oo)
		G.add_edge(oo,oo)

	print('Showing graph...')
	# The following labels replace "Infinity" with ∞ symbol
	labels = {v: (r"$\infty$" if v == Infinity else v) for v in G.vertices()}
	if len(S)==2:
		color=['lightgray', 'white']
	else:
		color=['lightgray','lightblue','steelblue','lightgreen','green', 'white']
	colors = {color[i] :[F0(P[0]^3) for P in prevertices if P in S[i]] for i in [0..len(S)-1]}
	colors['white'].append(oo)
	sage.graphs.graph_plot.DEFAULT_SHOW_OPTIONS['figsize'] = [100,100]
	G.show(loop_size=0.1, vertex_size=100, vertex_colors=colors, vertex_labels=labels, iterations=20000, title="Hess p=%d" % p)
	#GG = G.plot(loop_size=0.1, vertex_size=700, dpi=300)
	#GG.save("Hessian graph for p=%d.pdf" % p)
	sage.graphs.graph_plot.DEFAULT_SHOW_OPTIONS['figsize'] = [10,10]
	return G
	
def genFinder(F):
	"""Finds a generator of F*/F*^6"""
	q=F.cardinality()
	S.<X>=PolynomialRing(F)
	u=F.random_element()
	while u.is_square()==True or (mod(q,3)==1 and (X^3-u).is_irreducible()==False):
		u=F.random_element()
	return u
	
def iterHess(p,j,N):
	"""
	Computes the N-th iterated Hessian of j in Fp.
	"""
	F=GF(p,'a')
	S.<X>=PolynomialRing(F)
	u=genFinder(F)
	s=0
	while (X^3-j*u^(2*s)).is_irreducible():
		s+=1
	t=0
	if (u^(2*s)*(j-1728)).is_square()==False:
		t=1
	u=u^(2*s+3*t)
	x=(X^3-u*j).any_root()
	y=(X^2-u*(j-1728)).any_root()
	E=EllipticCurve([0,-u*1728])
	card=E.cardinality()
	N2=power_mod(-3, floor(N/2), card)
	(x,y,z)=N2*E([x,y,1])
	if z==0:
		return oo
	else:
		j=u^(-1)*x^3
		if mod(N,2)==1:
			if j==0:
				return oo
			else:
				return (6912-j)^3/(27*j^2)
		else:
			return j
	
