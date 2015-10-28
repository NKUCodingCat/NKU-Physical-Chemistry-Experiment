#coding=utf-8
import math
import numpy
import scipy
from scipy import optimize
from sklearn import linear_model

Data = [
	(0.002, 157.33),
	(0.004, 299.33),
	(0.006, 409.33),
	(0.007, 455.67),
	(0.008, 507.67),
	(0.009, 562.33),
	(0.010, 588.00),
	(0.012, 668.67),
	(0.014, 744.67),
	(0.016, 815.00),
	(0.018, 893.00),
	(0.020, 961.00),
]


def linear_Find(RAW_Data):
	G = []
	Data = numpy.array(sorted(RAW_Data, key = lambda x:x[0]))
	for i in range(3, len(RAW_Data)-3+1):
		V1 = scipy.stats.linregress(Data[:i, 0], Data[:i, 1])
		V2 = scipy.stats.linregress(Data[i:, 0], Data[i:, 1])
		#print V1[2]**2, V2[2]**2
		G.append((i, V1, V2, (V1[2]*V2[2])**2))
	return RAW_Data,G

class F:
	def __init__(self, Buff, Start):
		self.Buff = numpy.array(Buff)
		#self.Prev = (Buff[:Split, 0], Buff[:Split, 1])
		#self.Form = (Buff[Split:, 0], Buff[Split:, 1])
		self.args = (Start[0][0], Start[0][1], Start[1][0], Start[1][1], 10000) #a1, b1, a2, b2, Sigma
		self.res = 0
		self.cmc = 0
	
	def f(self, *x):
		#print self.args
		self.args = x[0]
		G = self.R_Squard(self.Buff, self.args)
		self.res = G
		return G
	
	def R_Squard(self, B, A):
		
		Const_x_CMC = float(A[3] - A[1])/float(A[0] - A[2])  #CMC = (b2-b1)/(a1-a2)
		g = lambda Sigma, delta:(1.0/(1.0+numpy.exp(float(Sigma)*float(delta))))
		Guess = lambda x: (A[0]*x+A[1])*g(A[4], x-Const_x_CMC) + (A[2]*x+A[3])*(1-g(A[4], x-Const_x_CMC))
		X, Y, Y_Guess = B[:, 0], B[:, 1], map(Guess, B[:, 0])
		R = lambda X_ser, Y_ser, Y_Guess: numpy.sum(numpy.array((Y_Guess-Y_ser)**2, dtype = numpy.float64))
		self.cmc = Const_x_CMC
		return R(X, Y, Y_Guess)
	
#G, K = numpy.array(sorted(Data, key = lambda x:x[0])), [(0, [60000,20,0,0,0], [30000,200,0,0,0], 0),]#linear_Find(Data)
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
G,K = linear_Find(Data)
for i in K:
	Func = F(G, i[1:3])
	plt.figure(i[0])
	
	RES = optimize.minimize(Func.f, Func.args, method='Nelder-Mead')#, options = {'maxiter': 15000, "eps":1e-15, 'ftol':1e-15, "xtol":1e-2})
	A = RES.x
	
	Const_x_CMC = float(A[3] - A[1])/float(A[0] - A[2])  #CMC = (b2-b1)/(a1-a2)
	CMC_in_exp = (lambda x:map(float, __import__("re").split("E", "%.2E"%x)))(Const_x_CMC)
	g = lambda Sigma, delta:(1.0/(1+numpy.exp(Sigma*delta)))
	Guess = lambda x: (A[0]*x+A[1])*g(A[4], x-Const_x_CMC) + (A[2]*x+A[3])*(1-g(A[4], x-Const_x_CMC))
	
	

	L = numpy.linspace(0, 0.025, 100)
	Gp = numpy.array(G)

	R2 = lambda y, y_guess: 1-(((y-y_guess)**2).sum())/(((y-y.mean())**2).sum())
	print "R^2 is ", R2(Gp[:, 1], numpy.array(map(Guess, Gp[:, 0])))

	plt.plot(Gp[:, 0], Gp[:, 1], "o")
	plt.plot(L, numpy.array(map(Guess, L)))
	Text = r"$y=(%d x+%d)\left(\frac{1}{1+e^{%d \left(x-%.2f\times10^{%d}\right)}}\right)+(%d x+%d) \left(1-\frac{1}{1+e^{%d \left(x-%.2f\times10^{%d}\right)}}\right)$"%(RES.x[0], RES.x[1], RES.x[-1], CMC_in_exp[0], CMC_in_exp[1], RES.x[2], RES.x[3], RES.x[-1], CMC_in_exp[0], CMC_in_exp[1])
	Tp, Tf, cmc_point = r"$y=%d x+%d$"%(RES.x[0], RES.x[1]), r"$y=%d x+%d$"%(RES.x[2], RES.x[3]), r"(%.6f, %.3f)"%(Const_x_CMC, Guess(Const_x_CMC)) 
	plt.plot(L[:40], numpy.array(map(lambda x:RES.x[0]*x+RES.x[1], L[:40])))
	plt.plot(L[20:], numpy.array(map(lambda x:RES.x[2]*x+RES.x[3], L[20:])))
	#plt.text(0.005, 200, Text, fontsize=24)
	plt.text(0.0005, 400, Tp, fontsize=18)
	plt.text(0.012, 600, Tf, fontsize=18)
	plt.text(0.015, 200,"$R^2=%.4f$"%R2(Gp[:, 1], numpy.array(map(Guess, Gp[:, 0]))) ,fontsize=18)
	#print cmc_point
	#plt.text(0.005, 700, cmc_point)
	print  Const_x_CMC, RES.fun
	CF = FontProperties(fname='C:\\Windows\\Fonts\\msyh.ttc')
	plt.xlabel(u"浓度$(mol\cdot L^{-1})$", fontproperties=CF, fontsize=20)
	plt.ylabel(u"电导率$(\mu S\cdot cm^{-1})$", fontproperties=CF, fontsize=20)
	#break
	
plt.show()
	
	