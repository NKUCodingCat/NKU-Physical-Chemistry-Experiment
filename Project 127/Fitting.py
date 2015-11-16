#coding=utf-8
import numpy
import scipy.optimize
import scipy.misc
import pylab
from matplotlib.font_manager import FontProperties

Solvent_Time = 104.43
Data = sorted([
	#(0.006, 174.00),
	(0.006*12/15.0, 154.66),
	(0.006*12/20.0, 140.06),
	(0.006*12/25.0, 132.16),
	(0.006*12/30.0, 126.69)
], key = lambda x:x[1], reverse=True)

R2 = lambda y, y_guess: (lambda y, y_guess: 1-(((y-y_guess)**2).sum())/(((y-y.mean())**2).sum()))(*map(numpy.array, (y, y_guess)))
_ = lambda a:map(str, list(a))

def Data_Fitting(Data, Solvent_Time):
	Sort_Data = sorted(Data, key = lambda x:x[1], reverse=True)
	D = numpy.array(Sort_Data)
	c, eta_r = D[:, 0], map(lambda y: y/Solvent_Time, D[:, 1])
	eta_sp = map(lambda y: y-1, eta_r)
	eta_sp_div_c, ln_eta_r_div_c = map(lambda x, y:y/x, c, eta_sp), map(lambda x, y:numpy.log(y)/x, c, eta_r)
	
	print "eta_r list is: \n", "\t".join(_(eta_r))
	print "eta_sp list is: \n", "\t".join(_(eta_sp))
	print "eta_sp/c list is: \n", "\t".join(_(eta_sp_div_c))
	print "ln(eta_r)/c list is: \n", "\t".join(_(ln_eta_r_div_c))
	
	return c, eta_sp_div_c, ln_eta_r_div_c


	
	
K = Data_Fitting(Data, Solvent_Time)
L1 = numpy.polyfit(K[0], K[1], 1)
L2 = numpy.polyfit(K[0], K[2], 1)
import os
with open(os.path.split(os.path.realpath(__file__))[0]+"/data1.dat", "w") as f:
	f.write("\n".join(map(lambda x:"%s\t%s"%(x[0], x[1]), (list(numpy.array([list(K[0])+list(K[0]), K[1]+K[2]]).T)))))
G = [
	["浓度$(g/mL)$", ]+list(K[0]),
	["$\\frac{\eta_{sp}}{c}$", ]+list(K[1]),
	["$\\frac{ln(\eta_{r})}{c}$", ]+list(K[2])
]
with open(os.path.split(os.path.realpath(__file__))[0]+"/data2.dat", "w") as f:
	f.write("\n".join(map(lambda x:"\t".join(map(str, x)), G)))


z1, z2 = map(numpy.poly1d, [L1, L2])
print "==================================\nif use two indepentant line"
print "y1 = %.4f*x+%.4f, y2 = %.4f*x+%.4f"%(L1[0], L1[1], L2[0], L2[1])
print "point is ", (z1(0)+z2(0))/2.0
#===========
class F_class:
	def __init__(self, x, y1, y2):
		self.x = x
		self.y1 = y1
		self.y2 = y2
		self.CanCha = lambda x, y, func: ((numpy.array(map(func, x))-numpy.array(y))**2).sum()
		self.Sta = [1e5, -1e3, 100] #a1, a2, b
	def Func(self, *arg):
		#print arg
		a1, a2, b = arg[0]
		zt1, zt2 = lambda x:a1*x+b, lambda x:a2*x+b
		return self.CanCha(self.x, self.y1, zt1) + self.CanCha(self.x, self.y2, zt2)

F = F_class(*K)	
RES = scipy.optimize.minimize(F.Func, F.Sta, method='Nelder-Mead', options = {'maxiter': 15000, 'ftol':1e-10, "xtol":1e-10})
A = RES.x			
#print numpy.array([K[1],K[2]]).shape
#popt, pcov = scipy.optimize.curve_fit(Anti_Fitting,numpy.array([K[1],K[2]]) , K[0], maxfev = 65536)
zn1, zn2=lambda x:A[0]*x+A[2], lambda x:A[1]*x+A[2]
print "=======================================\nif frocely let two lines cross on y axis"
print "y1 = %.4f*x+%.4f, y2 = %.4f*x+%.4f"%(A[0], A[2], A[1], A[2])
#print popt
#===========

#Cross = (L2[1]-L1[1])/(L1[0]-L2[0])
#print K[2], map(z2, K[0])
#print R2(K[1], map(z1, K[0])), R2(K[2], map(z2, K[0]))

import matplotlib.pyplot as plt
c = list(K[0])
plt.figure(0)
plt.plot(K[0], K[1], "o")
plt.plot([-0.001, ]+c, map(z1, [-0.001, ]+c), label="$\\frac{\eta_{sp}}{c}$")
plt.plot(K[0], K[2], "o")
plt.plot([-0.001, ]+c, map(z2, [-0.001, ]+c), "b--", label="$\\frac{ln(\eta_{r})}{c}$")
plt.legend(loc=7, bbox_transform=plt.gcf().transFigure, fontsize=24)

CF = FontProperties(fname='C:\\Windows\\Fonts\\msyh.ttc')

plt.text(K[0][2], K[1][2]+3, "$R^{2}=%.4f$"%R2(K[1], map(z1, K[0])), fontsize=20)
plt.text(K[0][2], K[2][2]+3, "$R^{2}=%.4f$"%R2(K[2], map(z2, K[0])), fontsize=20)

plt.text(K[0][2], K[1][2]+5, "$y_{1}=%.2fx+%.2f$"%(L1[0], L1[1]), fontsize=20)
plt.text(K[0][2], K[2][2]+5, "$y_{2}=%.2fx+%.2f$"%(L2[0], L2[1]), fontsize=20)
plt.xlabel(u"浓度$(g/mL)$", fontproperties=CF, fontsize=20)
plt.ylabel(u"$\\frac{\eta_{sp}}{c}$ or $\\frac{ln(\eta_{r})}{c}$", fontproperties=CF, fontsize=20)
plt.title(u"$\\frac{\eta_{sp}}{c}/\\frac{ln(\eta_{r})}{c}-c$图", fontproperties=CF, fontsize=24)





#==================================
plt.figure(1)
plt.plot(K[0], K[1], "o")
plt.plot([-0.001, ]+c, map(zn1, [-0.001, ]+c), label="$\\frac{\eta_{sp}}{c}$")
plt.plot(K[0], K[2], "o")
plt.plot([-0.001, ]+c, map(zn2, [-0.001, ]+c), "b--", label="$\\frac{ln(\eta_{r})}{c}$")
plt.legend(loc=7, bbox_transform=plt.gcf().transFigure, fontsize=24)

CF = FontProperties(fname='C:\\Windows\\Fonts\\msyh.ttc')

plt.text(0.003, 105, "$y_{1}=%.2fx+%.2f$"%(A[0], A[2]), fontsize=20)
plt.text(0.003, 75, "$y_{2}=%.2fx+%.2f$"%(A[1], A[2]), fontsize=20)
plt.text(0.003, 100, "$R^{2}=%.4f$"%R2(K[1], map(zn1, K[0])), fontsize=20)
plt.text(0.003, 85, "$R^{2}=%.4f$"%R2(K[2], map(zn2, K[0])), fontsize=20)

plt.xlabel(u"浓度$(g/mL)$", fontproperties=CF, fontsize=20)
plt.ylabel(u"$\\frac{\eta_{sp}}{c}$ or $\\frac{ln(\eta_{r})}{c}$", fontproperties=CF, fontsize=20)
plt.title(u"$\\frac{\eta_{sp}}{c}/\\frac{ln(\eta_{r})}{c}-c$图", fontproperties=CF, fontsize=24)



plt.show()