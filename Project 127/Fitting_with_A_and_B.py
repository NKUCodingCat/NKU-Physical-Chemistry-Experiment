#coding=utf-8
import numpy
import scipy.optimize
import scipy.misc
import pylab
from matplotlib.font_manager import FontProperties

Solvent_Time = 79.5
Data = sorted([
	(0.006, 130.8),
	(0.006*12/15.0, 117.853),
	(0.006*12/20.0, 106.813),
	(0.006*12/25.0, 100.787),
	(0.006*12/30.0, 96.98)
], key = lambda x:x[1], reverse=True)

R2 = lambda y, y_guess: (lambda y, y_guess: 1-(((y-y_guess)**2).sum())/(((y-y.mean())**2).sum()))(*map(numpy.array, (y, y_guess)))
_ = lambda a:map(str, list(a))

class F_class:
	def __init__(self, c, t):
		self.c = c
		self.t = t
		self.CanCha = lambda x, y, func: ((numpy.array(map(func, x))-numpy.array(y))**2).sum()
		self.Sta = [1e2, 0.1] #A, B
		self.eta_r = lambda t, A, B: (A*t - B/t)/(A*Solvent_Time - B/Solvent_Time)
		self.eta_sp = lambda t, A, B:(self.eta_r(t, A, B) - 1)

	def Func(self, *arg):
		A, B = arg[0]
		if abs(A/B)< 0.01:
			return 1e200
		self.y1 = map(lambda c, t:(self.eta_sp(t, A, B)/c), self.c, self.t)
		self.y2 = map(lambda c, t:(numpy.log(self.eta_r(t, A, B))/c), self.c, self.t)
		self.z1, self.z2 = numpy.polyfit(self.c, self.y1, 1), numpy.polyfit(self.c, self.y2, 1)
		self.zt1, self.zt2 = map(numpy.poly1d, [self.z1, self.z2])
		return self.CanCha(self.c, self.y1, self.zt1) + self.CanCha(self.c, self.y2, self.zt2)


K = numpy.array(Data)
F = F_class(K[:, 0], K[:, 1])	
RES = scipy.optimize.minimize(F.Func, F.Sta, method='Nelder-Mead', options = {'maxiter': 15000, 'ftol':1e-10, "xtol":1e-10})
A = RES.x
print A
import matplotlib.pyplot as plt
c = list(K[0])
F.c2 = list(F.c)
plt.figure(0)
plt.plot(F.c2, F.y1, "o")
plt.plot([0, ]+F.c2, map(F.zt1, [0, ]+F.c2), label="$\\frac{\eta_{sp}}{c}$")
plt.plot(F.c2, F.y2, "o")
plt.plot([0, ]+F.c2, map(F.zt2, [0, ]+F.c2), "b--", label="$\\frac{ln(\eta_{r})}{c}$")
plt.legend(loc=7, bbox_transform=plt.gcf().transFigure, fontsize=24)
CF = FontProperties(fname='C:\\Windows\\Fonts\\msyh.ttc')

plt.text(F.c2[3], F.y1[3]+2.5, "$R^{2}=%.4f$"%R2(map(F.zt1, F.c), F.y1), fontsize=20)
plt.text(F.c2[3], F.y2[3]+2.5, "$R^{2}=%.4f$"%R2(map(F.zt2, F.c), F.y2), fontsize=20)
plt.text(F.c2[3], F.y1[3]+5, "$y_{1}=%.2fx+%.2f$"%(tuple(F.z1)), fontsize=20)
plt.text(F.c2[3], F.y2[3]+5, "$y_{2}=%.2fx+%.2f$"%(tuple(F.z2)), fontsize=20)

plt.xlabel(u"浓度$(g/mL)$", fontproperties=CF, fontsize=20)
plt.ylabel(u"$\\frac{\eta_{sp}}{c}$ or $\\frac{ln(\eta_{r})}{c}$", fontproperties=CF, fontsize=20)
plt.title(u"$\\frac{\eta_{sp}}{c}/\\frac{ln(\eta_{r})}{c}-c$图", fontproperties=CF, fontsize=24)
plt.show()