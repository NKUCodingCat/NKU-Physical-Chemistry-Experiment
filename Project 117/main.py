#coding=utf-8
Data = [
	(  100,  1.2036),  #mA, V
	(   50,  1.0599),
	(   20,  0.9003),
	(   10,  0.8013),
	(    5,  0.7366),
	(    2,  0.6883),
	(    1,  0.6694),
	(  0.5,  0.6437),
	(  0.2,  0.6158),
	(  0.1,  0.5960),
	( 0.05,  0.5843),  
	( 0.02,  0.5742),
	( 0.01,  0.5682),
	(0.005,  0.5606),
	(0.002,  0.5581)
]

import numpy
from scipy.optimize import curve_fit
from scipy import optimize
import scipy
import random

def Fun(x, a1, b1, a2, b2):
	###
	# Function is y=(a1*x+b1)*g(x) + (a2*x+b2)*(1-g(x)) where g(x) = 1/(1+e^(point))
	###

	Point = (b2-b1)/(a1-a2) if a1-a2 else -12

	return numpy.array([(a1*i+b1 if i < Point else (a2*i+b2)) for i in x])

class F:
	def __init__(self, xdata, ydata):
		self.xdata = xdata
		self.ydata = ydata
		self.R2 = lambda y, y_guess: 1-(((y-y_guess)**2).sum())/(((y-y.mean())**2).sum())
	
	def Func(self, *x):
		self.args = x[0]
		Point = (self.args[3]-self.args[1])/(self.args[0]-self.args[2])
		T = Fun(self.xdata, *self.args)
		return numpy.sum((ydata - T)**2)


D = numpy.array(Data)
xdata, ydata = numpy.array(map(lambda x:numpy.log(x*2/1000), D[:, 0])), D[:, 1]

MM = F(xdata, ydata)


from sklearn import linear_model
V1 = scipy.stats.linregress(xdata[:-3], ydata[:-3])
V2 = scipy.stats.linregress(xdata[:3], ydata[:3])
# 此处猜测线性关系大致a1, a2, b2, b1

popt = optimize.minimize(MM.Func, V1[:2]+V2[:2],  method='SLSQP').x#, options = {'maxiter': 15000, 'ftol':1e-20, "xtol":1e-20}).x
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

y_pop = Fun(xdata, *popt)
print popt
print y_pop
print "\n".join(map(lambda x: '%.3f'%x, list(xdata)))


Point = (popt[3]-popt[1])/(popt[0]-popt[2])
plt.plot(xdata, ydata, 'o')

plt.plot(numpy.array(sorted(list(xdata)+[Point, ])), Fun(numpy.array(sorted(list(xdata)+[Point, ])), *popt))
CF = FontProperties(fname = 'C:\\Windows\\Fonts\\msyh.ttc')
plt.text(xdata[-5], ydata[-5]+0.5*(ydata[0]-ydata[-1]), "$R^2 = %.4f$"%(MM.R2(ydata, y_pop)), fontsize=18)
plt.text(xdata[-3], ydata[-3]+0.1*(ydata[0]-ydata[-1]), r"$\varphi = %.4f+%.4f\,log\,i$"%(popt[1], popt[0]),  fontsize=18)
plt.text(xdata[2], ydata[2]-0.1*(ydata[0]-ydata[-1]), r"$\varphi = %.4f+%.4f\,log\,i$"%(popt[3], popt[2]),  fontsize=18)
plt.xlabel(r"电流密度对数值 $log\,i$($A/cm^2$)".decode('utf-8'), fontproperties=CF, fontsize=20)
plt.ylabel(r"电位差 $\varphi$(V, vs. SCE)".decode('utf-8'), fontproperties=CF, fontsize=20)
plt.title(r"$\varphi$ - $log\,i$ 图".decode('utf-8'), fontproperties=CF, fontsize=24)
plt.show()