#coding=utf-8
import numpy
import scipy.optimize
import scipy.misc
import pylab
from matplotlib.font_manager import FontProperties

Sigma_Water = 0.07119
Temp = 30+273.15
c_to_g = lambda c, d, T:(-c/(8.31* T))*d
Data = sorted([
	(0.00, 741),
	(0.02, 684),
	(0.05, 605),
	(0.10, 524),
	(0.20, 451),
	(0.30, 389),
	(0.40, 345),
	(0.50, 304)

], key = lambda x:x[1], reverse=True)

def func_log(x, a, c, d):
	#a, b, c, d = p
	return a*numpy.log(c+x)+d

def Sigma_C_fitting(Data, Sigma_Water, Freedom=1):
	Sort_Data = sorted(Data, key = lambda x:x[1], reverse=True)
	#print Sort_Data
	if Sort_Data[0][0] != 0:
		raise ValueError("I need a data for pure water")
	
	D = numpy.array(Sort_Data)
	
	x, y = D[:, 0], map(lambda x:Sigma_Water*x/D[0][1], D[:, 1])
	#z = numpy.poly1d(numpy.polyfit(x, y, Freedom))
	#print z
	print "Sigma_list is:\n"
	print "\t".join(map(lambda x:"%.4f"%x, list(y)))
	print ""
	popt, pcov = scipy.optimize.curve_fit(func_log, x, y, maxfev = 65536)
	z = lambda x:func_log(x, *popt)
	#print x
	#print y, map(z, x)
	print "Func is y = %s*ln(%s+x)+%s\n"%tuple(popt)
	
	return (x, y), z
	

H =  numpy.array(Data)	
D, Func_Sigma_C = Sigma_C_fitting(Data, Sigma_Water, 3)
c_and_g= numpy.array((list(H[:, 0]), map(lambda x:scipy.misc.derivative(Func_Sigma_C, x0 = x, dx=1e-10),  H[:, 0]))).T



c_div_g_and_c = numpy.array([(i[0], i[0]/c_to_g(i[0], i[1], Temp)) for i in list(c_and_g) if i[0] != 0.0])

print "d\sigma/dc is:\n"
print "\t".join(map(lambda x:"%.4f"%(x), list(c_and_g[:, 1])))
print ""

print "Gamma is(times 10^6):"
print "\t".join(map(lambda x, y:"%.4f"%(x/y*1e6 if y!=0 else 0), [0, ]+list(c_div_g_and_c[:, 0]), [0, ]+list(c_div_g_and_c[:, 1])))
print ""

print "c/\Gamma is:\n"
print "\t".join(map(lambda x:"%.4f"%(x), [0, ]+list(c_div_g_and_c[:, 1])))
print ""

Tmp = numpy.polyfit(c_div_g_and_c[:, 0], c_div_g_and_c[:, 1], 1)
print "1/g_{inf} = ", Tmp[0]
z = numpy.poly1d(Tmp)
#print c_div_g_and_c[:, 1]
#print map(z, c_div_g_and_c[:, 0])


import matplotlib.pyplot as plt
plt.figure(0)
plt.plot(D[0], D[1], "o")
plt.plot(numpy.linspace(-0.01, 0.55), map(Func_Sigma_C, numpy.linspace(-0.01, 0.55)))
R2 = lambda y, y_guess: 1-(((y-y_guess)**2).sum())/(((y-y.mean())**2).sum())
print "R^2 is ", R2(numpy.array(map(Func_Sigma_C, D[0])), numpy.array(D[1]))
CF = FontProperties(fname='C:\\Windows\\Fonts\\msyh.ttc')
plt.text(0.4, 0.06, "$R^{2}=%.4f$"%R2(numpy.array(map(Func_Sigma_C, D[0])), numpy.array(D[1])), fontsize=20)
plt.xlabel(u"浓度$(mol\cdot L^{-1})$", fontproperties=CF, fontsize=20)
plt.ylabel(u"表面张力$(N\cdot m^{-1})$", fontproperties=CF, fontsize=20)
plt.title(u"$\sigma -c$图", fontproperties=CF, fontsize=24)

plt.figure(1)
plt.plot(c_div_g_and_c[:, 0], c_div_g_and_c[:, 1], "o")
plt.plot(c_div_g_and_c[:, 0], map(z, c_div_g_and_c[:, 0]))
plt.xlabel(u"浓度$(mol\cdot L^{-1})$", fontproperties=CF, fontsize=20)
plt.ylabel(u"$\\frac{c}{\Gamma}(N\cdot mol\cdot m^{-3})$", fontproperties=CF, fontsize=20)
plt.title(u"$\\frac{c}{\Gamma} -c$图", fontproperties=CF, fontsize=30)

plt.show()