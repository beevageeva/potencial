from scipy import integrate
import getopt
import numpy as np
import matplotlib.pyplot as plt
from math import e,pi, sqrt
from matplotlib.widgets import Slider, Button
import ast, inspect
import sys,dis, os.path

G = 6.6 * 10 **(-11)

DEFAULT_NUMBER_POINTS = 1000


def usage():
        print("Usage: %s [--rhodef=<filename>] [--type=p|v|m|mi] [--test=<filename>]\n Option --rhodef is the filename  from where the definition of rho is read(Default: rhodef).\n Option --type can have values : p(potential), v(velocity), m(mass), dp(integrated mass) and tells what to plot.\n Option --test if present tells to use this file (it must be a python file, the extension .py  in this option must be omitted) that MUST be in the same directory to plot the corresponding function calculated analytically\n Option --plotd plot the distribution function in another subplot" % sys.argv[0])



try:
				opts, args = getopt.getopt(sys.argv[1:], "", ["help",  "rhodef=", "type=", "test=", "plotd"])
except getopt.GetoptError as err:
				# print help information and exit:
				print(str(err)) # will print something like "option -a not recognized"
				usage()
				sys.exit(2)
ptype = "p" 
rhodefFilename = "rhodef"
testFilename = None
plotd = False
for o, a in opts:
				if o in("--type"):
								ptype = a
				if o in("--plotd"):
								plotd = True
				elif o in("--test"):
								testFilename = a
				elif o in ("--help"):
								usage()
								sys.exit()
				elif o in ("--rhodef"):
								rhodefFilename = a
				else:
								print("option %s not recognized " % o)
print("rhodefFile=%s, testFilename=%s, plotType=%s, plotd=%s" % (rhodefFilename, testFilename, ptype, str(plotd)))				
if testFilename: 
	if not os.path.isfile(testFilename + ".py") or os.path.dirname(testFilename + ".py")!="" :
		print("testFilename %s is defined but %s.py does not exist OR IT is NOT in thye same directory" % (testFilename, testFilename))
		sys.exit(0)
	exec("import %s as testplot"%testFilename)

K = None
calculateFunction = None
label = None
testFunction = None

maxX = 2


if ptype == 'p':

	label = "Potencial(J/kg)"
	K = [0,0]	

	#r here is only 1 element of r array we want to plot
	def calculateFunction(r, rho):	
		if(r==0):
			return 0
		epsilon = 0.0004
		def f(x):
			int1 = integrate.quad(lambda y: y**2 * rho(y), 0, x )
			return (1.0 / x**2) * int1[0]
		return 4* pi * G * integrate.quad(f, epsilon, r)[0] + K[0]/r + K[1]
	if testFilename:
		testFunction = testplot.calculateP

elif ptype == 'v': 
	label = "Vc(m/s)"
	K = [0]	

	def calculateFunction(r, rho):
		if(r==0):
			return 0
		int1 = integrate.quad(lambda y: y**2 * rho(y), 0, r )	
		t1 =  4.0 * pi * G * int1[0]
		pr1 = t1 + K[0]
		if(pr1<0):
			print("vc**2 negative=%4.2f t1=%4.2f, K=%4.2f, return 0"%(pr1, t1,K[0]))
			return 0
		return sqrt( pr1 / r )
	if testFilename:
		testFunction = testplot.calculateV

elif ptype == 'dp': 
	label = "dp(kg/m2)"
	K = []	

	def calculateFunction(s, rho):
		int1 = integrate.quad(lambda y: (y * rho(y))/sqrt(y**2-s**2), s, np.inf )	
		return  2 *  int1[0] 
	if testFilename:
		testFunction = testplot.calculateDp

elif ptype == 'm': 
	label = "m(kg)"
	K = [0]	

	def calculateFunction(r, rho):
		int1 = integrate.quad(lambda y: y**2 * rho(y), 0, r )	
		return  4.0 * pi *  int1[0] + K[0]
	if testFilename:
		testFunction = testplot.calculateM

else:
	usage()
	sys.exit(0)


#AST TREE
restricted_words = ["x", "e", "pi"] #la variable x de lambda function, las constantes  e y pi importadas de math
numvars = {}

rhof = open(rhodefFilename)
rhoexpr = rhof.read().strip(' \t\n\r')
print("Rho expression = "+ rhoexpr)
def parseLambdaString(expr):
	tree = ast.parse(expr)
	w = ast.walk(tree)
	index = 0
	stop = False
	while(not stop):
		try:
			#python2 - 3
			#http://www.diveinto.org/python3/porting-code-to-python-3-with-2to3.html
			if (sys.version_info[0]==2):
				n=w.next()
			else:
				n=w.__next__()
			print("index=" + str(index) + ", type = " + type(n).__name__)
			if((index==0) and (type(n).__name__!='Module')):
				return False
			if((index==1) and (type(n).__name__!='Expr')):
				return False
			if(type(n).__name__== 'Name'):
				print("nametype, id="+n.id)
				if(not n.id in restricted_words):
					print("NUMVAR")
					numvars[n.id] = 1
			index+=1
		except StopIteration:
			print("finish at : index="+str(index)) 
			stop = True
	return True

if not parseLambdaString(rhoexpr):
	sys.exit("rho expression not valid")

newrhoexpr = rhoexpr
for nv in numvars.keys():
	newrhoexpr = newrhoexpr.replace(nv,"1")

exec("rho = lambda x: " + newrhoexpr)

updateGraph = True


nplots = 1
if(testFunction):
	nplots+=1
if(plotd):
	nplots+=1

fig, ax =  plt.subplots(nplots,1,True)

if(not testFunction and not plotd):
	ax = [ax]

for a in ax:
	a.grid(True)
	a.autoscale(True)

ax0 = ax[0]
ax[nplots - 1].set_xlabel("r(m)")
ax0.set_ylabel(label)


plt.subplots_adjust(left=0.15, bottom=0.4)

print("Start")
print("K=" + str(K) +  ",Lambda:")
print (dis.dis(rho))



#r=np.linspace(float(100)/DEFAULT_NUMBER_POINTS,100,DEFAULT_NUMBER_POINTS)
r=np.linspace(0,100,DEFAULT_NUMBER_POINTS)
yvals = []
for relem in r:
	yvals.append(calculateFunction(relem, rho))
l, = ax0.plot(r,yvals, lw=2, color='red')
ax0.autoscale_view(True,True,True)

#PLOT test FUNCTION
if(testFunction):
	ax1 = ax[1]
	yvalstest = []
	for relem in r:
		yvalstest.append(testFunction(numvars,K,relem))
	l2, = ax1.plot(r,yvalstest, lw=2, color='red')
	ax1.set_ylabel(label)
	ax1.autoscale_view(True,True,True)

#PLOT distribution FUNCTION
if(plotd):
	if(testFunction):
		ax2 = ax[2]
	else:
		ax2 = ax[1]
	yvals = []
	for relem in r:
		yvals.append(rho(relem))
	l3, = ax2.plot(r,yvals , lw=2, color='red')
	ax2.set_ylabel("Densidad")
	ax2.autoscale_view(True,True,True)



def setSliderKVal(val):
	for i in range(len(K)):
		K[i] = eval("sliderK" + str(i) +".val")
	print("in function setSliderKVal SET SLIDER VAL: K=" + str(K) + ", updateGraph=" + str(updateGraph)) 
	if(updateGraph):
		updateG()
#TODO? to create a function for every slider?? in order not to set always for all K


ya=0.1
for i in range(len(K)):
	exec("axK"+ str(i) + " = plt.axes([0.25, " + str(ya)+ ", 0.65, 0.03], axisbg='white')")
	exec("sliderK"+ str(i) + " = Slider(axK" + str(i) + ", 'K" + str(i+1) + "', -30.0, 30.0, valinit=0)")
	exec("sliderK"+ str(i) + ".on_changed(setSliderKVal)")
	ya+=0.05





def setSliderNumVarVal(val):
	#global rho
	global updateGraph
	print("in setSliderNUmVarVal updateGraph="+ str(updateGraph))
	oldUpdateGraph = updateGraph
	for nv in numvars.keys():
		#I keep this in the hash because global vars are not working ?
		numvars[nv] = eval("nv_slider"+nv + ".val")
		#print("type of nv eval = " + str(type(eval(nv))))
	if(oldUpdateGraph):
		updateGraph = False
	for i in range(len(K)):
		exec("sliderK" + str(i) +".set_val(0)")
	if(oldUpdateGraph):
		updateGraph = True
		updateG()	


#TODO put in for before
for nv in numvars.keys():
	exec("nv_ax"+nv + " = plt.axes([0.25, " + str(ya)+ ", 0.65, 0.03], axisbg='white')")
	exec("nv_slider"+nv + " = Slider(nv_ax" + nv + ", '" + nv + "', -30.0, 30.0, valinit=1)")
	exec("nv_slider"+nv + ".on_changed(setSliderNumVarVal)")
	ya+=0.05
	

def setSliderMaxX(val):
	global maxX
	maxX = val
	updateG()

axMaxx = plt.axes([0.25,ya,0.65, 0.03], axisbg='white')
sliderMaxx = Slider(axMaxx, "MaxX(log)", 1, 30, valinit=2)
sliderMaxx.on_changed(setSliderMaxX)




def updateG():
	#print("updateG numvars")
	#print(numvars)
	newrhoexpr = rhoexpr
	for nv in numvars.keys():
		newrhoexpr = newrhoexpr.replace(nv,str(numvars[nv]))
	#print("updateG newrhoexpr")
	#print(newrhoexpr)
	#TODO PYTHON3 exec for lamda assignment is NOT working
	#exec("rho = lambda x: " + newrhoexpr, {"rho": rho}, None)
	rho = eval("lambda x: %s" % newrhoexpr)
	print("CALL updateG maxX = %d" % maxX)
	print("K=" + str(K) +  ",Rho:")
	print (dis.dis(rho))
	txt = plt.text(-5.5, 0.5, "recalculating")
	fig.canvas.draw()
	yvals = []
	r=np.linspace(float(10**maxX)/DEFAULT_NUMBER_POINTS ,int(10**maxX),DEFAULT_NUMBER_POINTS)
	#print("xvals")
	#print(r)
	for relem in r:
		yvals.append(calculateFunction(relem, rho))
	l.set_xdata(r)
	l.set_ydata(yvals)
	#print("yvals")
	#print(len(yvals))
	#TODO why necessary
	ax0.set_xlim(0,int(10**maxX))
	ax0.relim()
	ax0.autoscale_view(True,True,True)
	#replot test function
	if(testFunction):
		yvalstest = []
		for relem in r:
			yvalstest.append(testFunction(numvars,K,relem))
		l2.set_xdata(r)
		l2.set_ydata(yvalstest)
		#TODO why necessary
		ax1.set_xlim(0,int(10**maxX))
		ax1.relim()
		ax1.autoscale_view(True,True,True)
	#replot distr function TODO DO not replot it if K are changed
	if(plotd):
		yvalstest = []
		for relem in r:
			yvalstest.append(rho(relem))
		l3.set_xdata(r)
		l3.set_ydata(yvalstest)
		#TODO why necessary
		ax2.set_xlim(0,int(10**maxX))
		ax2.relim()
		ax2.autoscale_view(True,True,True)

	txt.remove()
	fig.canvas.draw()







resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color='white', hovercolor='0.975')
def reset(event):
	global updateGraph
	updateGraph = False
	for i in range(len(K)):
		exec("sliderK" + str(i) +".reset()")
		K[i] = 0
	for nv in numvars.keys():
		exec("nv_slider"+nv + ".reset()")
	updateGraph = True
	updateG()

button.on_clicked(reset)



plt.show(block=True)



	
