#for the potential there will be two integration constants K[0] and K[1]
#for the rest(velocity , masss, integrated mass) there is only one K[0]
#numvars is the hash containing variable parameters
def calculateP(numvars,K,x):
	return numvars["A"] * x + numvars["B"] + K[0] - K[1]

def calculateV(numvars, K , x):
	return x * numvars["A"]


def calculateM(numvars, K , x):
	return x * numvars["A"]

def calculateDp(numvars, K , x):
	return numvars["A"] + numvars["B"]

