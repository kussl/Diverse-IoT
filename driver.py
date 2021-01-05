from simulator import Simulator
import os, statistics as stats 
import numpy as np


S = Simulator()

# def exp_stepwise():
# 	OX = [] 
# 	labels = [] 
# 	O,i, IS = S.simulate_stepwise(.9,.1)
# 	OX.append(IS)
# 	labels.append('β=0.9 ω=0.1')

# 	O,i, IS = S.simulate_stepwise(.5,.5)
# 	OX.append(IS)
# 	labels.append('β=0.5 ω=0.5')

# 	O,i, IS = S.simulate_stepwise(.1,.9)
# 	OX.append(IS)
# 	labels.append('β=0.1 ω=0.9')
	
# 	return OX,labels,i

def exp_stepwise():
	OX = [] 
	labels = [] 
	P = .9

	N = 10

	def one_iteration(P,OX,labels,N,a,b):
		ISA = [] 
		i = 0 
		for j in range(N):
			O,i, IS = S.simulate_stepwise(a,b)
			ISA.append(np.array(IS))

		mean = np.mean(ISA, axis=0).tolist()
		OX.append(mean)
		labels.append('β={} ω={}'.format(a,b))
		return i 

	one_iteration(P,OX,labels,N,.9,.1)
	#P = round(P*.5,2)
	i = one_iteration(P,OX,labels,N,.5,.5)
	#P = round(P*.8,2)
	i = one_iteration(P,OX,labels,N,.1,.9)

	return OX,labels,i

def exp_inaccurate():
	OX = [] 
	labels = [] 
	P = .9

	N = 10

	def one_iteration(P,OX,labels,N):
		ISA = [] 
		i = 0 
		for j in range(N):
			O,i, IS = S.simulate_inaccurate(.9,.1,P)
			ISA.append(np.array(IS))

		mean = np.mean(ISA, axis=0).tolist()
		OX.append(mean)
		labels.append('P={}'.format(P))
		return i 

	one_iteration(P,OX,labels,N)
	#P = round(P*.5,2)
	i = one_iteration(.7,OX,labels,N)
	#P = round(P*.8,2)
	i = one_iteration(.5,OX,labels,N)

	return OX,labels,i

def exp_full():
	OX = [] 
	labels = [] 
	P = .9

	N = 10

	def one_iteration(P,OX,labels,N,a,b):
		ISA = [] 
		i = 0 
		for j in range(N):
			O,i, IS = S.simulate_basic(a,b)
			ISA.append(np.array(IS))

		mean = np.mean(ISA, axis=0).tolist()
		OX.append(mean)
		labels.append('β={} ω={}'.format(a,b))
		return i 

	one_iteration(P,OX,labels,N,.9,.1)
	#P = round(P*.5,2)
	i = one_iteration(P,OX,labels,N,.5,.5)
	#P = round(P*.8,2)
	i = one_iteration(P,OX,labels,N,.1,.9)

	return OX,labels,i


OX,labels,i = exp_stepwise()

S.plot_multi_results(OX,i,labels)