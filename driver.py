from simulator import Simulator,SimulationPlots
import os, statistics as stats 
import numpy as np
import sys,csv

iterations = 150 
N = 10
multiple = False 
c = 2**5

parameters = [(.9,.1),(.5,.5),(.1,.9)]

def exp_stepwise():
	OX = [] 
	labels = [] 
	P = .9
	S = Simulator(c=c)

	def one_iteration(P,OX,N,a,b):
		ISA = [] 
		i = 0 
		for j in range(N):
			O,i, IS = S.simulate_stepwise(a,b,iterations = iterations)
			ISA.append(np.array(IS))

		mean = np.mean(ISA, axis=0).tolist()
		OX.append(mean)
		return i 

	for p in parameters:
		i = one_iteration(P,OX,N,p[0],p[1])
		labels.append('β={},Ω={}'.format(p[0],p[1]))

	return OX,labels,i

def exp_inaccurate():
	OX = [] 
	labels = [] 
	P = .9

	def one_iteration(P,OX,N):
		ISA = [] 
		i = 0 
		for j in range(N):
			O,i, IS = S.simulate_inaccurate(.5,.5,P,iterations=iterations)
			ISA.append(np.array(IS))

		mean = np.mean(ISA, axis=0).tolist()
		OX.append(mean)
		return i 

	parameters = (2**3, 2**4, 2**5, 2**6, 2**7)
	for p in parameters:
		S = Simulator(c=p)
		i = one_iteration(P,OX,N)
		
		labels.append('C={}'.format(p))

	return OX,labels,i


def exp_random():
	OX = [] 
	labels = [] 
	P = .9

	def one_iteration(P,OX,N):
		ISA = [] 
		i = 0 
		for j in range(N):
			O,i, IS = S.simulate_randomized(.9,.1,P,iterations=iterations,inaccurate=True)
			ISA.append(np.array(IS))

		mean = np.mean(ISA, axis=0).tolist()
		OX.append(mean)
		return i 

	S = Simulator(c=c)
	i = one_iteration(P,OX,N)
	labels.append('Random'.format(c))

	return OX,labels,i

def exp_full():
	OX = [] 
	labels = [] 
	P = .9

	S = Simulator(c=c)
	def one_iteration(P,OX,N,a,b,iterations=iterations):
		ISA = [] 
		i = 0 
		for j in range(N):
			O,i, IS = S.simulate_basic(a,b)
			ISA.append(np.array(IS))

		mean = np.mean(ISA, axis=0).tolist()
		OX.append(mean)
		
		return i 

	for p in parameters:
		i = one_iteration(P,OX,N,p[0],p[1])
		labels.append('Complete (β={},Ω={}'.format(p[0],p[1]))

	return OX,labels,i

def process_res(OX,i,labels):
	S = SimulationPlots() 
	f = open('results.data','w')
	writer = csv.writer(f)
	for O,LB in zip(OX,labels): 
		writer.writerow(O+[LB])
	writer.writerow([i])
	f.close()
	S.plot_multi_results(OX,i,labels)

def load_res():
	f = open('results.data','r')
	reader = csv.reader(f)
	OX=[]
	labels = [] 
	i = 0 
	for row in reader: 
		if len(row) == 1: 
			i = eval(row[0])
			break 
		data = row[:-1]
		data = list(map(eval,data))
		# N = len(data)
		# avg = sum(data)/N
		# data = [round(x/avg,2) for x in data]
		OX.append(data)
		labels.append(row[-1])
	f.close() 
	return OX,labels,i

if sys.argv[1] == 'rand': 
	OX,labels,i = exp_random()
	process_res(OX,i,labels)

elif sys.argv[1] == 'inacc': 
	OX,labels,i = exp_inaccurate()
	process_res(OX,i,labels)

elif sys.argv[1] == 'basic': 
	OX,labels,i = exp_full()
	process_res(OX,i,labels)

elif sys.argv[1] == 'step': 
	OX,labels,i = exp_stepwise()
	process_res(OX,i,labels)

elif sys.argv[1] == 'mix': 
	OX,labels,i = exp_inaccurate()
	OX2,labels2,i2 = exp_random()
	OX += OX2 
	labels += labels2
	process_res(OX,i,labels)

elif sys.argv[1] == 'full-rand': 
	OX,labels,i = exp_full()
	OX2,labels2,i2 = exp_random()
	OX += OX2 
	labels += labels2
	process_res(OX,i,labels)

elif sys.argv[1] == 'all': 
	OX = [] 
	labels = [] 
	OX,labels,i = exp_full()
	OX2,labels2,i = exp_random()
	OX3,labels3,i = exp_inaccurate()
	OX4,labels4,i = exp_stepwise()
	OX += OX2 
	OX += OX3 
	OX += OX4
	labels += labels2 
	labels += labels3
	labels += labels4
	process_res(OX,i,labels)

elif sys.argv[1] == 'plot': 
	OX,labels,i = load_res()
	S = SimulationPlots() 
	S.plot_multi_results(OX,i,labels)



