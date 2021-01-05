import random,os,subprocess,time,csv
from math import ceil, floor 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sklearn.preprocessing import normalize
import statistics as stats 

class BARON:
	def __init__(self):
		pass 
	def options(self,findex='0',path=''):
		cplex_lib_path = os.environ['CPLEX_LIB']
		fpath = path+str(findex) 
		opt = """OPTIONS {{
				results: 1;
				PrLevel: 0;
				resname: \"{}.res\";
				TimName: \"{}.tim\";
				OptName: \"{}.opt\";
				MaxIter: 100; 
				CplexLibName: \"{}\";}}""".format(fpath,fpath,fpath,cplex_lib_path)
		return opt 

	def binary_vars(self,indexes,prefix='x'):
		if len(indexes) < 1: 
			return ''
		v = "BINARY_VARIABLES "
		for i in indexes: 
			v+= prefix+str(i)+","
		return v[:-1]+";"

	def positive_vars(self,indexes,prefix='x'):
		if len(indexes) < 1: 
			return ''
		v = "POSITIVE_VARIABLES "
		for i in indexes: 
			if len(prefix) == 0: 
				v+= i+","
			else: 
				v+= prefix+str(i)+","
		return v[:-1]+";"

	def declare_equations(self,indexes,prefix='e'):
		v = "EQUATIONS "
		for i in indexes: 
			if len(prefix) == 0: 
				v+= i+","
			else: 
				v+= prefix+str(i)+","
		return v[:-1]+";"

	def obj_fn(self, fn):
		return 'OBJ: '+fn+';'

	def out_file(self,program,findex,fpath):
		opts  = self.options(findex,'baron/')
		f = open(fpath,'w')
		print(opts,file=f)
		print('',file=f)
		print(program['binary_vars'],file=f)
		print('',file=f)
		print(program['positive_vars'],file=f)
		print('',file=f)
		print(program['equations'],file=f)
		print('',file=f)
		print(program['obj'],file=f)

		f.close()

	def run_exec(self,path):
		bpath = os.environ.get('BARON_PATH')
		os.execl(bpath, "baron", path)

	def run_subp(self,path):
		bpath = os.environ.get('BARON_PATH')
		cmd = [bpath, path]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE) 
		p.wait()

	def read_sol(self,index=0):
		path = '{}/baron/{}.res'.format(os.environ['IOTDIV_PATH'],index)
		f = open(path,'r')
		for line in f: 
			if line.find('The best solution found is') > -1: 
				next(f)
				next(f)
				break 
		V = [] 
		for line in f: 
			if len(line) < 2: 
				break 
			line = line.split('			')
			sol = eval(line[2])
			v = line[0].lstrip()
			V.append((v,sol))

		objective_value = 0 
		for line in f: 
			if line.find('The above solution') > -1: 
				s = line.split(':')
				s = s[1].lstrip().rstrip()
				objective_value = eval(s)
		V.append(('OBJ',objective_value))
		return V


class BARONiotdiv(BARON):
	def __init__(self):
		pass 
	def row_constraints(self,indexes,m,c):
		limit = 1  
		equations = ''
		for i in range(m):
			equations+='e{}: '.format(indexes[i])
			for j in range(c):
				equations+= 'x{}_{}'.format(i,j)+'+'
			equations=equations[:-1]
			equations+= '=={};\n'.format(limit) 
		return equations 

	def col_constraints(self,indexes,m,c,diversity):
		limit = diversity 
		equations = ''
		for i in range(c):
			equations+='e{}: '.format(indexes[i])
			for j in range(m):
				equations+= 'x{}_{}'.format(j,i)+'+'
			equations=equations[:-1]
			equations+= '<={};\n'.format(limit) 
		return equations 

	def obj_fun_constraints(self,indexes,R,m,c,index,var='ro'):
		#For each column add a coefficient. Repeat for each row. 
		equation = 'e{}: '.format(index) 
		for i in range(m):
			for j in range(c): 
				equation+= '{}*x{}_{}+'.format(R[j],i,j)
		equation=equation[:-1]
		equation+= '-{}==0;\n'.format(var) 
		return equation

	def process_sol(self,X): 
		V = self.read_sol()
		obj = 0 
		for v in V:
			s = v[1]
			v = v[0] 
			if v[0] == 'x':
				S = v.split('_')
				i = eval(S[0][1:])
				j = eval(S[1])
				X[i][j] = s 
			if v.find('OBJ')>-1: 
				obj = s 
		return X, obj 

	def gen_program(self,beta,omega,m,c,D,R,I):

		#Objective function:
		fn = "minimize {}*ro+{}*sigma".format(beta,omega)
		obj = self.obj_fn(fn)

		#Variables: 
		positive_vars = self.positive_vars(['ro','sigma'],prefix='')

		#Binary variables: 
		indexes = [] #an entry for each element of matrix X 
		for i in range(m):
			for j in range(c):
				indexes.append(str(i)+'_'+str(j))
		binary_vars   = self.binary_vars(indexes)
		
		#Equations: 
		no_eq = m 
		indexes = [i for i in range(no_eq)] #an entry for each row in X 
		#an entry for each column in X 
		cindexes = [i for i in range(no_eq,c+no_eq)]
		no_eq +=c 

		#declare row/col equations
		equations = self.declare_equations(indexes+cindexes+[no_eq,no_eq+1])

		#add row constraints
		equations+= '\n'+self.row_constraints(indexes,m,c)

		#add col constraints (controls diversity for each configuration)
		equations+= self.col_constraints(cindexes,m,c,D)

		equations+= self.obj_fun_constraints(indexes,R,m,c,no_eq,var='ro')
		equations+= self.obj_fun_constraints(indexes,I,m,c,no_eq+1,var='sigma')

		program = dict() 
		program['equations'] = equations
		program['binary_vars'] = binary_vars
		program['positive_vars'] = positive_vars
		program['obj'] = obj 

		return program 

class Simulator:
	def __init__(self):
		pass 

	def execute_baron(self,program,X):
		BR = BARONiotdiv()
		path = os.environ['IOTDIV_PATH']+'/baron/program.bar'
		start = time.time()
		BR.out_file(program,'0',path)
		BR.run_subp(path)
		end = time.time()

		#print('Solution produced in {} seconds.'.format(round(end-start,4))) 

		return BR.process_sol(X)

	def normalize(self,V): 
		V = np.array(V).reshape(1,-1)
		V = normalize(V)
		V = V.tolist()[0]
		return V  


	def update_R(self,X,R,c): 
		'''
		Examine each row, add one to the counter for the
		configurations that were used in the selected row. 
		'''
		for x in X: 
			R = [R[i] if x[i]==0 else R[i]+1 for i in range(c)]
		return R #self.normalize(R)

	def update_I(self,X,S,I,c): 
		'''
		Examine each attacker skill, add one to the counter for the
		configurations that were used AND compromised. 
		'''
		for x in X: 
			I = [I[i]+1 if S[i]==1 and x[i]==1 else I[i] for i in range(c)]
		return I #self.normalize(I)

	def count_compromised(self,X,S,I,c):
		'''
		Examine each attacker skill, add one to the compromised counter, 
		if the skill matches the chosen configuration. 
		'''
		compromised = 0 

		for x in X: 
			c = sum([1 if S[i]==1 and x[i]==1 else 0 for i in range(c)]) 
			if c > 0: 
				compromised+=1 
		return compromised

	def update_attacker(self, sample_space, n, c):
		chosen_indexes = random.sample(sample_space, n)
		S = [1 if i in chosen_indexes else 0 for i in range(c)]
		return S 

	def plot_results(self,R, I, S, O, iterations):
		fig, ax = plt.subplots()

		x = [i for i in range(iterations)]
		y = [item for item in O]
		
		line1, = ax.plot(x, y, label='Objective function value')
		#line1.set_dashes([2, 2, 10, 2]) 
		
		#y = [item for item in I]
		#line2, = ax.plot(x, y, label='Insecurity value')

		ax.legend()
		plt.savefig("fig.PDF",dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)

	def plot_multi_results(self, OX, iterations, labels):
		plt.rcParams['font.size'] = 16
		plt.rcParams['font.weight'] = 'bold'

		fig, ax = plt.subplots()
		k = 1
		markers = ['','*','o','+','-','x']
		for O,LB in zip(OX,labels): 
			x = [i for i in range(iterations)]
			y = [item for item in O]
			#f2 = interp1d(x, y, kind='cubic')
			#M = max(x)
			#xnew = np.linspace(0, M, endpoint=True)

			line1, = ax.plot(x, y, label=LB)
			#plt.plot(x, y, 'o', xnew, f2(xnew), '--')
			k+=1 
		
		

		ax.legend()
		ax.set_xlabel('Iterations', fontsize=18)
		ax.set_ylabel('Insecurity cost', fontsize=18)
		ax.grid(color='.95',linestyle='-', linewidth=.5)
		#ax.set_title('')
		fig.tight_layout()
		plt.savefig("fig.PDF",dpi=None, facecolor='w', edgecolor='w',
		orientation='portrait', format=None,
		transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)

	def record_attackers(self, S):
		f = open("attackers.csv","a+")
		writer = csv.writer(f)
		writer.writerow(S) 
		f.close() 

	def read_attackers(self): 
		try: 
			f = open("attackers.csv", "r")
		except Exception as e: 
			return None 
		reader = csv.reader(f)
		AS = [list(map(int,x)) for x in reader]
		return AS 

	def load_configs(self):
		#Number of configurations
		c = 32
		#Number of devices 
		d = 8
		#Machines per devices 
		k = 1 
		#Machines
		m = d*k 
		#Number of attackers 
		nA = 4

		#Configuration matrix
		X = [[0 for i in range(c)] for i in range(m)] 
		#Diversity parameter 
		D = c // m 
		if D < 1: 
			D = 1 
		#Configuration repetition counter
		R = [0 for i in range(c)]
		#Configuration insecurity counter 
		I = [0 for i in range(c)]
		#Attacker capabilities 
		p = .2 
		#Choose the space of samples for configurations.
		#If each config appears once, then they all have the same
		#probability of being vulnerable. We'll use a slightly
		#different weight for some of them. 
		sample_space = [i for i in range(c)]
		#Repeat the first four configs twice each.
		#That means, these first four are more vulnerable. 
		for i in range(2):
			sample_space += [i for i in range(0,4)]

		return c,d,k,m,nA,X,D,R,I,p,sample_space

	def normalized_mean(self, V): 
		X = np.array([x for x in V])
		X = normalize(X.reshape(1, -1))
		m = np.mean(X)
		print(m)
		raise Exception('nothing')
		return m 

	def normalized_max(self, V): 
		X = np.array([x for x in V])
		X = normalize(X.reshape(1, -1)).tolist()
		X = [sum(x) for x in X]
		#print(X[0])
		m = stats.mean(X)
		#print(m)
		#raise Exception('nothing')
		return m 


	def simulate_basic(self,beta,omega):
		c,d,k,m,nA,X,D,R,I,p,sample_space = self.load_configs()

		print('Simulating defense with beta {} and omega {}'.format(beta,omega))

		n = ceil(c*p)

		BR = BARONiotdiv()

		#record of objective function values. 
		O = []
		iterations = 150

		#record insecurity scores
		IS = [] 

		#Reset repetition rates every xr ticks. 
		xr = ceil(iterations/ 100000) 

		#If we have logged attacker, we'll use them. 
		AS = self.read_attackers()
		AS = [] 
		if AS: 
			l = len(AS)
		else: 
			l = 0  

		for i in range(iterations):
			#print("({})".format(i))
			#Start generating BARON program
			program = BR.gen_program(beta,omega,m,c,D,self.normalize(R),self.normalize(I))
			X,obj = self.execute_baron(program,X)
			O.append(obj) 

			#print('Obj: {}\nConfig: {}'.format(obj,X))
			
			#Repeat for the number of attackers.
			if l > 0:
				S = AS[l-1]
				l-=1  
			else: 
				S = self.update_attacker(sample_space, n, c)
				self.record_attackers(S)

			#print('Attacker skills:',S)
			R = self.update_R(X,R,c)
			I = self.update_I(X,S,I,c)

			# print(R)
			# print(I)
			# print(S)
			# print('X:')
			# for x in X:
			# 	print(x)
			# print(obj)

			IS.append(self.normalized_mean(I))
			
			# print('Repetition rates:',R)
			# print('Insecurity rates:', I)
			# print('Attacker, repetition, and insecurity rates updated.')

			xr-=1 
			if (xr == 0):
				xr = ceil(iterations/ 10) 
				R = [0 for i in range(c)]

		return O,i+1, IS

	def expand_attacker_skills(self, S): 
		T = [] 
		c = len(S)
		for i in range(c): 
			if S[i] == 1: 
				T.append([0 if j!= i else 1 for j in range(c)])
		return T 

	'''
	This type of the attacker is persistent and only reveals their skills 
	gradually. That is, the attacker can compromise two types of the OS, 
	they will not do so in one attack attempt. 
	The attacker starts with one configuration and moves to another 
	one at a time. This is expected to confuse the defense. 
	'''
	def simulate_stepwise(self,beta,omega): 
		c,d,k,m,nA,X,D,R,I,p,sample_space = self.load_configs()

		print('Simulating defense with beta {} and omega {}'.format(beta,omega))

		n = ceil(c*p)

		BR = BARONiotdiv()

		#record of objective function values. 
		O = []
		iterations = 150

		#record insecurity scores
		IS = [] 

		#Reset repetition rates every xr ticks. 
		xr = ceil(iterations/ 10) 

		#If we have logged attacker, we'll use them. 
		AS = self.read_attackers()
		if AS: 
			l = len(AS)-1
		else: 
			l = -1  

		ASX = [] 
		for i in range(l):
			T = self.expand_attacker_skills(AS[i])
			#T*=4 
			ASX += T 

		for i in range(iterations):
			#print("({})".format(i))
			#Start generating BARON program
			program = BR.gen_program(beta,omega,m,c,D,self.normalize(R),self.normalize(I))
			X,obj = self.execute_baron(program,X)
			O.append(obj) 

			if len(ASX) == 0: 
				break 

			S = ASX.pop() 

			R = self.update_R(X,R,c)
			I = self.update_I(X,S,I,c)

			# print(R)
			# print(I)
			# print(S)
			# print('X:')
			# for x in X:
			# 	print(x)
			# print(obj)

			IS.append(self.normalized_mean(I))

		return O, i+1, IS


	'''
	This simulator uses a defense that is not accurate in 
	detecting attacks. That is, for some probability P, the
	defense can detect if a specific configuration was compromised. 
	'''

	def simulate_inaccurate(self,beta,omega,P=.7):
		c,d,k,m,nA,X,D,R,I,p,sample_space = self.load_configs()

		recorded_I = I[:]

		print('Simulating defense with probability {}'.format(P))
		n = ceil(c*p)

		BR = BARONiotdiv()

		#record of objective function values. 
		O = []
		iterations = 150

		#record insecurity scores
		IS = [] 

		#Reset repetition rates every xr ticks. 
		xr = ceil(iterations/ 100000) 


		#If we have logged attacker, we'll use them. 
		AS = self.read_attackers()
		AS = [] 
		if AS: 
			print('Loaded logged attacker.')
			l = len(AS)
		else: 
			l = 0  

		for i in range(iterations):
			#print("({})".format(i))
			#Start generating BARON program
			program = BR.gen_program(beta,omega,m,c,D,self.normalize(R),self.normalize(I))
			X,obj = self.execute_baron(program,X)
			O.append(obj) 

			#print('Obj: {}\nConfig: {}'.format(obj,X))
			
			#Repeat for the number of attackers.
			if l > 0:
				S = AS[l-1]
				l-=1  
			else: 
				S = self.update_attacker(sample_space, n, c)
				self.record_attackers(S)

			#print('Attacker skills:',S)
			#Update repetition score
			R = self.update_R(X,R,c)

			'''Update insecurity score given probability P.
			We will keep two copies. One copy represents the true
			insecurity that will be used to examine the results 
			and one is used to compute the next configuration. '''

			choices = random.choices([1,0], weights = [P,1-P], k=1)
			choice = choices[0]
			if (choice == 1):
				I = self.update_I(X,S,I,c)
			recorded_I = self.update_I(X,S,recorded_I,c)

			IS.append(self.normalized_mean(recorded_I))
		

			#Reset repetition rate
			# xr-=1 
			# if (xr == 0):
			# 	xr = ceil(iterations/ 10) 
			# 	R = [0 for i in range(c)]

		return O,i+1, IS
		
		
		