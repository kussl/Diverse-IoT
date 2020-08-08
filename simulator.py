import random,os,subprocess
from math import ceil, floor 

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
		cmd = [bpath, 'baron', path]
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
				equations+= 'x{}{}'.format(i,j)+'+'
			equations=equations[:-1]
			equations+= '=={};\n'.format(limit) 
		return equations 

	def col_constraints(self,indexes,m,c,diversity):
		limit = diversity 
		equations = ''
		for i in range(c):
			equations+='e{}: '.format(indexes[i])
			for j in range(m):
				equations+= 'x{}{}'.format(j,i)+'+'
			equations=equations[:-1]
			equations+= '<={};\n'.format(limit) 
		return equations 

	def obj_fun_constraints(self,indexes,R,m,c,index,var='ro'):
		#For each column add a coefficient. Repeat for each row. 
		equation = 'e{}: '.format(index) 
		for i in range(m):
			for j in range(c): 
				equation+= '{}*x{}{}+'.format(R[j],i,j)
		equation=equation[:-1]
		equation+= '-{}==0;\n'.format(var) 
		return equation

	def process_sol(self,X): 
		V = self.read_sol()
		for v in V:
			s = v[1]
			v = v[0] 
			if v[0] == 'x':
				i = eval(v[1])
				j = eval(v[2])  
				X[i][j] = s 
		return X 

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
				indexes.append(str(i)+str(j))
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

	def next_config(self,program,X):
		BR = BARONiotdiv()
		path = os.environ['IOTDIV_PATH']+'/baron/program.bar'
		BR.out_file(program,'0',path)
		BR.run_subp(path)
		print('OK') 
		return BR.process_sol(X)


	def simulate_basic(self,beta,omega):
		#Number of configurations
		c = 8
		#Number of devices 
		d = 4 
		#Machines per devices 
		k = 1 
		#Machines
		m = d*k 
		#Configuration matrix
		X = [[0 for i in range(c)] for i in range(m)] 
		#Diversity parameter 
		D = m // c 
		if D < 1: 
			D = 1 
		#Configuration repetition counter
		R = [0 for i in range(c)]
		#Configuration insecurity counter 
		I = [0 for i in range(c)]
		#Attacker capabilities 
		p = .2 
		chosen_indexes = random.sample([i for i in range(c)], ceil(c*p))
		S = [1 if i in chosen_indexes else 0 for i in range(c)]

		#Start generating BARON program
		BR = BARONiotdiv()

		program = BR.gen_program(beta,omega,m,c,D,R,I)
		X = self.next_config(program,X)

		
		
		