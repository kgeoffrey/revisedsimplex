### Revised Simplex Method
### © Geoffrey Kasenbacher

import numpy as np
import scipy.linalg as la
## Division by zero will yield inf - only care about the minimum value in vector 
np.seterr(divide='ignore', invalid='ignore')

class RevisedSimplex:
    
    def __init__(self, cn, An, b):
        ## as of now Basis will always be I, have to write function 
        ## that will catch arbitrary A and split it into An and B
        ## --> LU decomposition allows for this
        self.Initializer(cn, An, b)
        
    @classmethod
    def Initializer(cls, cn, An, b): 
        cls.An = An
        cls.B = (np.eye(len(An)))
        cls.b = b
        cls.cn = cn
        cls.cb = np.zeros(len(An))
        ## variable list constructor
        def lis(x, length):
            L = []
            for i in range(x):
                L.append("x" + str(i + length))
            return L
        cls.nonbasic = lis(len(An.T), 1)
        cls.basic = lis(len(An), 1 + len(cls.nonbasic))
        cls.Andic = dict(zip(cls.nonbasic, An.T[:,]))

        ## Decompose B!
        cls.inves = cls.TriangularFactorization(cls.B)
        
        ## Optimality Flag
        cls.optimal = False
        
        ## disabled trinagular fact. if initial basis is I
        '''
        one = np.matrix([1])
        zeros = np.zeros(len(An)-1).reshape(len(An)-1,1)
        cls.inves = [(0,np.vstack((one, zeros)))]
        '''
    def Enter(self):
         
        self.basic = RevisedSimplex.basic
        self.nonbasic = RevisedSimplex.nonbasic
        self.b = RevisedSimplex.b
        self.An = RevisedSimplex.An
        self.cb = RevisedSimplex.cb
        self.cn = RevisedSimplex.cn
        self.inves = RevisedSimplex.inves
        self.Andic = dict(zip(self.nonbasic, self.An.T[:,]))
        self.B = (np.eye(len(An)))
        
        def BTRAN(x):
            y = np.matrix(x) 
            for i in self.inves[::-1]:
                position, col = i
                invE = np.matrix(np.eye(len(An)))
                invE[:,position] = np.matrix(col)
                y = np.dot(y , invE)
            return y

        y = BTRAN(self.cb)
        
        enter = self.cn - np.dot(y, self.An)
        entering = int(np.argmax(enter, axis = 1)) 
        entering_variable = self.nonbasic[entering]
        optimal = np.all(enter <= 0)

        return entering, entering_variable, optimal 
    
    def Exit(self):
        entering, entering_variable, optimal = self.Enter(self)
        
        def FTRAN(x):
            d = x.reshape(len(x),-1)
            for i in self.inves:
                position, col = i
                invE = np.matrix(np.eye(len(An)))
                invE[:,position] = np.matrix(col)
                d = np.dot(invE , d)
            return d
        
        d = FTRAN(self.Andic[entering_variable])

        calculate_t = np.divide(self.b, d) 

        exiting = np.argmin(calculate_t)
        
        t = np.min(calculate_t)
        bnew = self.b - t*d 
        
        ## if checks if the problem is unbounded
        if np.any(bnew < 0) == True:
            raise Exception('This problem is unbounded!')
            
        bnew[bnew == 0] = t
        return entering, exiting, bnew, d, optimal
    
    @classmethod 
    def Update(self):
        ### --- Helper functions ---
        def changecollumn(matrix,position,d):
            E = matrix.copy()
            E[:,position] = np.matrix(d)
            return E
        def swap(A,B,i,j):
            TEMP_B = B.T[j].copy()
            B.T[j] = A.T[i]
            A.T[i] = TEMP_B
            return A,B
        
        entering, exiting, bnew, d, optimal = self.Exit(self)
        ## This terminates Algorithm, when the optimal dictionary is found
        RevisedSimplex.optimal = optimal
        
        # This file refactors the ETA file if it gets too long!
        if len(RevisedSimplex.inves) > len(RevisedSimplex.inves) + 20:
            self.refactorize()
        
        if RevisedSimplex.optimal == True:
            return
        
        print(entering, exiting)
        
        self.nonbasic[entering], self.basic[exiting] = self.basic[exiting], self.nonbasic[entering]
        RevisedSimplex.nonbasic = self.nonbasic
        RevisedSimplex.basic = self.basic
        RevisedSimplex.cn, RevisedSimplex.cb = swap(self.cn,self.cb,entering,exiting)
        RevisedSimplex.An[:,entering] = self.B[:,exiting] 
        RevisedSimplex.b = bnew 
        
        def Etainverse(matrix,position,a):
            E = matrix
            E[:,position] = a.reshape(-1,1)
            test = E[position,position]
            E[:,position] = -E[:,position]/test
            E[position,position] = 1/test
            Epack = tuple((position, E[:,position]))
            return Epack
        
        Epack = Etainverse(np.matrix(np.eye(len(An))), exiting, d)
        self.inves.append(Epack)
        RevisedSimplex.inves = self.inves


    def Pivot(self):
        # counts the number of pivtos 
        count = 0
        while RevisedSimplex.optimal != True:
            self.Update()
            count += 1
        return count
    
    def solve(self):
        count = self.Pivot()
        if RevisedSimplex.optimal == True:
            z = self.maximum()
            print("The current dictionary is optimal! The number of iterations to solve was:", count) 
            print("maximum value is:", z)

    
    ### --- Helper functions --- ###
    
    def tabulate():
        # will write function that puts results in table
        # also want to show slack per variable
        # small table
        pass
    
    def cycling():
        # cycling is rare but is guaranteed to never occur under "bland's rule"
        # -> chose smallest subscript variable if there is a choice 
        # do this later, cycling almost never occurs
        pass
    
    def maximum(self):
        if RevisedSimplex.optimal == True: 
            z = np.dot(RevisedSimplex.cb,  RevisedSimplex.b)
            return z
        
    @staticmethod
    def TriangularFactorization(X):
        (P, L, U) = la.lu(X) #B matrix does not have to be I :-)
        
        def InvEtafact(matrix):
            d = np.matrix(matrix)
            Ilist = []
            for i in range(len(matrix)):
                col = d[:,i]
                I = np.matrix(np.eye(len(matrix)))
                I[:,i] = col
                test = I[i,i]
                I[:,i] = -I[:,i]/test
                I[i,i] = 1/test
                Ipack = tuple((i, I[:,i]))
                Ilist.append(Ipack)
            return Ilist
        l = InvEtafact(L)
        p = InvEtafact(P)
        u = InvEtafact(U)
        
        def createETA(p,l,u):
            c = [j for i in zip(p,l) for j in i]
            Um = u[::-1]
            new = c + Um
            return new
        initiate = createETA(p,l,u)
        return initiate
    
    @classmethod
    def refactorize(self):
        ## This function will refactorize B the current B when the ETA file
        ## or "inves" gets too long, maybe after 20 iterations? 
        def reconstruct():
            y = np.eye(len(RevisedSimplex.An))
            for i in RevisedSimplex.inves:
                pos, col = i 
                I = np.matrix(np.eye(len(RevisedSimplex.An)))
                # I[:,i] = col
                test = I[pos,pos]
                I[:,pos] = -I[:,pos]/test
                I[pos,pos] = 1/test
                y = np.dot(y , I)
            return y
        newB = reconstruct()
        RevisedSimplex.inves = RevisedSimplex.TriangularFactorization(newB)
