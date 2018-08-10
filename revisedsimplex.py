### Revised Simplex Method! 
import numpy as np
import scipy.linalg as la

import numpy as np

class RevisedSimplex:
    
    def __init__(self, cn, An, b):
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
        cls.inves = cls.TriangularFactorization(cls)
        
        ## Optimality Flag
        cls.optimal = False
    
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
        
        #print(self.inves)
        
        def BTRAN(x):
            y = np.matrix(x)
            # print(self.inves)
            ## used to solve y 
            for i in self.inves[::-1]:
                position, col = i
                invE = np.matrix(np.eye(len(An)))
                invE[:,position] = np.matrix(col)
                #print(invE)
                y = np.dot(y , invE)
            return y

        y = BTRAN(self.cb)
        #print("y is:",y)
        # y = np.dot(self.cb, inv(self.Bnew))
        
        enter = self.cn - np.dot(y, self.An)
        entering = int(np.argmax(enter, axis = 1)) 
        entering_variable = self.nonbasic[entering]
        optimal = np.all(enter <= 0)

        print(enter, RevisedSimplex.cn)

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
                #print("the d vector is :",d )
                #print(d)
            return d
        
        d = FTRAN(self.Andic[entering_variable])
        #print("the d vector is :",d )
        #d = (np.dot( inv(self.Bnew), self.Andic[entering_variable] ))
        print("the d vector is :",d )
        calculate_t = np.divide(self.b, d) 

        exiting = np.argmin(calculate_t)
        #exiting_variable = basic[exiting]
        # this calculates the new b vector, need to func later to check if problem is unbounded!
        t = np.min(calculate_t) #np.min(self.b.T / d)
        bnew = self.b - t*d 

        #print("the shape of b is", b)
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
        RevisedSimplex.optimal = optimal
        
        # First swap entering and exiting in lists
        
        print(entering, exiting)
        
        self.nonbasic[entering], self.basic[exiting] = self.basic[exiting], self.nonbasic[entering]
        RevisedSimplex.nonbasic = self.nonbasic
        RevisedSimplex.basic = self.basic
        # Compute new cn and cb:
        RevisedSimplex.cn, RevisedSimplex.cb = swap(self.cn,self.cb,entering,exiting)
        # compute new An
        RevisedSimplex.An[:,entering] = self.B[:,exiting] # here B cannot change between iterations (B_new and B_old can)
        # compute new b vector!
        RevisedSimplex.b = bnew #.T

        #print("shape of d is,",d.shape)
        
        def Etainverse(matrix,position,a):
            E = matrix
            E[:,position] = a.reshape(-1,1)
            test = E[position,position]
            E[:,position] = -E[:,position]/test
            E[position,position] = 1/test
            Epack = tuple((position, E[:,position]))
            ## need to append tuple with position and column
            return Epack
        
        #print("in update d vector is", d)
        Epack = Etainverse(np.matrix(np.eye(len(An))), exiting, d)
        #print("in Update Epack is:", Epack)
        self.inves.append(Epack)
        RevisedSimplex.inves = self.inves

        print(RevisedSimplex.optimal)

    def Pivot(self):
        # counts the number of pivtos 
        count = 0
        while RevisedSimplex.optimal != True:
            self.Update()
            count += 1    
        print("The current dictionary is optimal! The number of iterationrs to solve was:", count) 
        
       #print("current cb is",RevisedSimplex.cn)
    
    def solve():
        # cetral function: calls and organizes all other functions
        # returns optimal coefficients, optimal Z, number of iterations to solve, maybe time?
        pass

    
    ### --- Helper functions --- ###
    
    def cycling():
        # intended to catch cycling behaviour in problems
        # maybe keep track of variables and check if it cycling through the same variables
        # do this later, cycling almost never occurs
        pass
    
    def maximum(self):
        if RevisedSimplex.optimal == True:
            ## calculation of z 
           # z =
            pass
        pass
    # @staticmethod
    def TriangularFactorization(cls):
        (P, L, U) = la.lu(cls.B) #B matrix does not have to be I :-)
        
        def Etafact(matrix):
            d = np.matrix(matrix)
            Ilist = []
            for i in range(len(matrix)):
                col = d[:,i]
                I = np.matrix(np.eye(len(matrix)))
                I[:,i] = col
                Ipack = tuple((i, I[:,i]))
                Ilist.append(Ipack)
            return Ilist
        l = Etafact(L)
        p = Etafact(P)
        u = Etafact(U)
        
        def createETA(p,l,u):
            c = [j for i in zip(p,l) for j in i]
            Um = u[::-1]
            new = c + Um
            return new
        initiate = createETA(p,l,u)
        return initiate