# -*- coding: utf-8 -*-

from   math import exp, factorial
import numpy as np
from   numpy.linalg import matrix_power
from   scipy.stats  import poisson
from   scipy.linalg import norm, null_space, solve, solve_sylvester, expm, inv
import matplotlib.pyplot as plt
import sys, warnings
from   tqdm import tqdm
'''
̄W combining macron, treated as a special character by python; may not be correctly rendered by pyplot.
W̅ combining overline, treated as a letter by python; may not be correctly rendered by pyplot.
'''

def QuadraticMatrixEquationMinimalPositiveSolution(A=1, B=1, C=1, *, dim:int=None, method:str='Bernoulli', epsilon:float=1.e-10, maxiter:int=None):
    """QuadraticMatrixEquationMinimalPositiveSolution

    Calculates the minimal elementwise-positive solvent the quadratic system `A X^2 + B X + C == O`, where `A, B, C, X, O` are all square matrices.
    The existence and uniqueness of such solvent is demonstrated in https://dx.doi.org/10.1016/j.amc.2011.08.070

    Args:
        A (matrix, optional): square matrix or real number. Defaults to 1.
        B (matrix, optional): square matrix or real number. Defaults to 1.
        C (matrix, optional): square matrix or real number. Defaults to 1.
        the provided A, B, C must be of the same order or a real number; if some of them are provided as real numbers, the dim of the system must be given explicitly.
        dim (int, optional): the size of the square matrices. Defaults to None.
        method (str, optional): 'Bernoulli' or 'Newton'. Defaults to 'Bernoulli'.
        epsilon (float, optional): Tolerance of the norm of the residue. Defaults to 1.e-10.
        maxiter (int, optional): Maximum rounds of iteration. Defaults to None.

    Raises:
        ValueError: Raised when the coefficients are not coherent square matrices.
        NotImplementedError: When a method other than 'Bernoulli' or 'Newton' is specified.

    Returns:
        np.ndarray: the minimal positive solution to the quadratic system
    """
    if isinstance(A, (int, float)):
        if type(dim) is int: A = A*np.eye(dim)
        else: raise TypeError("The quadratic coefficient is provided as a real number, but the dimensionality is not correctly specified.")
    elif type(dim) is int and A.shape != (dim, dim): raise ValueError("The quadratic coefficient is not of specified shape.")
    if type(dim) is not int: dim = A.shape[0]
    if isinstance(A, (int, float)): A = A*np.eye(dim)
    elif type(dim) is int and A.shape != (dim, dim): raise ValueError("The quadratic coefficient is not of specified shape.")
    if type(dim) is not int: dim = A.shape[0]
    if isinstance(A, (int, float)): A = A*np.eye(dim)
    elif type(dim) is int and A.shape != (dim, dim): raise ValueError("The quadratic coefficient is not of specified shape.")
    if type(dim) is not int: dim = A.shape[0]
    try:
        B = solve(A, B)
        C = solve(A, C)
    except:
        warnings.warn(f"Warning: Exceptions encountered normalizing the linear system.")
    def BernoulliIteration(X)->np.ndarray: return solve(          X+B,       -C)
    def    NewtonIteration(X)->np.ndarray: return solve_sylvester(X+B, X, X@X-C)
    if   method=='Bernoulli': iter = BernoulliIteration
    elif method==   'Newton': iter =    NewtonIteration
    else : raise NotImplementedError(f"the specified iteration scheme {method} is not implemented.")
    X = np.zeros(A.shape)
    if type(epsilon) is float:
        if maxiter==None:
            while norm(iter(X)-X)>epsilon: X=iter(X)
        else:
            cnt = 0
            while norm(iter(X)-X)>epsilon and cnt<maxiter: X=iter(X); cnt=cnt+1
    else:
        cnt = 0
        while cnt<maxiter: X=iter(X); cnt=cnt+1
    return X
solve_qme = QuadraticMatrixEquationMinimalPositiveSolution

class RandomVariable:...
class PointProcess(RandomVariable):
    def __init__(self, rate:float=1):
        if rate<=0: raise ValueError(f"The rate should be a positive real, but {rate} is provided.")
        else: self.rate = rate
class MarkovianArrivalProcess(PointProcess):
    def __init__() -> None: raise NotImplementedError
class ArrivalProcess(PointProcess):
    def __init__(self, rate:float=1, *, dist=None):
        super(ArrivalProcess, self).__init__(rate)
        self.λ = self.lbd = self.arrival_rate = self.rate
        if    dist == None : self._dist = lambda x: (1/self.λ)*exp(-x/self.λ)
        elif callable(dist): self._dist = dist
        else: raise TypeError("The provided parameter is not a callable object.")
    def dist(self, x:float)->float: return self._dist(x)
    def moment(self, k:int  )->float: return NotImplementedError("The moment computation method for this process has not been implemented.")

    def as_MAP(self)->MarkovianArrivalProcess:
        """ Convert the arrival process into a Markovian one."""
        raise NotImplementedError("The MAP for this distribution is currently unknown.")
class ServiceProcess(PointProcess):
    def __init__(self, rate:float=1):
        super(ServiceProcess, self).__init__(rate)
        self.μ = self.mu = self.service_rate = self.rate
    def   dist(self, x:float)->float: return (1/self.μ)*exp(-x/self.μ)
    def moment(self, k:int  )->float: return self.μ**k
class QueuingProcess(PointProcess):
    def __init__(self, arrival_process:ArrivalProcess, service_process:ServiceProcess, *, service_policy:str="FIFO"):
        self.arrival_process = arrival_process
        self.service_process = service_process
        self.λ = self.lbd = self.arrival_rate = self.arrival_process.rate
        self.μ = self.mu     = self.service_rate = self.service_process.rate
        if self.service_rate<=0 :
            raise ValueError(f"the specified service rate is not a positive real number")
        else:
            self.μ = self.mu = self.service_rate = self.service_rate
            if self.service_rate<=self.arrival_rate:
                warnings.warn(f"Service rate {self.service_rate:6f} is slower than the arrival rate {self.arrival_rate:6f}. The system is not positive-recurrent.",stacklevel=3)
        self.ρ = self.rho = self.utilization_factor = self.arrival_rate/self.service_rate
        self.policy = self.service_policy = service_policy
    def is_positive_recurrent(self)->bool:
        return True if self.service_rate > self.arrival_process.rate else False
    def as_MAP(self):
        queue = QueuingProcess(self.arrival_process.as_MAP(), self.service_process, service_policy=self.service_policy)
        queue.M = queue.arrival_process.M
        return queue
class ExponentialServiceProcess(ServiceProcess):...

class MarkovianArrivalProcess(ArrivalProcess):
    def __init__(self,D0:np.ndarray,D1:np.ndarray):
        """MarkovianArrivalProcess

        Args:
            D0 (ndarray): matC, the transmission matrix
            D1 (ndarray): matD, the emission matrix
        """
        matC=np.asarray(D0)
        matD=np.asarray(D1)
        if   matC.ndim != 2 or matC.ndim != 2:
            raise ValueError(f"Matrices required. The matC provided is {self.matC.shape}-d, but the matD is {self.matD.shape}-d")
        elif matC.shape    != matD.shape:
            raise ValueError(f"Shape mismatch. the matC provided has shape {self.matC.shape}, but the matD provided has {self.matD.shape}")
        elif matC.shape[0] != matC.shape[1]:
            raise ValueError(f"Not square matrices. the matC provided has shape {self.matC.shape}, but the matD provided has {self.matD.shape}")
        else:
            self.M = self.dim = self.number_of_states = matC.shape[0]
            if not np.array_equal(matD, np.abs(matD)):
                warnings.warn(f"The emission matrix is not elementwise positive. Trying to modify it to make the system legal.")
                matD = np.abs(matD)
            if not np.array_equal(matC-np.diag(np.diag(matC)), np.abs(matC-np.diag(np.diag(matC)))):
                warnings.warn(f"The transmission matrix is not off-diagonal elementwise positive. Trying to modify it to make the system legal.")
                matC = np.diag(np.diag(matC))+np.abs(matC-np.diag(np.diag(matC)))
            if np.abs(rowsums:=np.sum(matC+matD,axis=1)).sum() != 0:
                warnings.warn(f"The generator has rowsums {rowsums} instead of 0's. Trying to modify the diagonal of the transmission matrix to make the system legal.",stacklevel=3)
                matC = matC-np.diag(rowsums)
            if (pi:=null_space(matC+matD)).shape != (self.M,1): raise ValueError(f"No unique steady solution for the markov chain; the generator has null_space {pi.tolist()}")
            self.D0=self.matC=matC
            self.D1=self.matD=matD
            self.π = self.pi     = self.arrival_steady_distribution = self.steady_distribution = pi.reshape((self.M,))/pi.sum()
            super(MarkovianArrivalProcess, self).__init__(rate=self.pi @ self.matD @ np.ones(self.M))
            self.λ = self.lbd = self.arrival_rate
            self.θ = self.theta  = np.amax(np.abs(self.matC))
            self._M_factorial = factorial(self.M)

    @classmethod # trick: use class method to implement method overloading
    def from_p(cls, p:float):
        """Use the original example used by the paper

        Args:
            p (float): an adjusting parameter for the emission matrix
            matC=np.array([[-2  ,      0],[0        ,-0.5  ]])
            matD=np.array([[ 2*p,2*(1-p)],[0.5*(1-p), 0.5*p]])
        """
        if p<0 or p>1:
            raise ValueError(f"the provided parameter p={p} does not lie within [0,1]")
        matC=np.array([[-2  , 0      ],[ 0        ,-0.5  ]])
        matD=np.array([[ 2*p, 2*(1-p)],[ 0.5*(1-p), 0.5*p]])
        queue = cls(matC, matD)
        return queue

    def   dist(self, x:float)->float: return                   self.π.dot(self.D0).dot( expm(x*self.D0)).dot(np.ones(self.M))
    def moment(self, k:int  )->float: return self._M_factorial*self.π.dot(matrix_power(-inv(self.D0),k)).dot(np.ones(self.M))
MAP = MarkovianArrivalProcess

class PoissonArrivalProcess(ArrivalProcess):
    def __init__(self, rate:float=1):
        super(PoissonArrivalProcess, self).__init__(rate)
        self.dist = poisson(rate)
    def as_MAP(self)->MarkovianArrivalProcess:
        return MarkovianArrivalProcess(np.array([-self.arrival_rate]),np.array([self.arrival_rate]))

class PoissonArrivalExponentialServiceFIFOQueue(QueuingProcess):
    def __init__(self, arrival_rate:float=1, service_rate:float=1):
        super(PoissonArrivalExponentialServiceFIFOQueue, self).__init__(PoissonArrivalProcess(arrival_rate), ExponentialServiceProcess(service_rate))
MM1Queue = PoissonArrivalExponentialServiceFIFOQueue

class PoissonArrivalExponentialServiceProcessorSharingQueue(QueuingProcess):
    def __init__(self, arrival_rate:float=1, service_rate:float=1):
        super(PoissonArrivalExponentialServiceProcessorSharingQueue, self).__init__(PoissonArrivalProcess(arrival_rate), ExponentialServiceProcess(service_rate), service_policy="PS")
MM1PSQueue =  PoissonArrivalExponentialServiceProcessorSharingQueue

class MarvovianArrivalExponentialServiceFIFOQueue(QueuingProcess):
    """
    A MAP/M/1-FIFO queue
    """
    def __init__(self, MAP:MarkovianArrivalProcess, service_rate:float):
        super(MarvovianArrivalExponentialServiceFIFOQueue, self).__init__(MAP, ExponentialServiceProcess(service_rate))
        self.MAP = self.arrival_process
        self.M   = self.arrival_process.M
    def as_MAP(self)->QueuingProcess:
        return self
MAPM1Queue = MarvovianArrivalExponentialServiceFIFOQueue

class MarvovianArrivalExponentialServiceProcessorSharingQueue(QueuingProcess):
    """
    A MAP/M/1-PS queue
    """
    def _get_sojourn_time_matrix(self)->np.ndarray:
        return np.transpose(solve_qme(self.mu, np.transpose(self.MAP.matC-self.mu*np.eye(self.M)), np.transpose(self.MAP.matD), dim=self.M))
    def __init__(self, MAParr:MarkovianArrivalProcess, service_rate:float):
        super(MarvovianArrivalExponentialServiceProcessorSharingQueue, self).__init__(MAParr, ExponentialServiceProcess(service_rate), service_policy="PS")
        self.MAP = self.arrival_process
        self.M   = self.arrival_process.M
        self.R   = self.sojourn_time_matrix = self._get_sojourn_time_matrix()
        self.π0  = self.pi0 = self.empty_steady_distriution = self.MAP.pi.dot(np.eye(self.M)-self.R)
    def as_MAP(self)->QueuingProcess:
        return self
MAPM1PSQueue = MarvovianArrivalExponentialServiceProcessorSharingQueue
