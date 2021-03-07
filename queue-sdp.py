"""
    This script contains three parts:
    - temporal iteration method
    - matrix-geometrical method
    - parameter generation for the corresponding SDP formulation.
        The generation of SDP parameters follows the notations of (Bertsima et Natarajan, 2007).
    The script is only partially OOP; refinements should be done in the future.
"""
"""
    By default this script uses the following example:
    A two-queue FQ system with MAP arrival

    - Arrival:
        Using a two-node hidden markov chain to generate arrivals.
        Transmission matrix (D0):
        [[-2.00, 0.  ],
            [ 0.  ,-0.50]]
        Emission matrix     (D1):
        [[ 0.60, 1.40],
            [ 0.35, 0.15]]
        The coefficients are chosen arbitrarily.
        Note that the model is in continuous time, so the matrices are infinitesimal generators;
        To acquire matrices for temporal iteration, one should use exp(M δt) or (1+M δt), where δt is a suitably chosen segment of time.

        Generated arrivals are send to either queue with equal probability.
    - Service:
        The server have a simple exponential distribution of service time;
        The queues are served with equal probability (fair-queuing, aka weighed-fair-queuing with even chances)
        The mean service rate is 5.

        Note that this rate surpasses D1, and the system is positive-recurrent.
    - SDP:
    For a certain costumer at the k-th position in a queue, his waiting time equals the total service time of all precedent costumers;
    As long as service times are idd. and its PDF is f(t), his merely f(kt)
    For a newly arrived costumer, his waiting time have the distribution ∑_{k>0} P(queue length=k)f(kt)
    Hence, it is easy to move between the waiting-time distribution and the queue-length distribution.

    Following the notions in Bertsima 2007, the system have the following characteristics:
    -- the concerning stochastic variables are (W1, X1, W2, X2), where W1, W2 are waiting times of the two queues,
        and (X1, X2) is determined by the arrival and service distributions.
        Note that the MAP arrival is not idd., however in the steady regime temporal auto-correlation does not matter;
        thus only the steady distribution of arrivals is used.

    -- the decision-variables are
        x(α1,β1,α2,β2;k)=E(W1^α1,X1^β1,W2^α2,X2^β2, in the k-th geometrical region S_k)
    -- there are four geometrical regions:
        A. queue 1 empty, queue 2 empty  ↔  W1≥0, W1+X1≤0, W2≥0, W2+X2≤0    g(W1, W2) = (  0  ,   0  )
        B. queue 1 empty, queue 2  busy  ↔  W1≥0, W1+X1≤0, W2≥0, W2+X2≥0    g(W1, W2) = (  0  , W2+X2)
        C. queue 1  busy, queue 2 empty  ↔  W1≥0, W1+X1≥0, W2≥0, W2+X2≤0    g(W1, W2) = (W1+X1,   0  )
        D. queue 1  busy, queue 2  busy  ↔  W1≥0, W1+X1≥0, W2≥0, W2+X2≥0    g(W1, W2) = (W1+X1, W2+X2)

    with s(W1, X1, W2, X2) = W1 + W2, where {(W1, X1, W2, X2)|s(W1, X1, W2, X2)>0} contains the set {(W1, X1, W2, X2)|(W1, X1, W2, X2)>0}
    and r = 2, d = 1

    M_r = [[(0,0,0,0), (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1), (2,0,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1), (0,2,0,0), (0,1,1,0), (0,1,0,1), (0,0,2,0), (0,0,1,1), (0,0,0,2)],
           [(1,0,0,0), (2,0,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1), (3,0,0,0), (2,1,0,0), (2,0,1,0), (2,0,0,1), (1,2,0,0), (1,1,1,0), (1,1,0,1), (1,0,2,0), (1,0,1,1), (1,0,0,2)],
           [(0,1,0,0), (1,1,0,0), (0,2,0,0), (0,1,1,0), (0,1,0,1), (2,1,0,0), (1,2,0,0), (1,1,1,0), (1,1,0,1), (0,3,0,0), (0,2,1,0), (0,2,0,1), (0,1,2,0), (0,1,1,1), (0,1,0,2)],
           [(0,0,1,0), (1,0,1,0), (0,1,1,0), (0,0,2,0), (0,0,1,1), (2,0,1,0), (1,1,1,0), (1,0,2,0), (1,0,1,1), (0,2,1,0), (0,1,2,0), (0,1,1,1), (0,0,3,0), (0,0,2,1), (0,0,1,2)],
           [(0,0,0,1), (1,0,0,1), (0,1,0,1), (0,0,1,1), (0,0,0,2), (2,0,0,1), (1,1,0,1), (1,0,1,1), (1,0,0,2), (0,2,0,1), (0,1,1,1), (0,1,0,2), (0,0,2,1), (0,0,1,2), (0,0,0,3)],
           [(2,0,0,0), (3,0,0,0), (2,1,0,0), (2,0,1,0), (2,0,0,1), (4,0,0,0), (3,1,0,0), (3,0,1,0), (3,0,0,1), (2,2,0,0), (2,1,1,0), (2,1,0,1), (2,0,2,0), (2,0,1,1), (2,0,0,2)],
           [(1,1,0,0), (2,1,0,0), (1,2,0,0), (1,1,1,0), (1,1,0,1), (3,1,0,0), (2,2,0,0), (2,1,1,0), (2,1,0,1), (1,3,0,0), (1,2,1,0), (1,2,0,1), (1,1,2,0), (1,1,1,1), (1,1,0,2)],
           [(1,0,1,0), (2,0,1,0), (1,1,1,0), (1,0,2,0), (1,0,1,1), (3,0,1,0), (2,1,1,0), (2,0,2,0), (2,0,1,1), (1,2,1,0), (1,1,2,0), (1,1,1,1), (1,0,3,0), (1,0,2,1), (1,0,1,2)],
           [(1,0,0,1), (2,0,0,1), (1,1,0,1), (1,0,1,1), (1,0,0,2), (3,0,0,1), (2,1,0,1), (2,0,1,1), (2,0,0,2), (1,2,0,1), (1,1,1,1), (1,1,0,2), (1,0,2,1), (1,0,1,2), (1,0,0,3)],
           [(0,2,0,0), (1,2,0,0), (0,3,0,0), (0,2,1,0), (0,2,0,1), (2,2,0,0), (1,3,0,0), (1,2,1,0), (1,2,0,1), (0,4,0,0), (0,3,1,0), (0,3,0,1), (0,2,2,0), (0,2,1,1), (0,2,0,2)],
           [(0,1,1,0), (1,1,1,0), (0,2,1,0), (0,1,2,0), (0,1,1,1), (2,1,1,0), (1,2,1,0), (1,1,2,0), (1,1,1,1), (0,3,1,0), (0,2,2,0), (0,2,1,1), (0,1,3,0), (0,1,2,1), (0,1,1,2)],
           [(0,1,0,1), (1,1,0,1), (0,2,0,1), (0,1,1,1), (0,1,0,2), (2,1,0,1), (1,2,0,1), (1,1,1,1), (1,1,0,2), (0,3,0,1), (0,2,1,1), (0,2,0,2), (0,1,2,1), (0,1,1,2), (0,1,0,3)],
           [(0,0,2,0), (1,0,2,0), (0,1,2,0), (0,0,3,0), (0,0,2,1), (2,0,2,0), (1,1,2,0), (1,0,3,0), (1,0,2,1), (0,2,2,0), (0,1,3,0), (0,1,2,1), (0,0,4,0), (0,0,3,1), (0,0,2,2)],
           [(0,0,1,1), (1,0,1,1), (0,1,1,1), (0,0,2,1), (0,0,1,2), (2,0,1,1), (1,1,1,1), (1,0,2,1), (1,0,1,2), (0,2,1,1), (0,1,2,1), (0,1,1,2), (0,0,3,1), (0,0,2,2), (0,0,1,3)],
           [(0,0,0,2), (1,0,0,2), (0,1,0,2), (0,0,1,2), (0,0,0,3), (2,0,0,2), (1,1,0,2), (1,0,1,2), (1,0,0,3), (0,2,0,2), (0,1,1,2), (0,1,0,3), (0,0,2,2), (0,0,1,3), (0,0,0,4)]]

    We have

    M_{r-d}
        = [[(1,0,0,0)+(0,1,0,0), (2,0,0,0)+(1,1,0,0), (1,1,0,0)+(0,2,0,0), (1,0,1,0)+(0,1,1,0), (1,0,0,1)+(0,1,0,1)],
           [(2,0,0,0)+(1,1,0,0), (3,0,0,0)+(2,1,0,0), (2,1,0,0)+(1,2,0,0), (2,0,1,0)+(1,1,1,0), (2,0,0,1)+(1,1,0,1)],
           [(1,1,0,0)+(0,2,0,0), (2,1,0,0)+(1,2,0,0), (1,2,0,0)+(0,3,0,0), (1,1,1,0)+(0,2,1,0), (1,1,0,1)+(0,2,0,1)],
           [(1,0,1,0)+(0,1,1,0), (2,0,1,0)+(1,1,1,0), (1,1,1,0)+(0,2,1,0), (1,0,2,0)+(0,1,2,0), (1,0,1,1)+(0,1,1,1)],
           [(1,0,0,1)+(0,1,0,1), (2,0,0,1)+(1,1,0,1), (1,1,0,1)+(0,2,0,1), (1,0,1,1)+(0,1,1,1), (1,0,0,2)+(0,1,0,2)]]
"""

from   math import prod
from   itertools import combinations, chain
from   typing import Tuple
import numpy as np
from   scipy.special import comb
import torch
'''
    cuda = True if torch.cuda.is_available() else False
    if cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)
        TORCH = torch.cuda
        Tensor     = (lambda x: torch.cuda.Tensor(     x).to(torch.device('cuda:0')))
        ByteTensor = (lambda x: torch.cuda.ByteTensor( x).to(torch.device('cuda:0')))
        FloatTensor= (lambda x: torch.cuda.FloatTensor(x).to(torch.device('cuda:0')))
    else:
        TORCH = torch
        Tensor     = torch.Tensor
        ByteTensor = torch.ByteTensor
        FloatTensor= torch.FloatTensor
'''
from reprod_masuyama2003 import MarkovianArrivalProcess, ExponentialServiceProcess
class RandomVariable:
    def dist_func():    raise NotImplementedError
    def moment():       raise NotImplementedError

def to_tuple(lst:list)->Tuple:
    return tuple(to_tuple(i) if isinstance(i, list) else i for i in lst)

def Partition(lst:list, n:int, d:int=0):
    """MMA like list partition generator

    Args:
        lst (list): list to be partitioned
        n (int): length of sublists
        d (int, optional): offset. Defaults to n.

    Yields:
        list: sublist of length n, offset from the previous one by d.
    """
    d = n if d==0 else d
    for i in range(0, len(lst)-d, d):
            yield lst[i:i + n]
def IntegerCompositions(n:int, k:int):
    return [[y[1]-y[0]-1 for y in Partition([0]+list(x)+[n+k],2,1)] for x in combinations(range(1,n+k),k-1)]
class INTegerCompositions:
    def __init__(self): self._integerCompositions = {}
    def __call__(self, n:int, k:int):
        if (n,k) not in self._integerCompositions: self._integerCompositions[(n,k)] = list(reversed(IntegerCompositions(n,k)))
        return self._integerCompositions[(n,k)]
integerCompositions=INTegerCompositions()
def edg2mat(lst:list)->list:
    return [[x+y for y in lst] for x in lst]
def edg2mat(lst:np.ndarray)->np.ndarray:
    return np.array([lst+x for x in lst])
'''
    def edg2mat(lst:Tensor)->Tensor:
        return Tensor([lst+x for x in lst])
'''

D0 = np.array([[-2.00, 0.  ],
               [ 0.  ,-0.50]])
D1 = np.array([[ 0.60, 1.40],
               [ 0.35, 0.15]])

mu = 5
K  = 2 #* 2 queues
r  = 2 #* rank of relaxation


arrival =   MarkovianArrivalProcess(D0, D1) #* provides the moments of the variable T
service = ExponentialServiceProcess(mu)     #* provides the moments of the variable S
inprob  = np.ones(K)/K #* input evenly classified

def m(β:np.ndarray)->np.ndarray:
    """Calculating the m^β as defined in Bertsima 2007
    m^β := E[X1^β1*X2^β2] = E[X1^β1]*E[X2^β2]

    Returns:i
        np.ndarray
    """
    return arrival.moment(β[0])*arrival.moment(β[1])/(inprob[0]**(β[0]+1))/(inprob[1]**(β[1]+1))

def g(occupation:Tuple[bool, ...], α_out:np.ndarray, α_in:np.ndarray, β_in:np.ndarray)->int:
    """Expansion coefficient of E[g(W+X)]^{α,β} in terms of E[(W+X)]^{α,β}
        g_{k,α_out}^{α_in,β_in}

    Parameters:
        occupation: a tuple of booleans, showing whether each queue is empty.
        α, β: multi-indices

    Returns:
        float
    """
    if α_in + β_in == α_out:
        ans = 1
        for n, oc in enumerate(occupation):
            if oc == False: continue
            else: ans = ans*comb(α_out[n], α_in[n])
    else: ans=0
    return ans

'''
_η={(0,0): np.array([[0, 0],   #* α
                     [0, 0]])} #* β
def η(i:int, j:int)->Tuple[np.ndarray,np.ndarray]:
    """η(i, j) superscript in M_r(x)
    M_r[i][j] = x^η(i, j)

    Returns:
        np.ndarray: (2, 0) array providing the multi-index
    """
    if (i,j) in _η: return _η[(i,j)]
    #// elif i==0 and j==0: return (np.array([0, 0]), np.array([0, 0]))
    elif i> 0 and j==0:
        #// n = int(sqrt(8*i**2+1)/2)
        #// _η[(i,j)] = (np.array([(1+n)*n/2-i, (1-n)*n/2+i]), np.array([(1+n)*n/2-i, (1-n)*n/2+i]))
        if _η(i-1, 0)[0] == np.array([0, 0]) and _η(i-1, 1)[0][1] == 0:
            pass
    elif i==0 and j> 0: _η[(i,j)] = η(j, 0)
    elif i> 0 and j> i: _η[(i,j)] = η(i, 0) + η(j, 0)
    return _η[(i,j)]

def s(α:np.ndarray, β:np.ndarray)->int:
    """s_γ
    Using s(W1, X1, W2, X2) = W1 + W2
    alternatives:
        s(W1, X1, W2, X2) = (W1**2)*W2 + W1*(W2**2)
        s(W1, X1, W2, X2) = W1 + W2 + (W1**2)*W2 + W1*(W2**2)
    Any polynomial of (W1, X1, W2, X2) will do, as long as it is positive over {(W1, W2)| W1>0 && W2>0}.

    Args:
        α (np.ndarray)
        β (np.ndarray)

    Returns:
        int: the coefficient of the corresponding γ in s
    """
    ans = 0
    for k in α:
        if   ans == 0 and k == 1: ans = 1
        elif ans == 1 and k != 0: return 0
    return ans
# s = dict([((np.eye(2*K)[k],np.zeros(K)), 1) for k in range(K)]) # an equivalent but somewhat better way

def localized_matrix_element(i:int, j:int)->Dict[np.ndarray, float]:
    """Localized M_{r-d} using polynomial s

    Returns:
        Dict[np.ndarray, float]: multi-indices as keys pointing to their coefficients as values
    """
    ...
'''

edgIC = np.array(list(chain.from_iterable([[np.array([ic[:K],ic[K:]],np.int16) for ic in integerCompositions(rr, 2*K)] for rr in range(  r+1)])))
matIC = edg2mat(edgIC)
matICdim = (len(matIC),len(matIC[0]))
αIC   = np.array(list(chain.from_iterable([[np.array( ic[:K]        ,np.int16) for ic in integerCompositions(rr,   K)] for rr in range(2*r+1)])))

class LocInMatIC:
    def __init__(self, key:str='byte'): self._key = key; self._loc = {}
    def __call__(self, multi_index):
        if   self._key == 'byte':
            _index = multi_index.tobytes()
        elif self._key == 'tuple':
            _index = to_tuple(multi_index.tolist())
        else: raise NotImplementedError
        if _index not in self._loc:
            for nrow, row in enumerate(matIC):
                for ncol, some_index in enumerate(row[nrow:]):
                    if not (multi_index-some_index).any():
                        self._loc[_index] = (nrow, ncol)
                        break
        try: return self._loc[_index]
        except KeyError:
            print("existing locs:\n", self._loc)
            raise
loc=LocInMatIC()

indieA_i1  = torch.IntTensor([list((mm,)+loc(multi_index                )) for mm in range(2**K) for multi_index in chain.from_iterable([row[ind:] for ind,row in enumerate(matIC)])])
indieA_v1  = torch.Tensor([              1                                 for _  in range(2**K) for      _      in chain.from_iterable([row[ind:] for ind,row in enumerate(matIC)])])
indieA_im  = torch.IntTensor([list((mm,)+loc(multi_index*np.array([1,0]))) for mm in range(2**K) for multi_index in chain.from_iterable([row[ind:] for ind,row in enumerate(matIC)])])
indieA_vm  = torch.Tensor([m(multi_index[1]) for _ in range(2**K) for multi_index in chain.from_iterable([row[ind:] for ind,row in enumerate(matIC)])])
indieA     = torch.sparse_coo_tensor(torch.cat((indieA_i1,indieA_im)).t(),torch.cat((indieA_v1,indieA_vm)),(2**K,)+matICdim)
print(indieA)
# from operator import mul
from sympy import Poly, symbols
combieA_dicts = [Poly(prod(sum(symbols(f'w{k},x{k}'))**idx for k, idx in enumerate(tuple(multi_index))), symbols(f'w:{K},x:{K}')).as_dict() for multi_index in αIC]
combieA_i1 = [[loc(np.array([list(idx[:K]),list(idx[K:])])) for idx in dict_] for _, dict_ in enumerate(combieA_dicts)]
print(combieA_i1)
# combieA_i1 = torch.IntTensor([Poly(prod(sum(symbols(f'w{k},x{k}'))**idx for k, idx in enumerate(tuple(multi_index))), symbols(f'w:{K},x:{K}')).as_dict()
#                                                                            for multi_index in chain.from_iterable([row[ind:] for ind,row in enumerate(matIC)])])
# print(combieA_i1)
