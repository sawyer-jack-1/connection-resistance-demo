import numpy as np

def schur_comp(M, indices, indices_comp):
    A = M[np.ix_(indices, indices)]
    B = M[np.ix_(indices, indices_comp)]
    C = M[np.ix_(indices_comp, indices)]
    D = M[np.ix_(indices_comp, indices_comp)]
    Dinv = np.linalg.inv(D)
    return A - np.dot(B, np.dot(Dinv, C))

def nonzero_eigenvalue_product(A, thr=10e-2):
    eval, v = np.linalg.eig(A)

    idx = eval>= thr
    eval_nonzero = eval[eval>= thr]
    wtf = np.prod(eval_nonzero)
    return np.prod(eval_nonzero)