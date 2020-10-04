import mpy4py
import random

def generate(sz, rg, sl):
    """
    :type int sz: size of square matrix(sz * sz)
    :type int rg, int sl: range and scale of data. For example if the value field will be [-0.2, 0.2], than
    rg = 2, sl = 10
    :rtype: (matrix List[List[float]], vector List[float]), float in Python equal to double in C

    Way for generate matrix in parallel is already shown in C version, I will not do that again in Python
    """
    random.seed(0)
    matrix = []
    vector = []
    for _ in range(sz):
        ls = []
        for i in range(sz):
            ele = (-rg + random.random() * 2 * rg) / sl
            ls.append(ele)
        matrix.append(ls)
        vector.append((-rg + random.random() * 2 * rg) / sl)
    return matrix, vector


def re_arrange(M, V):
    """
    :type M: List[List[float]] generated matrix
    :type V: Lists[float] generated vector. len(vector) == len(matrix)
    :rtype (M List[List[float]], V List[float]) rearranged matrix and vector. Ax = b => PAx = Pb
    """
    def find_mx(col):
        mx = 0
        idx = 0
        for i in range(col, len(M)):
            cur = abs(M[i][col])
            if cur > mx:
                mx = cur
                idx = i
        if mx == 0:
            print("Invalid Matrix")
            exit(0)
        return idx
    for i in range(M):
        target = find_mx(i)
        M[i], M[target] = M[target], M[i]
        V[i], V[target] = V[target], V[i]
    return M, V

def lu(M):
    

