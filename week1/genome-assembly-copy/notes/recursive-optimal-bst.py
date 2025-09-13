import numpy as np

n = 7
p = [0.00, 0.04, 0.06, 0.08, 0.02, 0.10, 0.12, 0.14]
q = [0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05]

e = np.zeros((8,8))
w = np.zeros((8,8))
root = np.zeros((8,8))

def get_e(i, j):
    if j < i:
        return q[j]
    if e[i][j] != 0:
        return e[i][j]
    if i == j:
        root[i][j] = i
        e[i][j] = get_w(i,i) + q[i-1] + q[i]
        return e[i][j]
    m = 66666
    for p in range(j-i+1):
        x = get_e(i, i + p - 1) + get_e(i + p + 1, j) + get_w(i, j)
        if x < m:
            m = x
            root[i][j] = p + i
    e[i][j] = m
    return e[i][j]

def get_w(i, j):
    if j < i:
        return q[j]
    if w[i][j] != 0:
        return w[i][j]
    else:
        w_temp = get_w(i, j-1) + p[j] + q[j]
        w[i][j] = w_temp
        return w_temp

get_e(1, 7)
print(e)
print(w)
print(root)
