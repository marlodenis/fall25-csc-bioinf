import numpy as np

n = 7
p = [0.00, 0.04, 0.06, 0.08, 0.02, 0.10, 0.12, 0.14]
q = [0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05]

e = np.zeros([9, 8])
w = np.zeros([9, 8])
root = np.zeros([8, 8])
for i in range(1, 9):
	e[i, i - 1] = q[i - 1]
	w[i, i - 1] = q[i - 1]
for l in range(1, 8):
	for i in range(1, 9 - l):
		j = i + l - 1
		e[i, j] = 100000000
		w[i, j] = w[i, j - 1] + p[j]+ q[j]
		for r in range(i, j + 1):
			t = e[i, r - 1] + e[r + 1, j] + w[i, j]
			if t < e[i, j]:
				e[i, j] = t
				root[i, j] = r

print(e)
print(w)
print(root)
