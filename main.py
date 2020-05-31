#test
from collections import Counter
import numpy as np
from random import shuffle
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
n = 1000
complexity = []

def bubble(arr):
    for i in range(len(arr) - 1):
        for j in range(len(arr) - i - 1):
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
    return arr

for i in tqdm(range(1, n + 1)):
    a = list(range(1, i+1))[::-1]
    s = time.time()
    a = bubble(a)
    complexity.append(time.time() - s)
complexity = np.array([100000*v  for v in complexity])

X = list(range(1, n+1))
X_ = np.array([[1, x, x**2] for x in X])
print(X_.shape)
print(complexity.shape)
m = Ridge().fit(X_, complexity)
x_plot = np.linspace(0, n, n)

y_pred = m.predict(X_)


plt.figure()
plt.plot(complexity, color='blue')
plt.plot(y_pred, color='red')
plt.xlabel('n')
plt.ylabel('time')
plt.show()
