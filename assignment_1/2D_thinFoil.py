import numpy as np
from matplotlib import pyplot as plt


NO_POINTS = 20
C_LEN = 1.
EPS = 0.1 * C_LEN

Xc  = [C_LEN/NO_POINTS * (i-0.25) for i in range(1,NO_POINTS+1)]
Zc  = [0 * EPS * Xc[i]/C_LEN * (1 - Xc[i]/C_LEN) for i in range(NO_POINTS)]

X  = [C_LEN/NO_POINTS * (i-0.75) for i in range(1,NO_POINTS+1)]
Z  = [0 * EPS * X[i]/C_LEN * (1 - X[i]/C_LEN) for i in range(NO_POINTS)]

plt.plot(Xc,Zc,"-ro",)
plt.plot(X,Z,"-go",)
plt.xlim(0,1)
plt.ylim(0,1)
plt.show()