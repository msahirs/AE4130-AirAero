import numpy as np
import math
from matplotlib import pyplot as plt

# Base Inputs
NO_POINTS = 10
C_LEN = 1
# EPS = 0.1 * C_LEN # Use if thickness considered
V_INF = 1
RHO_INF = 1
ALPHA = 0  # in degrees

# Intermediate "Inputs"
ALPHA/=math.pi
Q_INF = 0.5 * RHO_INF * V_INF**2
U_INF = V_INF * math.cos(ALPHA)
W_INF = V_INF * math.sin(ALPHA)


contour_x = [(x/NO_POINTS/10) * C_LEN for x in range(NO_POINTS*10 + 1)]
contour_y = [0 for _ in range(NO_POINTS*10 + 1)]

# Collocation Points
Xc  = [C_LEN/NO_POINTS * (i-0.25) for i in range(1,NO_POINTS+1)]
# Zc  = [0 * EPS * Xc[i]/C_LEN * (1 - Xc[i]/C_LEN) for i in range(NO_POINTS)] # Use this entry if non-planar assumption used
Zc  = [0 for _ in range(NO_POINTS)]

# Vortex Points
X  = [C_LEN/NO_POINTS * (i-0.75) for i in range(1,NO_POINTS+1)]
# Z  = [0 * EPS * X[i]/C_LEN * (1 - X[i]/C_LEN) for i in range(NO_POINTS)] # Use this entry if non-planar assumption used
Z  = [0 for _ in range(NO_POINTS)]

plt.plot(contour_x,contour_y,"k-", label = "Flat Plate")
plt.plot(Xc,Zc,"ro", label = "Control Points", markersize = 3)
plt.plot(X,Z,"go", label = "Vortex Points", markersize = 3)


plt.xlim(-0.5, 1.5)
plt.ylim(-1, 1)
plt.legend()

plt.show()