import numpy as np
import math
from matplotlib import pyplot as plt


def vort_mag_2d(x,z,x_j,z_j):
    
    r_sq = (x-x_j)**2 + (z-z_j)**2

    curl_mat = np.array([[0 , 1],
              [-1, 0]])
    disp = np.array([[x-x_j],
                     [z-z_j]])
    
    return 1/(2*np.pi * r_sq) * curl_mat @ disp
# Base Inputs
NO_POINTS = 5
C_LEN = 1.
# EPS = 0.1 * C_LEN # Use if thickness considered
V_INF = 1.
RHO_INF = 1.
ALPHA = 5.  # in degrees

# Intermediate "Inputs"
ALPHA/=math.pi
Q_INF = 0.5 * RHO_INF * V_INF**2
U_INF = V_INF * math.cos(ALPHA)
W_INF = V_INF * math.sin(ALPHA)

# Plot plate for visuals (No. of points usage don't really matter for calculation of lift/pressure)
contour_x = np.array([(x/NO_POINTS/10) * C_LEN for x in range(NO_POINTS * 10 + 1)])
contour_y = np.array([0 for _ in range(NO_POINTS*10 + 1)])

# Collocation Points
Xc = (np.arange(1,NO_POINTS+1) - 0.25) * C_LEN/NO_POINTS
# Xc  = [C_LEN/NO_POINTS * (i-0.25) for i in range(1,NO_POINTS+1)]
# Zc  = [0 * EPS * Xc[i]/C_LEN * (1 - Xc[i]/C_LEN) for i in range(NO_POINTS)] # Use this entry if non-planar assumption used
Zc  = np.zeros(NO_POINTS)

# Vortex Points
X  = (np.arange(1,NO_POINTS+1) - 0.75) * C_LEN/NO_POINTS
# Z  = [0 * EPS * X[i]/C_LEN * (1 - X[i]/C_LEN) for i in range(NO_POINTS)] # Use this entry if non-planar assumption used
Z  = np.zeros(NO_POINTS)

geom_slope = ((Zc-Z)/(Xc-X)) [:,np.newaxis]

normals = np.hstack((-geom_slope, np.ones((NO_POINTS,1)))) / np.sqrt(geom_slope**2 + 1)
tangents = np.hstack((normals[:,1], -normals[:,0]))

A = np.zeros((NO_POINTS,NO_POINTS))

for i in range(NO_POINTS):
    for j in range(NO_POINTS):
        vel = vort_mag_2d(Xc[i],Zc[i],X[j],Z[j])
        term = np.dot(vel[:,0],normals[i,:])
        A[i,j] = term

rhs_vel = np.vstack((np.ones(NO_POINTS)*U_INF, np.ones(NO_POINTS)*W_INF)).T

rhs = - np.sum(rhs_vel * normals, axis = 1)

circs = np.linalg.solve(A,rhs)

d_lift = RHO_INF * V_INF * circs

d_c_p = d_lift/ (Q_INF * (C_LEN/NO_POINTS))

plt.plot(contour_x,contour_y,"k-", label = "Flat Plate")
plt.plot(Xc,Zc,"ro", label = "Control Points", markersize = 4)
plt.plot(X,Z,"go", label = "Vortex Points", markersize = 4)


plt.xlim(-0.5, 1.5)
plt.ylim(-1, 1)
plt.legend()

plt.show()