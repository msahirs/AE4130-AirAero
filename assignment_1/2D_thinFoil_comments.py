import numpy as np
import math
from matplotlib import pyplot as plt


# Curl matrix for 2D vortex calculations
curl_mat = np.array([[0 , 1],
        [-1, 0]])

def vort_mag_2d(x,z,x_j,z_j):
    """
    Calculate the magnitude of the vortex velocity at a point (x, z) due to a vortex located at (x_j, z_j).
    """
    r_sq = (x-x_j)**2 + (z-z_j)**2

    disp = np.array([[x-x_j],
                     [z-z_j]])
    
    return 1/(2*np.pi * r_sq) * curl_mat @ disp

# Base Inputs
NO_POINTS = 50 # Number of collocation points (N+1 Panels)
C_LEN = 1. # UNit chord used
EPS = 0.1 * C_LEN # Use if thickness considered
V_INF = 1. # Unit vel used
RHO_INF = 1. # unit density used
ALPHA = 10  # in degrees

# Intermediate "Inputs"
ALPHA = math.radians(ALPHA)
Q_INF = 0.5 * RHO_INF * V_INF**2 #calc flow dyn pressure
U_INF = V_INF * math.cos(ALPHA) # induced horiz vel
W_INF = V_INF * math.sin(ALPHA) # induced vert vel

# Plot plate for visuals (No. of points usage don't really matter for calculation of lift/pressure)
contour_x = np.array([(x/NO_POINTS/40) * C_LEN for x in range(NO_POINTS * 40 + 1)])
# contour_y = np.array([4 * EPS * (i/C_LEN) * (1-(i/C_LEN))  for i in contour_x]) # for non planar assumption
contour_y = np.zeros_like(contour_x)

# Collocation Points
# Xc = (np.arange(1,NO_POINTS+1) - 0.25) * C_LEN/NO_POINTS
Xc  = [C_LEN/NO_POINTS * (i-0.25) for i in range(1,NO_POINTS+1)]
# Zc  = 4 * EPS * (Xc/C_LEN) * (1 - Xc/C_LEN) # Use this entry if non-planar assumption used
Zc  = np.zeros(NO_POINTS)

# Vortex Points
X  = (np.arange(1,NO_POINTS+1) - 0.75) * C_LEN/NO_POINTS
# Z  = 4 * EPS * (X/C_LEN) * (1 - X/C_LEN) # Use this entry if non-planar assumption used
Z  = np.zeros(NO_POINTS)

# Calculate the slope of the airfoil geometry
geom_slope = ((Zc-Z)/(Xc-X)) [:,np.newaxis]

# Compute normals and tangents
normals = np.hstack((-geom_slope, np.ones((NO_POINTS,1)))) / np.sqrt(geom_slope**2 + 1)
tangents = np.hstack((normals[:,1], -normals[:,0]))

# Initialize the coefficient matrix A
A = np.zeros((NO_POINTS,NO_POINTS))

# Populate matrix A with influence coefficents
for i in range(NO_POINTS):
    for j in range(NO_POINTS):
        vel = vort_mag_2d(Xc[i],Zc[i],X[j],Z[j])
        term = np.dot(vel[:,0],normals[i,:])
        A[i,j] = term

# Set up the right-hand side (RHS) velocity vector
rhs_vel = np.vstack((np.ones(NO_POINTS)*U_INF, np.ones(NO_POINTS)*W_INF)).T

# Compute the RHS of the linear system
rhs =  - np.sum(rhs_vel * normals, axis = 1)

# Solve for circulation strengths
circs = np.linalg.solve(A,rhs)

# Calculate lift distribution
d_lift = RHO_INF * V_INF * circs

# Compute element-wise pressure coefficient
d_cp = d_lift/ (C_LEN/NO_POINTS) / Q_INF

# Analytical pressure coefficient distribution
d_cp_A = 4 * np.sqrt((C_LEN - contour_x)/contour_x) * ALPHA # +  \
        # 32 * (EPS/C_LEN) * np.sqrt((contour_x/C_LEN) * (1 - contour_x/C_LEN))

# PLOTTTING ADDITIONS

#Plot geometry and control/vortex points
plt.figure()
plt.plot(contour_x,contour_y,"k-", label = "Planar Airfoil Contour")
plt.plot(Xc,Zc,"ro", label = "Control Points", markersize = 2)
plt.plot(X,Z,"gx", label = "Vortex Points", markersize = 4)
plt.title(f"Geometrical Discretisation, N = {NO_POINTS}")
plt.legend()

plt.savefig(f"assignment_1/plots/disc_N_{NO_POINTS}_AoA_{math.degrees(ALPHA):1.0f}.pdf")


# Plot cp and cl

plt.figure()
plt.plot(X,d_cp,"bo-", label = "Numerical $\Delta C_p$", markersize = 3)
plt.plot(contour_x,d_cp_A,"r-", label = "Analytical $\Delta C_p$", markersize = 1, alpha = 0.75)


plt.xlabel(r"Chord Length $[x/c]$")
plt.ylabel(r"Element-Wise $\Delta C_p$")

# plt.title(f"Numerical $C_L$: {np.sum(d_lift)/(Q_INF * C_LEN):.4f} vs. Analytical: {2*math.pi * (ALPHA + 2 * EPS/C_LEN):.4f} at $Re_\infty$ = {1}")
plt.title(f"Numerical $C_L$: {np.sum(d_lift)/(Q_INF * C_LEN):.4f} vs. Analytical: {2*math.pi * (ALPHA):.4f}, $Re_\infty$ = {1}, AoA = {math.degrees(ALPHA):1.0f}")

# plt.xlim(-0.5, 1.5)
# plt.ylim(0, 5)
plt.legend()
plt.savefig(f"assignment_1/plots/cp_N_{NO_POINTS}_AoA_{math.degrees(ALPHA):1.0f}.pdf")

# plt.show()