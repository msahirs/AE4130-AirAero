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
NO_POINTS = 200 # Number of collocation points (N+1 Panels)
C_LEN = 1. # Unit chord used
# EPS = 0.1 * C_LEN # Use if thickness considered
V_INF = 1. # Unit vel used
RHO_INF = 1. # unit density used
ALPHA = 0  # in degrees

NACA = "4421" # For Task 2, you will need a reference airfoil. You can get your reference airfoil by adding up the three last digits of your student number

p = float(NACA[1])/10
m = float(NACA[0])/100

print(p)

alphas_deg = np.linspace(-5, 15, 21)  # from -5 to 15 degrees
alphas_rad = np.radians(alphas_deg)

# Initialize arrays for lift coefficient, moment coefficient, and lift gradient
C_L_values = []
C_m_c4_values = []

# Intermediate "Inputs"
ALPHA = math.radians(ALPHA)
Q_INF = 0.5 * RHO_INF * V_INF**2 #calc flow dyn pressure
U_INF = V_INF * math.cos(ALPHA) # induced horiz vel
W_INF = V_INF * math.sin(ALPHA) # induced vert vel

def get_NACA_camber(x,p,m):
    
    x/=C_LEN

    print(x)

    mask1 = (x <= p) & (x >= 0)
    mask2 = (x >= p) & (x <= 1)
    
    z = np.empty_like(x/C_LEN)

    z[mask1] = m / p**2 * (2 * p * x[mask1] - x[mask1]**2)
    z[mask2] = m / (1 - p)**2 * (1 - 2 * p + 2 * p * x[mask2] - x[mask2]**2)

    return z

# Plot plate for visuals (No. of points usage don't really matter for calculation of lift/pressure)
contour_x = np.array([(x/NO_POINTS/20) * C_LEN for x in range(NO_POINTS * 20 + 1)])
# contour_y = np.array([4 * EPS * (i/C_LEN) * (1-(i/C_LEN))  for i in contour_x]) # for non planar assumption
contour_y = get_NACA_camber(contour_x,p,m)

# Collocation Points
Xc = (np.arange(1,NO_POINTS+1) - 0.25) * C_LEN/NO_POINTS
# Xc  = [C_LEN/NO_POINTS * (i-0.25) for i in range(1,NO_POINTS+1)]
# Zc  = 4 * EPS * (Xc/C_LEN) * (1 - Xc/C_LEN) # Use this entry if non-planar assumption used
Zc  = get_NACA_camber(Xc,p,m)

# Vortex Points
X  = (np.arange(1,NO_POINTS+1) - 0.75) * C_LEN/NO_POINTS
# Z  = 4 * EPS * (X/C_LEN) * (1 - X/C_LEN) # Use this entry if non-planar assumption used
Z  = get_NACA_camber(X,p,m)

# Calculate the slope of the airfoil geometry
geom_slope = ((Zc-Z)/(Xc-X))[:,np.newaxis]

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

for alpha in alphas_rad:
    # Update flow conditions for the current angle of attack
    U_INF = V_INF * math.cos(alpha)
    W_INF = V_INF * math.sin(alpha)
    rhs_vel = np.vstack((np.ones(NO_POINTS) * U_INF, np.ones(NO_POINTS) * W_INF)).T
    rhs = -np.sum(rhs_vel * normals, axis=1)
    circs = np.linalg.solve(A, rhs)

    # Compute lift coefficient
    d_lift = RHO_INF * V_INF * circs
    C_L = np.sum(d_lift) / (Q_INF * C_LEN)
    C_L_values.append(C_L)

    # Compute quarter-chord moment coefficient
    moment_c4 = -np.sum(circs * (Xc - 0.25 * C_LEN)) * RHO_INF * V_INF / (Q_INF * C_LEN**2)
    C_m_c4_values.append(moment_c4)

# Estimate lift gradient (dC_L/dalpha) using linear fit
C_L_fit = np.polyfit(alphas_deg, C_L_values, 1)  # Linear fit
lift_gradient = C_L_fit[0]  # Slope of the fit

# Calculate lift distribution
d_lift = RHO_INF * V_INF * circs

# Compute element-wise pressure coefficient
d_cp = d_lift/ (C_LEN/NO_POINTS) / Q_INF

# Analytical pressure coefficient distribution
d_cp_A = 4 * np.sqrt((C_LEN - contour_x)/contour_x) * alphas_rad[-1] # +  \
        # 32 * (EPS/C_LEN) * np.sqrt((contour_x/C_LEN) * (1 - contour_x/C_LEN))

# PLOTTTING ADDITIONS

#Plot geometry and control/vortex points
plt.figure()
plt.plot(contour_x,contour_y,"k-", label = "Planar Airfoil Contour",alpha = 0.8)
plt.plot(Xc,Zc,"ro", label = "Control Points", markersize = 4)
plt.plot(X,Z,"gx", label = "Vortex Points", markersize = 5)
plt.title(f"Geometrical Discretisation, N = {NO_POINTS}")
plt.legend()
plt.ylim(-0.05,0.2)
plt.xlim(-0.1,1.1)
plt.grid()
ax = plt.gca()
# ax.set_aspect('equal')
plt.tight_layout()

plt.savefig(f"assignment_1/plots/disc_N_{NO_POINTS}_AoA_{math.degrees(ALPHA):1.0f}.pdf")


# Plot cp and cl
plt.figure()
plt.plot(X,d_cp,"bo-", label = "Numerical $\Delta C_p$", markersize = 3)
plt.plot(contour_x,d_cp_A,"r-", label = "Analytical $\Delta C_p$", markersize = 1, alpha = 0.75)
plt.xlabel(r"Chord Length $[x/c]$")
plt.ylabel(r"Element-Wise $\Delta C_p$")
# plt.title(f"Numerical $C_L$: {np.sum(d_lift)/(Q_INF * C_LEN):.4f} vs. Analytical: {2*math.pi * (ALPHA + 2 * EPS/C_LEN):.4f} at $Re_\infty$ = {1}")
plt.title(f"Numerical $C_L$: {np.sum(d_lift)/(Q_INF * C_LEN):.4f} vs. Analytical plate: {2*math.pi * (alphas_rad[-1]):.4f}, AoA = {math.degrees(ALPHA):1.0f}")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f"assignment_1/plots/cp_N_{NO_POINTS}_AoA_{math.degrees(alphas_rad[-1]):1.0f}.pdf")

# plt.show()

# Plot Lift Coefficient vs. Angle of Attack
plt.figure()
plt.plot(alphas_deg, C_L_values, "bo-", label="Numerical $C_L$", markersize=4)
plt.plot(alphas_deg, 2 * math.pi * np.radians(alphas_deg), "r-", label="Flate plate $C_L$", alpha=0.7)
plt.xlabel(r"Angle of Attack $\alpha$ [degrees]")
plt.ylabel(r"Lift Coefficient $C_L$")
plt.title(f"Lift Coefficient vs. Angle of Attack")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f"assignment_1/plots/lift_coefficient_vs_alpha.pdf")

# Plot Quarter-Chord Moment Coefficient vs. Angle of Attack
plt.figure()
plt.plot(alphas_deg, C_m_c4_values, "go-", label="Numerical $C_{m,c/4}$", markersize=4)
plt.xlabel(r"Angle of Attack $\alpha$ [degrees]")
plt.ylabel(r"Quarter-Chord Moment Coefficient $C_{m,c/4}$")
plt.title(f"Quarter-Chord Moment Coefficient vs. Angle of Attack")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f"assignment_1/plots/moment_coefficient_vs_alpha.pdf")

# Plot Lift Gradient vs. Alpha (indirectly as slope annotation)
plt.figure()
plt.plot(alphas_deg, C_L_values, "bo-", label="Numerical $C_L$", markersize=4)
plt.plot(alphas_deg, np.polyval(C_L_fit, alphas_deg), "r--", label=f"Linear Fit: $dC_L/d\\alpha = {lift_gradient*180/np.pi:.3f}$", alpha=0.7)
plt.xlabel(r"Angle of Attack $\alpha$ [degrees]")
plt.ylabel(r"Lift Coefficient $C_L$")
plt.title("Lift Gradient Estimation")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f"assignment_1/plots/lift_gradient_vs_alpha.pdf")