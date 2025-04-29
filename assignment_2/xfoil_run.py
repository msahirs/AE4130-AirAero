import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

#ENV Variables
XFOIL_PATH = "C:\\Rec\\UniStuff\\MSc\\AE4130-22_AircraftAero\\XFOIL6.99\\xfoil.exe"
plot_path = "assignment_2\\plots\\"
data_path = "assignment_2\\data\\"

def run_xfoil_2_4():

    try:
        os.remove(data_path+"polar_base.dat")
    except:
        next
    
    # commands to drive XFOIL
    cmds = [
        "NACA 2009",      # load NACA 2009
        "OPER",           # enter OPER menu
        "VPAR",           # set amplification factor
        "N 12",
        "",
        "VISC 0.7e6",       # set Reynolds number
        "MACH 0",         # set Mach number
        "PACC",           # start polar accumulation
        data_path+"polar_base.dat",      # output file
        "",               # no dump file
        "ASEQ -2 8 0.9999",    # α from -2 to 8 in steps of 1°
        "PACC",           # stop polar accumulation
        "",               # no polar‐file arg
        "QUIT"            # exit XFOIL
    ]
    
    # launch XFOIL and feed it all commands at once
    proc = subprocess.Popen(
        [XFOIL_PATH],
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    input_str = "\n".join(cmds) + "\n"
    # feed XFOIL and wait, encode as ASCII
    proc.communicate(input_str.encode('ascii'))

    # skip the 12‐line XFOIL header before the numeric table
    data = np.loadtxt(data_path+"polar_base.dat", skiprows=12)
    alpha, Cl, Cd = data[:,0], data[:,1], data[:,2]

    # plot Cl vs alpha
    plt.figure()
    plt.plot(alpha, Cl, "o-", lw=1.5)
    plt.xlabel("α (deg)")
    plt.ylabel("C_l")
    plt.title("Lift Curve for NACA 2009 (Re=0.7e6, VPAR n=12)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path+"Cl_vs_alpha.pdf", dpi=300)

    # plot Cd vs alpha
    plt.figure()
    plt.plot(alpha, Cd, "s-", color="r", lw=1.5)
    plt.xlabel("α (deg)")
    plt.ylabel("C_d")
    plt.title("Drag Curve for NACA 2009 (Re=0.7e6, VPAR n=12)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path+"Cd_vs_alpha.pdf", dpi=300)

    # plot drag polar (C_d vs C_l)
    plt.figure()
    plt.plot(Cl, Cd, "d-", color="m", lw=1.5)
    plt.xlabel("C_l")
    plt.ylabel("C_d")
    plt.title("Drag Polar for NACA 2009 (Re=0.7e6, VPAR n=12)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path+"Cd_vs_Cl.pdf", dpi=300)

    # plot transition point on top & bottom vs alpha
    plt.figure()
    plt.plot(alpha, data[:,5], '^-', lw=1.5, label="Top Xtr")
    plt.plot(alpha, data[:,6], 'v-', lw=1.5, label="Bot Xtr")
    plt.xlabel("α (deg)")
    plt.ylabel("Transition location x/c")
    plt.title("Boundary‐Layer Transition vs Angle of Attack")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path+"Xtr_vs_alpha.pdf", dpi=300)

    # clean up if desired
    # os.remove(data_path+"polar.dat")

if __name__ == "__main__":
    run_xfoil_2_4()