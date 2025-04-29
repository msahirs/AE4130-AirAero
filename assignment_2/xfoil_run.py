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
    plt.savefig(plot_path+"Cl_vs_alpha.pdf", dpi=150)

    # plot Cd vs alpha
    plt.figure()
    plt.plot(alpha, Cd, "s-", color="r", lw=1.5)
    plt.xlabel("α (deg)")
    plt.ylabel("C_d")
    plt.title("Drag Curve for NACA 2009 (Re=0.7e6, VPAR n=12)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path+"Cd_vs_alpha.pdf", dpi=150)

    # plot drag polar (C_d vs C_l)
    plt.figure()
    plt.plot(Cl, Cd, "d-", color="m", lw=1.5)
    plt.xlabel("C_l")
    plt.ylabel("C_d")
    plt.title("Drag Polar for NACA 2009 (Re=0.7e6, VPAR n=12)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path+"Cd_vs_Cl.pdf", dpi=150)

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
    plt.savefig(plot_path+"Xtr_vs_alpha.pdf", dpi=150)

    # clean up if desired
    # os.remove(data_path+"polar.dat")


def run_xfoil_fixed_transition_sweep():
    """
    Force the boundary‐layer transition to specified chordwise locations,
    run α from -2 to 8°, then plot:
      1) actual transition (should be flat at forced xtr)
      2) drag coefficient vs forced transition location (at α=0°)
    """
    # list of forced transition locations (x/c)
    data_dir = "transition_variation\\"
    
    xtr_list = [0.2, 0.4, 0.6, 0.8]
    cd0_list = []
    # collect (xtr, alpha, Cd) for combined plotting
    sweep_results = []

    for xtr in xtr_list:
        
        fname = f"{data_path + data_dir}polar_xtr_{xtr:.2f}.dat"
        try:
            os.remove(fname)
        except OSError:
            pass
        print(fname)
        # build XFOIL command deck
        cmds = [
            "NACA 2009",
            "OPER",
            "VPAR",                   # enter VPAR menu
            f"XTR {xtr:.2f} {xtr:.2f}",# force transition top & bottom
            "N 12",                   # amplification factor
            "",
            "VISC 0.7e6",
            "MACH 0",
            "PACC",
            fname,                    # polar output
            "",
            "ASEQ -2 8 0.9999",            # α from -2 to 8° in 1° steps
            "PACC",
            "",
            "QUIT"
        ]

        # launch XFOIL
        proc = subprocess.Popen(
            [XFOIL_PATH],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        input_str = "\n".join(cmds) + "\n"
        # feed XFOIL and wait, encode as ASCII
        proc.communicate(input_str.encode('ascii'))

        # load polar, skip header
        data = np.loadtxt(fname, skiprows=12)
        alpha = data[:,0]
        Cd    = data[:,2]
        xtr_top = data[:,5]
        xtr_bot = data[:,6]

        # store for combined drag‐vs‐alpha plot
        sweep_results.append((xtr, alpha, Cd))
        
        # record Cd at α=0°
        idx0 = np.argmin(np.abs(alpha - 0.0))
        cd0_list.append(Cd[idx0])

        # plot forced‐transition vs α
        plt.figure()
        plt.plot(alpha, xtr_top, '^-', lw=1.5, label="Top Xtr")
        plt.plot(alpha, xtr_bot, 'v-', lw=1.5, label="Bot Xtr")
        plt.xlabel("α (deg)")
        plt.ylabel("Transition x/c")
        plt.title(f"Forced Transition @ x/c={xtr:.2f}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{plot_path}Xtr_vs_alpha_{xtr:.2f}.pdf", dpi=150)
        plt.close()

    # plot Cd at α=0° versus forced transition location
    plt.figure()
    plt.plot(xtr_list, cd0_list, "o-", lw=1.5)
    plt.xlabel("Forced Transition Location x/c")
    plt.ylabel("C_d @ α=0°")
    plt.title("Drag Coefficient vs. Forced Transition Point")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path + "Cd_vs_forced_xtr.pdf", dpi=150)
    plt.close()

    # plot drag coefficient vs angle of attack for all forced xtr in one figure
    plt.figure()
    for xtr, alpha, Cd in sweep_results:
        plt.plot(alpha, Cd, lw=1.5, label=f"xtr = {xtr:.2f}", marker = 'o')
    plt.xlabel("α (deg)")
    plt.ylabel("C_d")
    plt.title("Drag Coefficient vs Angle of Attack\nfor Forced Transition Locations")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path + "Cd_vs_alpha_forced_xtr.pdf", dpi=150)
    plt.close()


def run_xfoil_cf_distribution():
    """
    Compute and plot skin‑friction coefficient C_f vs. x/c on upper surface
    for α = 0° and α = 4°, Re=0.7e6, M=0, VPAR n=12.
    """
    alphas = [0.0, 4.0]
    cf_results = []

    for alpha in alphas:
        fname = f"{data_path}cf_alpha_{int(alpha)}.dat"
        try:
            os.remove(fname)
        except OSError:
            pass

        cmds = [
            "NACA 2009",
            "OPER",
            "VPAR",
            "N 12",
            "",
            "VISC 0.7e6",
            "MACH 0",
            "PACC",
            "",
            "",
            f"ALFA {alpha:.1f}",
            f"DUMP",    # write boundary‑layer data (x, C_f, ...)
            fname,
            "",
            "QUIT"
        ]

        # launch XFOIL
        proc = subprocess.Popen(
            [XFOIL_PATH],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        input_str = "\n".join(cmds) + "\n"
        # feed XFOIL and wait, encode as ASCII
        proc.communicate(input_str.encode('ascii'))

        # load C_f data; skip header line
        data = np.loadtxt(fname, skiprows=1)
        x_c   = data[:,1]
        cf_all = data[:,6]
        # assume first half of rows = upper surface
        mid = len(x_c)//2
        cf_results.append((alpha, x_c[:mid], cf_all[:mid]))

    # plot C_f vs x/c for both α
    plt.figure()
    for alpha, x_c, cf in cf_results:
        plt.plot(x_c, cf, lw=1, marker='.', label=f"$\\alpha$ = {alpha:.0f}°")
    plt.xlabel("x/c")
    plt.ylabel("C_f")
    plt.title("Skin‑Friction Coefficient vs Chord Position\nUpper Surface")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path + "Cf_vs_x_upper.pdf", dpi=150)
    plt.close()


def run_xfoil_cl_cd(cruise_cl=0.5):
    """
    Run XFOIL at fixed lift coefficient CL = cruise_cl,
    compute and print the Cl/Cd ratio.
    """
    fname = f"{data_path}polar_CL_{cruise_cl:.2f}.dat"
    # remove old file if present
    try:
        os.remove(fname)
    except OSError:
        pass

    # drive XFOIL: load airfoil, set VPAR, Re, Mach, then PACC + CL + stop
    cmds = [
        "NACA 2009",
        "OPER",
        "VPAR",
        "N 12",           # amplification factor
        "",
        "VISC 0.7e6",
        "MACH 0",
        "PACC", fname, "",        # start polar accumulation
        f"CL {cruise_cl:.2f}",    # fix lift coefficient
        "PACC", "", "",           # stop accumulation
        "QUIT"
    ]

    proc = subprocess.Popen(
        [XFOIL_PATH],
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    input_str = "\n".join(cmds) + "\n"
    # feed XFOIL and wait, encode as ASCII
    proc.communicate(input_str.encode('ascii'))

    # read the single‑point polar, skip header
    data = np.loadtxt(fname, skiprows=12)
    Cl = data[1]
    Cd = data[2]
    alpha = data[0]
    ratio = Cl / Cd
    print(f"At CL = {cruise_cl:.2f}:  alpha = {alpha:.4f}, Cd = {Cd:.5f},  Cl/Cd = {ratio:.1f}")
    return ratio

if __name__ == "__main__":
    # run_xfoil_2_4()
    # run_xfoil_fixed_transition_sweep()
    # run_xfoil_cf_distribution()
    run_xfoil_cl_cd(0.5)
    print("Done")