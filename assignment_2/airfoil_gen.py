import numpy as np
import matplotlib.pyplot as plt

def gen_from_csv():

    rescale_fac = 100

    airfoil_data_path = "assignment_2/naca_2009.csv"

    data = np.genfromtxt(airfoil_data_path,
                        dtype=float) 

    data /= rescale_fac

    x_data = np.append(data[:,0][::-1],data[:,2])
    y_data = np.append(data[:,1][::-1],data[:,3])

    coords = np.hstack((x_data[:,np.newaxis],y_data[:,np.newaxis],))

    with open('assignment_2/naca_2009.dat', 'w') as f:

        f.write('NACA_2009')
        f.write('\n')

        for i in coords:
            f.write(f'    {i[0]:.7f}   {i[1]:.7f}')
            f.write('\n')




    plt.plot(x_data,
            y_data, label = "online table ref")
    

def gen_from_xfoil_label():
    
    rescale_fac = 1

    airfoil_data_path = "assignment_2/b.dat"

    data = np.genfromtxt(airfoil_data_path,
                        dtype=float, skip_header=1) 

    x_data = data[:,0]
    y_data = data[:,1]

    plt.plot(x_data,
        y_data, label = "x_foil_generated")
    

gen_from_csv()
gen_from_xfoil_label()

plt.grid()
plt.legend()
plt.xlim(-0.2,1.2)
plt.ylim(-0.7,0.7)
plt.show()

