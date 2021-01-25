from dolfin import *
import numpy as np
import scipy.interpolate

# %% Parameters
#fiber diameters are taken from a normal distribution
fiber_mean = 8 #Mean fiber diameter (um)
fiber_std = 1 #Fiber diameter standard deviation (um)
lower_diam = 2 #Minimum fiber diameter (um)
n_axons = 500 #Number of fibers
num_fascicles = 1 #From create_mesh.py
k_fascicle = [0.08, 0.08, 0.5] #From solve_fenics.py
fasc_diam_x = [0.62] #From create_mesh.py
fasc_diam_y = [0.62] #From create_mesh.py
fasc_pos_x = [0] #From create_mesh.py
fasc_pos_y = [0] #From create_mesh.py
fasc_rotation = [0] #From create_mesh.py
nerve_pos_x = 0 #From create_mesh.py
nerve_pos_y = 0 #From create_mesh.py
contact_z = [-1, 0, 1] #From create_IF_mesh.py
waveform_type = 'Biphasic' #or Monophasic
phase_width = 0.1 #Width of each phase in waveform (ms)
cuff_length = 3 #From create_mesh.py
model_type = 'Cuff Electrode' #or 'Intrafascicular'
run_mode = 'voltage' #Which variable to use for passive fiber model. Only matters whether or not it is 'activating_function', which requires a 2nd order solution
#Options: 'axial_current_density', 'activating_function', 'total_current_density', 'radial_current_density', 'voltage'
fasc_rotation = np.radians(fasc_rotation).round(6)
    
if model_type == 'Cuff Electrode':
    model_length = max(20, cuff_length * 2 + 5)
else:
    unique_z = np.unique(contact_z).tolist()
    model_length = max(20, (unique_z[-1] - unique_z[0]) * 2 + 5)

# %% Read in solution    
mesh = Mesh()                                                                   #Create empty mesh
with XDMFFile('mesh.xdmf') as infile:                                           #Using main mesh file
    infile.read(mesh)                                                           #Read in the mesh
mvc = MeshValueCollection("size_t", mesh, 3)                                    #Create collection of cell markers
with XDMFFile('mesh_physical_region.xdmf') as infile:                           #Using mesh cell file    
    infile.read(mvc, "name_to_read")                                            #Read in the cell markers
markers = cpp.mesh.MeshFunctionSizet(mesh, mvc)                                 #Extract markers into variable

#Read solution
if run_mode == 'activating_function':
    V = FunctionSpace(mesh, 'P', 2);
    Q = FunctionSpace(mesh, 'DG', 0);
else:
    V = FunctionSpace(mesh, 'P', 1);
    Q = FunctionSpace(mesh, 'DG', 0);
u = Function(V)
f = HDF5File(mesh.mpi_comm(), 'solution.h5', "r")
f.read(u, "voltage")
f.close()

# %% Solve Fibers
fasc_markers = np.arange(2, num_fascicles+2)                                    #Determine the physical region markers assigned by gmsh for each fascicles                
if run_mode == 'total_current_density':
    E_x = project(u.dx(0), Q, solver_type="gmres", preconditioner_type="ilu")   #Calculate value of electric field in each direction
    E_y = project(u.dx(1), Q, solver_type="gmres", preconditioner_type="ilu")
    E_z = project(u.dx(2), Q, solver_type="gmres", preconditioner_type="ilu")
elif run_mode == 'axial_current_density':
    E_z = project(u.dx(2), Q, solver_type="gmres", preconditioner_type="ilu")   #Calculate electric field in z direction
elif run_mode == 'radial_current_density':
    E_x = project(u.dx(0), Q, solver_type="gmres", preconditioner_type="ilu")   #Calculate electric field in x and y directions
    E_y = project(u.dx(1), Q, solver_type="gmres", preconditioner_type="ilu")
if run_mode == 'activating_function':
    AF = project(u.dx(2).dx(2), Q, solver_type='gmres', preconditioner_type='ilu') #Calculate activating function
np.random.seed(0)                                                               #So results are the same every time
x_nodes_all = np.zeros(n_axons*num_fascicles)                                   #Create empty array for x positions of fibers
y_nodes_all = np.zeros(n_axons*num_fascicles)                                   #y positions            
z_nodes_all = np.zeros(n_axons*num_fascicles)                                   #z positions
diam_nodes_all = np.zeros(n_axons*num_fascicles)                                #diameters of each fiber
fascicle_nodes_all = np.empty(n_axons*num_fascicles, dtype='<U256')             #which fascicle each fiber is in
threshold_nodes_all = np.zeros(n_axons*num_fascicles)                           #Threshold current required to activate each fiber
max_val_nodes_all = np.zeros(n_axons*num_fascicles)                             #Max value of "run_mode" variable along each fiber at 1 mA
node_max = np.zeros((n_axons, num_fascicles))                                   #Max value of "run_mode" variable along each fiber at 1 mA
threshold = np.zeros((n_axons, num_fascicles))                                  #Threshold current required to activate each fiber
for k in range(len(fasc_markers)):                                              #Loop through each fascicle
    fascicle_elements = [val for val in markers.array() if val == fasc_markers[k]]  #Compute number of elements in the fascicle
    x_coord = np.zeros(len(fascicle_elements))                                  #Initialize vector for x coordinates of each element
    y_coord = np.zeros(len(fascicle_elements))                                  #Initialize vector for y coordinates
    z_coord = np.zeros(len(fascicle_elements))                                  #Initialize vector for z coordinates
    values = np.zeros(len(fascicle_elements))                                   #Initialize vector for voltage
    i = 0;                                                                      #Initialize counter
    for cell in cells(mesh):                                                    #Loop through all cells in the mesh
        if markers[cell.index()] == fasc_markers[k]:                            #If it is in the fascicle
            x_coord[i] = cell.midpoint().x()                                    #Save its x-coordinate    
            y_coord[i] = cell.midpoint().y()                                    #Save its y-coordinate
            z_coord[i] = cell.midpoint().z()                                    #Save its z-coordinate        
            if run_mode == 'axial_current_density':
                values[i] = E_z(cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()) * k_fascicle[2] * 1000 #Calculate axial current density
            elif run_mode == 'activating_function':
                values[i] = AF(cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()) #Calculate activating function
            elif run_mode == 'total_current_density':
                x_temp = E_x(cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()) * k_fascicle[0] * 1000
                y_temp = E_y(cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()) * k_fascicle[1] * 1000
                z_temp = E_z(cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()) * k_fascicle[2] * 1000
                values[i] = np.sqrt(x_temp**2 + y_temp**2 + z_temp**2)          #Calculate total current density
            elif run_mode == 'radial_current_density':
                x_temp = E_x(cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()) * k_fascicle[0] * 1000
                y_temp = E_y(cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()) * k_fascicle[1] * 1000
                values[i] = np.sqrt(x_temp**2 + y_temp**2)                      #Calculate radial current density
            elif run_mode == 'voltage':
                values[i] = (u(cell.midpoint().x(), cell.midpoint().y(),        #Calculate voltage
                             cell.midpoint().z()))
            i = i + 1;                                                          #Update counter

    coords = np.column_stack((x_coord, y_coord, z_coord))                       #Combine x,y,z coordinates of fascicle elements
    nn = scipy.interpolate.LinearNDInterpolator(coords, values)                 #Create interpolant for values in the fascicle
    x_nodes = np.zeros(n_axons)                                                 #Initialize x-coordinates of fibers
    y_nodes = np.zeros(n_axons)                                                 #Initialize y-coordinates of fibers
    z_nodes = np.zeros(n_axons)                                                 #Initialize z-coordinates of fibers
    node_diams = np.zeros(n_axons)                                              #Initialize diameters of fibers
    i = 0                                                                       #Set counter to 0
    while i < n_axons:                                                          #Loop until all fibers are created    
        while True:
            x_temp = np.random.uniform(-fasc_diam_x[k]/2, fasc_diam_x[k]/2, 1)  #Randomize x-position
            y_temp = np.random.uniform(-fasc_diam_y[k]/2, fasc_diam_y[k]/2, 1)  #Randomize y-position
            if (x_temp**2)/((fasc_diam_x[k]/2)**2)+(y_temp**2)/((fasc_diam_y[k]/2)**2) < 1: #Check that it is inside the fascicle
                break                                                           #If succesful, break the loop
        x = x_temp*np.cos(fasc_rotation[k]) - y_temp*np.sin(fasc_rotation[k])   #Rotate the position of the fiber along with the fascicle's rotation
        y = x_temp*np.sin(fasc_rotation[k]) + y_temp*np.cos(fasc_rotation[k])
        x = x + fasc_pos_x[k] + nerve_pos_x                                     #Translate position of fiber along with fascicle
        y = y + fasc_pos_y[k] + nerve_pos_y
        z = np.random.uniform(-0.5, 0.5, 1)                                     #Randomize z-coordinate (for the first node in the fiber)
        while True:
            diameter = np.random.normal(fiber_mean, fiber_std, 1)               #Randomize diameter of fiber
            if diameter > lower_diam:
                break
        z_step = (92.7652 * diameter + 108.9688) / 1000                         #Computer distance between nodes based on diameter
        z_temp = np.arange(z, model_length, z_step)                             #Create list of all nodes z-positions
        x_temp = np.ones(len(z_temp)) * x                                       #Create list of all nodes x-positions
        y_temp = np.ones(len(z_temp)) * y                                       #Create list of all nodes y-positions
        interp_coords = np.column_stack((x_temp, y_temp, z_temp))               #Combine x,y,z coordinates of nodes
        node_vals = nn(interp_coords)                                           #Interpolate value of "run_mode" variable to all nodes
        node_max_temp = np.nanmax(np.abs(node_vals))                            #Calculate max value of "run_mode" variable along the fiber    
        if not np.isnan(node_max_temp):
            x_nodes[i] = x                                                      #Save the x position
            y_nodes[i] = y                                                      #Save the y position
            z_nodes[i] = z                                                      #Save the z position
            node_diams[i] = diameter                                            #Save the diameter
            if run_mode == 'axial_current_density':                             #Calculate threshold variables based on the waveform and fiber diameter
                if waveform_type == 'Biphasic':
                    rheobase = 26.4426 * np.exp(-0.1584 * diameter) + 9.8285
                    chronaxie = 0.3695 * np.exp(-0.4630 * diameter) + 0.1206
                else:
                    rheobase = 38.3077 * np.exp(-0.2576 * diameter) + 13.1113
                    chronaxie = 0.4504 * np.exp(-1.2904 * diameter) + 0.1236
            elif run_mode == 'activating_function':
                if waveform_type == 'Biphasic':
                    rheobase = 0.1078 * np.exp(-0.1566 * diameter) + 0.0235
                    chronaxie = 0.3004 * np.exp(-0.3997 * diameter) + 0.1218
                else:
                    rheobase = 0.1338 * np.exp(-0.2043 * diameter) + 0.0325
                    chronaxie = 0.3666 * np.exp(-1.1497 * diameter) + 0.1240
            elif run_mode == 'voltage':
                if waveform_type == 'Biphasic':
                    rheobase = 0.0351 * np.exp(-0.4736 * diameter) + 0.0190
                    chronaxie = 0.2694 * np.exp(-0.2310 * diameter) + 0.0892
                else:
                    rheobase = 0.0553 * np.exp(-0.3958 * diameter) + 0.0173
                    chronaxie = 0.8375 * np.exp(-2.2482 * diameter) + 0.1323
            elif run_mode == 'total_current_density':
                if waveform_type == 'Biphasic':
                    rheobase = 26.1664 * np.exp(-0.1639 * diameter) + 10.3065
                    chronaxie = 0.3661 * np.exp(-0.4592 * diameter) + 0.1204
                else:
                    rheobase = 38.6359 * np.exp(-0.2644 * diameter) + 13.3910
                    chronaxie = 0.4518 * np.exp(-1.3096 * diameter) + 0.1239
            elif run_mode == 'radial_current_density':
                if waveform_type == 'Biphasic':
                    rheobase = 5.0102 * np.exp(-0.1616 * diameter) + 0.7789
                    chronaxie = 0.3559 * np.exp(-0.7352 * diameter) + 0.1518
                else:
                    rheobase = 5.1350 * np.exp(-0.1939 * diameter) + 1.4459
                    chronaxie = 0.7359 * np.exp(-1.5834 * diameter) + 0.1228
            threshold[i][k] = rheobase * (1 + chronaxie / phase_width)          #Calculate threshold for the fiber    
            node_max[i][k] = node_max_temp                                      #Save max value of "run_mode" variable along the fiber
            x_nodes_all[i+k*n_axons] = x                                        #Save fiber x-position
            y_nodes_all[i+k*n_axons] = y                                        #Save fiber y-position
            z_nodes_all[i+k*n_axons] = z                                        #Save fiber starting z-position
            diam_nodes_all[i+k*n_axons] = diameter                              #Save fiber diameter
            fascicle_nodes_all[i+k*n_axons] = 'fascicle'+str(k+1)               #Save the identity of the fascicle that the fiber is in
            threshold_nodes_all[i+k*n_axons] = threshold[i][k]                  #Save threshold of the fiber
            max_val_nodes_all[i+k*n_axons] = node_max_temp                      #Save max value of "run_mode" variable along the fiber
            i = i + 1                                                           #Update counter
current_thresholds = np.divide(threshold, node_max)                             #Calculate current necessary to activate each fiber
current_thresholds[~np.isfinite(current_thresholds)] = 0                        #Remove any thresholds that are "infinite"
column_max = np.round(np.nanmax(current_thresholds, axis=0) * 1.2, 6)           #Calculate current range for each fascicle. Will go from 0 mA to 1.2 times the max threshold within each fascicle
column_max = np.sort(column_max)                                                #Sort the max thresholds
prev_max = 0
column_current = np.array([])
for i in range(len(column_max)):                                                #Loop through each fascicle and break the range for that fascicle into 100 intensities
    column_current = np.append(column_current, np.arange(prev_max, column_max[i], column_max[i] / 100)) #For each subsequent fascicle, start at the end of the previous fascicle's range
    prev_max = column_max[i]
current = column_current
recruitment = np.zeros((len(current), num_fascicles))                           #Initialize recruitment vector
for j in range(num_fascicles):                                                  #Loop through each fascicle
    for i in range(len(current)):                                               #Loop through the current intensities
        recruitment[i][j] = np.sum((node_max[:,j] * current[i]) > threshold[:,j])*100/n_axons #Compute percent of fibers activated at each intensity    

array = np.array(np.column_stack((current, recruitment)))                       #Combine current and recruitment into matrix
np.savetxt('recruitment.txt', array, fmt='%1.4e')
fiber_array = np.array(np.column_stack((np.round(x_nodes_all, 8), np.round(y_nodes_all, 8), #Combine all fiber values into an array
            np.round(z_nodes_all, 8), fascicle_nodes_all, np.round(diam_nodes_all, 8), 
            np.round(threshold_nodes_all, 8), np.round(max_val_nodes_all, 8))))
np.savetxt('fibers.txt', fiber_array, fmt='%s')