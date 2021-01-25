from dolfin import *
import numpy as np

def outer(x, on_boundary):                                                      #Define the surface
    global model_length
    return (on_boundary and not near(x[2], 0, DOLFIN_EPS)                       #Return all points on outside that are not on either end
            and not near(x[2], model_length, DOLFIN_EPS))                       #The two ends have Neumann BC by default

electrode_currents = [-1, 1]
k_saline = [2, 2, 2] #Conductivity in x,y,z direction
k_fascicle = [0.08, 0.08, 0.5]
k_epineurium = [0.008, 0.008, 0.008]
k_perineurium = [0.00336, 0.00336, 0.00336]
k_cuff = [2e-10, 2e-10, 2e-10]
k_electrodes = 9e6
num_fascicles = 1 #Based on create_mesh.py
intra_fasc = 0 #Set to 1 if using create_IF_mesh.py
contact_x = [0, 0, 0] #Only used if intra_fasc=1
contact_y = [0, 0, 0] #Only used if intra_fasc=1
contact_z = [-1, 0, 1] #Only used if intra_fasc=1
run_mode = 'voltage' #Which variable to use for passive fiber model. Only matters whether or not it is 'activating_function', which requires a 2nd order solution
#Options: 'axial_current_density', 'activating_function', 'total_current_density', 'radial_current_density', 'voltage'

# %% Load Mesh   
mesh = Mesh()                                                                   #Create empty mesh
with XDMFFile('mesh.xdmf') as infile:                                           #Using main mesh file
    infile.read(mesh)                                                           #Read in the mesh
mvc = MeshValueCollection("size_t", mesh, 3)                                    #Create collection of cell markers
with XDMFFile('mesh_physical_region.xdmf') as infile:                           #Using mesh cell file    
    infile.read(mvc, "name_to_read")                                            #Read in the cell markers
markers = cpp.mesh.MeshFunctionSizet(mesh, mvc)                                 #Extract markers into variable
if intra_fasc:                                                                  #For intra_fasc=1, electrode faces will be identified manually
    faces = MeshFunction('size_t', mesh, 2);
    faces.set_all(0)
else:
    mvc = MeshValueCollection("size_t", mesh, 2)                                #Create collection of surface markers
    with XDMFFile('mesh_facet_region.xdmf') as infile:                          #Using mesh facet file
        infile.read(mvc, "name_to_read")                                        #Read in the facet markers
    faces = cpp.mesh.MeshFunctionSizet(mesh, mvc)                               #Extract markers into variable

if run_mode == 'activating_function':                                           #Activating function requires a 2nd order solution (much slower)
    V = FunctionSpace(mesh, 'P', 2);
    W = VectorFunctionSpace(mesh, 'P', 1)
else:
    V = FunctionSpace(mesh, 'P', 1);
    W = VectorFunctionSpace(mesh, 'P', 1)
    
# %% Ground Outer Surface of Model
global model_length
model_length = 0
for v in vertices(mesh):                                                        #Find length of model
    z = v.point().z()
    if z > model_length:
        model_length = z
                
bc_outer = DirichletBC(V, 0, outer)                                             #Create the boundary condition

# %% Define Equations
marker_array = markers.array()                                                  #List of all facet labels in mesh
marker_array = np.unique(marker_array)                                          #Unique values in the list
marker_array = [val for val in marker_array if val != 0]                        #Remove the value 0

#Basic equation is "dot(grad(u), grad(v)) * dx == f * v * dx"
dx = Measure("dx", domain=mesh, subdomain_data=markers)                         #Create "dx" measure
u = TrialFunction(V)                                                            #Define trial function
v = TestFunction(V)                                                             #Define test function    
f = Constant(0)                                                                 #Set f to zero as the net charge in the model is 0    

u_electrodes = (as_vector([k_electrodes[0]*grad(u)[0],                          #Adjust grad(u) for electrodes
                k_electrodes[1]*grad(u)[1], k_electrodes[2]*grad(u)[2]]))
a_electrodes = dot(u_electrodes, grad(v)) * dx(1)                               #Define dot product for all electrodes

u_fascicle = (as_vector([k_fascicle[0]*grad(u)[0],                              #Adjust grad(u) for the fascicle based on material properties
              k_fascicle[1]*grad(u)[1], k_fascicle[2]*grad(u)[2]]))
a_fascicle = dot(u_fascicle, grad(v))*dx(2)                                     #Define dot product for fascicle
for i in range(num_fascicles-1):
    a_fascicle = a_fascicle + dot(u_fascicle, grad(v))*dx(3+i)
    
u_perineurium = (as_vector([k_perineurium[0]*grad(u)[0],                        #Adjust grad(u) for perineurium
                 k_perineurium[1]*grad(u)[1], k_perineurium[2]*grad(u)[2]]))
a_perineurium = dot(u_perineurium, grad(v)) * dx(num_fascicles+2)               #Define dot product for perineurium

u_epineurium = (as_vector([k_epineurium[0]*grad(u)[0],                          #Adjust grad(u) for epineurium
                k_epineurium[1]*grad(u)[1], k_epineurium[2]*grad(u)[2]]))
a_epineurium = dot(u_epineurium, grad(v)) * dx(num_fascicles+3)                 #Define dot product for epineurium

u_cuff = (as_vector([k_cuff[0]*grad(u)[0],                                      #Adjust grad(u) for cuff
          k_cuff[1]*grad(u)[1], k_cuff[2]*grad(u)[2]]))
a_cuff = dot(u_cuff, grad(v)) * dx(num_fascicles+4)                             #Define dot product for cuff

u_saline = (as_vector([k_saline[0]*grad(u)[0],                                  #Adjust grad(u) for saline
            k_saline[1]*grad(u)[1], k_saline[2]*grad(u)[2]]))
a_saline = dot(u_saline, grad(v)) * dx(num_fascicles+5)                         #Define dot product for saline

a = (a_fascicle + a_perineurium + a_epineurium                                  #Combine into left side of the equation    
     + a_electrodes + a_cuff + a_saline)
L = f*v*dx                                                                      #Define right side of the equation
u = Function(V)                                                                 #Define u to be solved for
ds = Measure('ds', domain=mesh, subdomain_data=faces)                           #Create "ds" measure (outer boundary)
n = FacetNormal(mesh)                                                           #Define normal vector to each facet (for computing flux)

# %% Calibrate Voltage on each Electrode
if intra_fasc:                                                                  #For IF mesh, must loop through faces and manually find electrodes and label them
    faces.set_all(0)
    dim = mesh.topology().dim();                                                #3-dimensional
    mesh.init(dim - 1, dim);                                                    #Connect 2d to 3d (faces of cells)
    contact_z = [i + model_length / 2 for i in contact_z]                       #All contact z-positions    
    for f in facets(mesh):                                                      #Loop through all facets in the mesh
        if 1 in markers.array()[f.entities(dim)]: #1=electrodes                 #If it's a facet on an electrode
            vert = f.entities(0)                                                #Read in all vertices of the facet
            x = np.mean([Vertex(mesh, i).point().x() for i in vert])            #Calculate x-position of facet using the vertices
            y = 0                                                               #y-position doesn't matter
            z = np.mean([Vertex(mesh, i).point().z() for i in vert])            #Calculate z-position of facet using the vertices
            distance = model_length                                             #Initialize distance variable
            for i in range(len(contact_x)):                                     #Loop through each electrode and find out which one is closest to the current facet
                contact_point = Point(contact_x[i], 0, contact_z[i])
                current_point = Point(x,y,z)
                if contact_point.distance(current_point) < distance:
                    distance = contact_point.distance(current_point)
                    closest_index = i
            faces[f.index()] = closest_index + 1                                #Label the current facet with whichever electrode is nearest
                
face_array = faces.array()                                                      #List of all facet labels in mesh
face_array = np.unique(face_array)                                              #Unique values in the list
face_array = [val for val in face_array if val != 0]                            #Remove the value 0
all_BCs = [bc_outer]                                                            #List for all boundary conditions (only contains outer boundary so far)

indices = [i for i,j in enumerate(electrode_currents) if j==-1]                 #Find all electrodes set to -1 mA
if indices:
    bcs = [bc_outer]                                                            #Initialize BC list
    for i in indices:                                                           #Loop through all anodic electrodes
        bc = DirichletBC(V, 1, faces, face_array[i])                            #Apply unit voltage
        bcs = bcs + [bc]                                                        #Update BC list
    (solve(a == L, u, bcs,                                                      #Solve system with voltage applied to the anode only
     solver_parameters={"linear_solver": "gmres", "preconditioner": "petsc_amg"}, 
     form_compiler_parameters={"optimize": True}))
    grad_u = (project(grad(u), W, solver_type='gmres',                          #Project voltage gradient into function space
      preconditioner_type='ilu'))
    flux_outer = assemble(k_saline[0] * dot(grad_u, n) * ds) / 1000            #Compute flux exiting model
    anode_voltage = 0.001 / flux_outer                                            #Determine new voltage for anodic electrodes
    for i in indices:                                                           #Loop through anodic electrodes again
        BC_temp = DirichletBC(V, anode_voltage, faces, face_array[i])             #Apply updated voltage to each electrode
        all_BCs = all_BCs + [BC_temp]                                           #Update BC list     
    

indices = [i for i,j in enumerate(electrode_currents) if j==1]                  #Find all electrodes set to +1 mA
if indices:
    bcs = [bc_outer]                                                            #Initialize BC list
    for i in indices:                                                           #Loop through all cathodic electrodes
        bc = DirichletBC(V, 1, faces, face_array[i])                            #Apply unit voltage
        bcs = bcs + [bc]                                                        #Update BC list    
    (solve(a == L, u, bcs,                                                      #Solve system with voltage applied to the single electrode
     solver_parameters={"linear_solver": "gmres", "preconditioner": "petsc_amg"},     
     form_compiler_parameters={"optimize": True}))
    grad_u = (project(grad(u), W, solver_type='gmres',                          #Project voltage gradient into function space
      preconditioner_type='ilu'))
    flux_outer = assemble(k_saline[0] * dot(grad_u, n) * ds) / 1000            #Compute flux exiting model
    cathode_voltage = -0.001 / flux_outer                                           #Determine new voltage for cathodic electrodes   
    for i in indices:                                                           #Loop through cathodic electrodes again
        BC_temp = DirichletBC(V, cathode_voltage, faces, face_array[i])             #Apply updated voltage to each electrode
        all_BCs = all_BCs + [BC_temp]                                           #Update BC list
  
                
# %% Solve Model
(solve(a == L, u, all_BCs,                                                      #Solve full model with BCs applied to all electrodes
 solver_parameters={"linear_solver": "gmres", "preconditioner": "petsc_amg"}, 
 form_compiler_parameters={"optimize": True}))
        
# %% Save Solution
f = HDF5File(mesh.mpi_comm(), "solution.h5", "w")
f.write(u, "voltage")
f.close()