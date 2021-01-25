import pygmsh
import numpy as np
from scipy import optimize
import os
import meshio

# %% Define Functions
def overlap(params, a, b, x0, y0, theta):
    sum1 = (np.cos(theta)*(params[0]-x0)+np.sin(theta)*(params[1]-y0))**2/(a**2)
    sum2 = (np.sin(theta)*(params[0]-x0)-np.cos(theta)*(params[1]-y0))**2/(b**2)
    return sum1+sum2

def no_overlap(params, a, b, x0, y0, theta):
    sum1 = (np.cos(theta)*(params[0]-x0)+np.sin(theta)*(params[1]-y0))**2/(a**2)
    sum2 = (np.sin(theta)*(params[0]-x0)-np.cos(theta)*(params[1]-y0))**2/(b**2)
    return -(sum1+sum2)

def get_id(element):
    return element.id[0:element.id.index('[')]+'[1]';

# %% Parameters
nerve_diameter = 0.8
fascicle_diameter_x = [0.62]
fascicle_diameter_y = [0.62]
fascicle_position_x = [0]
fascicle_position_y = [0]
fascicle_rotation = [0]
contact_diam = 0.01
contact_x = [0, 0, 0]
contact_y = [0, 0, 0]
contact_z = [-1, 0, 1]
tip_length = 0.01
        
contact_pos = [] #Create list of all contact coordinates
for i in range(len(contact_x)):
    contact_pos = contact_pos + [[contact_x[i],contact_y[i],contact_z[i]]]

contact_width = np.sqrt(np.pi * (contact_diam / 2) ** 2) #Adjust square width to have same area as circle with given diameter
fascicle_rotation = np.radians(fascicle_rotation) #Convert to radians
perineurium_diameter_x = np.multiply(fascicle_diameter_x, 1.06) #Calculate perineurium outer diameter (thickness = 3% fascicle diameter)
perineurium_diameter_y = np.multiply(fascicle_diameter_y, 1.06)
unique_z = np.unique(contact_z).tolist() #Unique list of z positions

model_length = max(20, (unique_z[-1] - unique_z[0]) * 2 + 5) #Total length of model
model_diameter = nerve_diameter * 4 #Total diameter of model

segs = [0.0] #List of segments in z-direction
for i in range(len(unique_z)):
    segs = segs + [model_length/2 + unique_z[i] - contact_width/2]
    segs = segs + [model_length/2 + unique_z[i] + contact_width/2]
segs = segs + [model_length]
segs.sort()

seg_starts = segs[0:-1] #Starting point of each segment
seg_lengths = np.diff(segs) #Length of each segment
segments = len(seg_starts) #Number of segments

# %% MESH PARAMS
lc_outer = 2*model_diameter/30 #Mesh resolution outside of nerve
lc_nerve = nerve_diameter/15 #Mesh resolution inside nerve
    
min_x = min(np.subtract(perineurium_diameter_x,fascicle_diameter_x))
min_y = min(np.subtract(perineurium_diameter_y,fascicle_diameter_y))
abs_min = min(min_x,min_y,contact_width,tip_length,lc_nerve) #Absolute minimum size of each mesh element

#Create geometry object
geom = pygmsh.opencascade.Geometry( 
  characteristic_length_min=abs_min,
  characteristic_length_max=lc_outer,
)
geom.add_raw_code('Geometry.OldNewReg=0;') #To keep track of geometry tag numbers

# %% GEOMETRY
#Contacts
num_fascicles = len(fascicle_diameter_x)
unique_x = np.unique(contact_x).tolist()
upper_y = np.add(contact_y, tip_length).tolist()
max_y = nerve_diameter / 2 * 1.25 #Make electrodes stick up out of the nerve
unique_y = np.unique(contact_y+upper_y).tolist()+[max_y]
fasc_cont = []
fasc_cont_coords = []
peri_cont = []
peri_cont_coords = []
epi_cont = []
epi_cont_coords = []
outer_cont = []
outer_cont_coords = []
fasc_cont_index = []
for k in range(len(unique_x)):
    x1 = unique_x[k] - contact_width / 2
    x2 = unique_x[k] + contact_width / 2
    for j in range(len(unique_y) - 1):
        y1 = unique_y[j]
        y2 = unique_y[j+1]
        non_results_peri = []
        for i in range(num_fascicles):
            a_fasc=fascicle_diameter_x[i]/2
            b_fasc=fascicle_diameter_y[i]/2
            x0_fasc=fascicle_position_x[i]
            y0_fasc=fascicle_position_y[i]
            theta_fasc=fascicle_rotation[i]
            result = optimize.minimize(overlap,[np.mean([x1,x2]),np.mean([y1,y2])],args=(a_fasc,b_fasc,x0_fasc,y0_fasc,theta_fasc),bounds=[[x1,x2],[y1,y2]])
            non_result_fasc = optimize.minimize(no_overlap,[np.mean([x1,x2]),np.mean([y1,y2])],args=(a_fasc,b_fasc,x0_fasc,y0_fasc,theta_fasc),bounds=[[x1,x2],[y1,y2]])
            if result.fun < 1: #If electrode overlaps with fascicle
                if fascicle_diameter_x[i] < fascicle_diameter_y[i]:
                    temp_fasc = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_y[i]/2,fascicle_diameter_x[i]/2)
                    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(temp_fasc.id)+'};}')
                else:
                    temp_fasc = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_x[i]/2,fascicle_diameter_y[i]/2)
                geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(temp_fasc.id)+'};}')
                temp_cont = geom.add_rectangle([x1, y1, 0], x2 - x1, y2 - y1)
                fasc_cont = fasc_cont + [geom.boolean_intersection([temp_cont, temp_fasc])]
                fasc_cont_coords = fasc_cont_coords + [[unique_x[k], unique_y[j]]]
                fasc_cont_index = fasc_cont_index + [i]
            a_peri=perineurium_diameter_x[i]/2
            b_peri=perineurium_diameter_y[i]/2
            result = optimize.minimize(overlap,[np.mean([x1,x2]),np.mean([y1,y2])],args=(a_peri,b_peri,x0_fasc,y0_fasc,theta_fasc),bounds=[[x1,x2],[y1,y2]])
            non_result_peri = optimize.minimize(no_overlap,[np.mean([x1,x2]),np.mean([y1,y2])],args=(a_peri,b_peri,x0_fasc,y0_fasc,theta_fasc),bounds=[[x1,x2],[y1,y2]])
            non_results_peri = non_results_peri + [non_result_peri.fun]
            if result.fun < 1 and -non_result_fasc.fun > 1: #If electrode overlaps with perineurium and isn't solely contained in fascicle's ellipse
                if fascicle_diameter_x[i] < fascicle_diameter_y[i]:
                    peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_y[i]/2,perineurium_diameter_x[i]/2)
                    temp_fasc = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_y[i]/2,fascicle_diameter_x[i]/2)
                    temp_peri = geom.boolean_difference([peri_disk],[temp_fasc],delete_first=True,delete_other=True)
                    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(temp_peri.id)+'};}')
                else:
                    peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_x[i]/2,perineurium_diameter_y[i]/2)
                    temp_fasc = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_x[i]/2,fascicle_diameter_y[i]/2)
                    temp_peri = geom.boolean_difference([peri_disk],[temp_fasc],delete_first=True,delete_other=True)
                geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(temp_peri.id)+'};}')
                temp_cont = geom.add_rectangle([x1, y1, 0], x2 - x1, y2 - y1)
                peri_cont = peri_cont + [geom.boolean_intersection([temp_cont, temp_peri])]
                peri_cont_coords = peri_cont_coords + [[unique_x[k], unique_y[j]]]
        a=nerve_diameter/2
        b=nerve_diameter/2
        x0=0
        y0=0
        theta=0
        result = optimize.minimize(overlap,[np.mean([x1,x2]),np.mean([y1,y2])],args=(a,b,x0,y0,theta),bounds=[[x1,x2],[y1,y2]])
        non_result_nerve = optimize.minimize(no_overlap,[np.mean([x1,x2]),np.mean([y1,y2])],args=(a,b,x0,y0,theta),bounds=[[x1,x2],[y1,y2]])
        if result.fun < 1 and -np.max(non_results_peri) > 1: #Overlaps with nerve and not solely contained in perineurium's ellipse
            epineurium_disk = geom.add_disk([0,0,0],nerve_diameter/2,nerve_diameter/2)
            temp_fasc = [0 for i in range(len(fascicle_diameter_x))]
            temp_peri = [0 for i in range(len(fascicle_diameter_x))]
            for i in range(len(fascicle_diameter_x)):
                if fascicle_diameter_x[i] < fascicle_diameter_y[i]:
                    peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_y[i]/2,perineurium_diameter_x[i]/2)
                    temp_fasc[i] = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_y[i]/2,fascicle_diameter_x[i]/2)
                    temp_peri[i] = geom.boolean_difference([peri_disk],[temp_fasc[i]],delete_first=True,delete_other=False)
                    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(temp_peri[i].id)+'};}')
                    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(temp_fasc[i].id)+'};}')
                else:
                    peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_x[i]/2,perineurium_diameter_y[i]/2)
                    temp_fasc[i] = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_x[i]/2,fascicle_diameter_y[i]/2)
                    temp_peri[i] = geom.boolean_difference([peri_disk],[temp_fasc[i]],delete_first=True,delete_other=False)
                geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(temp_peri[i].id)+'};}')
                geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(temp_fasc[i].id)+'};}')
            temp_epi = geom.boolean_difference([epineurium_disk],temp_peri+temp_fasc,delete_first=True,delete_other=True)
            temp_cont = geom.add_rectangle([x1, y1, 0], x2 - x1, y2 - y1)
            epi_cont = epi_cont + [geom.boolean_intersection([temp_cont, temp_epi])]
            epi_cont_coords = epi_cont_coords + [[unique_x[k], unique_y[j]]]
        if -non_result_nerve.fun > 1:
            epineurium_disk = geom.add_disk([0,0,0],nerve_diameter/2,nerve_diameter/2)
            temp_cont = geom.add_rectangle([x1, y1, 0], x2 - x1, y2 - y1)
            outer_cont = outer_cont + [geom.boolean_difference([temp_cont],[epineurium_disk],delete_first=True,delete_other=True)]
            outer_cont_coords = outer_cont_coords + [[unique_x[k], unique_y[j]]]
    
#Nerve
epineurium_disk = geom.add_disk([0,0,0],nerve_diameter/2,nerve_diameter/2)    
fasc = [0 for i in range(len(fascicle_diameter_x))]    
peri = [0 for i in range(len(fascicle_diameter_x))]
for i in range(len(fascicle_diameter_x)):
    if fascicle_diameter_x[i] < fascicle_diameter_y[i]:
        peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_y[i]/2,perineurium_diameter_x[i]/2)
        fasc_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_y[i]/2,fascicle_diameter_x[i]/2)
        geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(peri_disk.id)+'};}')
        geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(fasc_disk.id)+'};}')
    else:
        peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_x[i]/2,perineurium_diameter_y[i]/2)
        fasc_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_x[i]/2,fascicle_diameter_y[i]/2)
    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(peri_disk.id)+'};}')
    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(fasc_disk.id)+'};}')
    peri[i] = geom.boolean_difference([peri_disk],[fasc_disk]+peri_cont,delete_other=False)
    if len(fasc_cont) > 0:
        fasc[i] = geom.boolean_difference([fasc_disk],fasc_cont,delete_other=False)
    else:
        fasc[i] = fasc_disk
epi = geom.boolean_difference([epineurium_disk],peri+fasc+epi_cont+fasc_cont+peri_cont,delete_first=False,delete_other=False)
    
#Saline
sal_disk = geom.add_disk([0,0,0],model_diameter/2,model_diameter/2)
sal = geom.boolean_difference([sal_disk],[epineurium_disk],delete_other=True)
sal = geom.boolean_difference([sal],outer_cont,delete_other=False)
    
geom.boolean_fragments([sal],fasc_cont+peri_cont+epi_cont+outer_cont+[epi]) #Run boolean fragments to clean up any tiny issues

# %% Extrude
char_len = lc_outer
sal_extrude = [0 for i in range(segments)]
fasc_cont_extrude = [[0 for i in range(len(fasc_cont))] for j in range(segments)]
num_raw = len(peri_cont+epi_cont+outer_cont+fasc+peri+[epi])
for i in range(segments):
    if i == 0 or i == segments - 1:
        char_len = lc_outer
    elif i%2 == 1: #odd numbers are electrode rows
        char_len = min(contact_width,lc_nerve)
    elif i%2 == 0: #even numbers are between electrode rows
        char_len = lc_nerve               
    if seg_lengths[i] < char_len:
        char_len = seg_lengths[i]        
    temp_length = seg_lengths[i]
    temp_layers = int(round(seg_lengths[i]/char_len))
    
    for j in range(num_fascicles):
        new_id = 1000 + j + num_raw * i
        if i == 0:
            surface = fasc[j].id
        else:
            surface = 'fasc_surfaces_'+str(i-1)+'_'+str(j)+'[]'
        geom.add_raw_code(
'''start_v = newv;
ex'''+str(new_id)+'''[] = Extrude {0,0,'''+str(temp_length)+'''} {Surface{'''+surface+'''}; Layers{'''+str(temp_layers)+'''}; };
end_v = newv;
N = #ex'''+str(new_id)+'''[];
fasc_volumes_'''+str(i)+'''_'''+str(j)+''' = {};
fasc_surfaces_'''+str(i)+'''_'''+str(j)+''' = {};
For i In {start_v : end_v - 1}
    For j In {0 : N - 1}
        If(ex'''+str(new_id)+'''[j] == i)
            fasc_surfaces_'''+str(i)+'''_'''+str(j)+'''[#fasc_surfaces_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j-1];
            fasc_volumes_'''+str(i)+'''_'''+str(j)+'''[#fasc_volumes_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j];
        EndIf
    EndFor
EndFor
'''
        )
    
    for j in range(num_fascicles):
        new_id = 1000 + num_fascicles + j + num_raw * i
        if i == 0:
            surface = peri[j].id
        else:
            surface = 'peri_surfaces_'+str(i-1)+'_'+str(j)+'[]'
        geom.add_raw_code(
'''start_v = newv;
ex'''+str(new_id)+'''[] = Extrude {0,0,'''+str(temp_length)+'''} {Surface{'''+surface+'''}; Layers{'''+str(temp_layers)+'''}; };
end_v = newv;
N = #ex'''+str(new_id)+'''[];
peri_volumes_'''+str(i)+'''_'''+str(j)+''' = {};
peri_surfaces_'''+str(i)+'''_'''+str(j)+''' = {};
For i In {start_v : end_v - 1}
    For j In {0 : N - 1}
        If(ex'''+str(new_id)+'''[j] == i)
            peri_surfaces_'''+str(i)+'''_'''+str(j)+'''[#peri_surfaces_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j-1];
            peri_volumes_'''+str(i)+'''_'''+str(j)+'''[#peri_volumes_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j];
        EndIf
    EndFor
EndFor
'''
        )
        
    new_id = 1000 + num_fascicles * 2 + num_raw * i
    if i == 0:
        surface = epi.id
    else:
        surface = 'epi_surfaces_'+str(i-1)+'[]'
    geom.add_raw_code(
'''start_v = newv;
ex'''+str(new_id)+'''[] = Extrude {0,0,'''+str(temp_length)+'''} {Surface{'''+surface+'''}; Layers{'''+str(temp_layers)+'''}; };
end_v = newv;
epi_volumes_'''+str(i)+''' = {};
epi_surfaces_'''+str(i)+''' = {};
N = #ex'''+str(new_id)+'''[];
For i In {start_v : end_v - 1}
    For j In {0 : N - 1}
        If(ex'''+str(new_id)+'''[j] == i)
            epi_surfaces_'''+str(i)+'''[#epi_surfaces_'''+str(i)+'''[]] = ex'''+str(new_id)+'''[j-1];
            epi_volumes_'''+str(i)+'''[#epi_volumes_'''+str(i)+'''[]] = ex'''+str(new_id)+'''[j];
        EndIf
    EndFor
EndFor
'''
        )
    
    #fasc_cont cannot get split into multiple regions
    for j in range(len(fasc_cont)):
        if i == 0:
            fasc_cont_extrude[i][j] = geom.extrude(
                fasc_cont[j],
                translation_axis=[0,0,temp_length],
                num_layers=temp_layers
            )
        else:
            fasc_cont_extrude[i][j] = geom.extrude(
                fasc_cont_extrude[i-1][j][0],
                translation_axis=[0,0,temp_length],
                num_layers=temp_layers
            )
        id_temp = get_id(fasc_cont_extrude[i][j][1])
        id_temp = id_temp[0:id_temp.index('[')]
    
    #peri_cont and epi_cont can get split
    for j in range(len(peri_cont)):
        new_id = 1000 + num_fascicles * 2 + 1 + j + num_raw * i
        if i == 0:
            surface = peri_cont[j].id
        else:
            surface = 'peri_cont_surfaces_'+str(i-1)+'_'+str(j)+'[]'
        geom.add_raw_code(
'''start_v = newv;
ex'''+str(new_id)+'''[] = Extrude {0,0,'''+str(temp_length)+'''} {Surface{'''+surface+'''}; Layers{'''+str(temp_layers)+'''}; };
end_v = newv;
N = #ex'''+str(new_id)+'''[];
peri_cont_volumes_'''+str(i)+'''_'''+str(j)+''' = {};
peri_cont_surfaces_'''+str(i)+'''_'''+str(j)+''' = {};
For i In {start_v : end_v - 1}
    For j In {0 : N - 1}
        If(ex'''+str(new_id)+'''[j] == i)
            peri_cont_surfaces_'''+str(i)+'''_'''+str(j)+'''[#peri_cont_surfaces_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j-1];
            peri_cont_volumes_'''+str(i)+'''_'''+str(j)+'''[#peri_cont_volumes_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j];
        EndIf
    EndFor
EndFor
'''
        )
        
    for j in range(len(epi_cont)):
        new_id = 1000 + num_fascicles * 2 + 1 + len(peri_cont) + j + num_raw * i
        if i == 0:
            surface = epi_cont[j].id
        else:
            surface = 'epi_cont_surfaces_'+str(i-1)+'_'+str(j)+'[]'
        geom.add_raw_code(
'''start_v = newv;
ex'''+str(new_id)+'''[] = Extrude {0,0,'''+str(temp_length)+'''} {Surface{'''+surface+'''}; Layers{'''+str(temp_layers)+'''}; };
end_v = newv;
N = #ex'''+str(new_id)+'''[];
epi_cont_volumes_'''+str(i)+'''_'''+str(j)+''' = {};
epi_cont_surfaces_'''+str(i)+'''_'''+str(j)+''' = {};
For i In {start_v : end_v - 1}
    For j In {0 : N - 1}
        If(ex'''+str(new_id)+'''[j] == i)
            epi_cont_surfaces_'''+str(i)+'''_'''+str(j)+'''[#epi_cont_surfaces_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j-1];
            epi_cont_volumes_'''+str(i)+'''_'''+str(j)+'''[#epi_cont_volumes_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j];
        EndIf
    EndFor
EndFor
'''
        )
    
    #outer_cont can get split    
    for j in range(len(outer_cont)):
        new_id = 1000 + num_fascicles * 2 + 1 + len(peri_cont) + len(epi_cont) + j + num_raw * i
        if i == 0:
            surface = outer_cont[j].id
        else:
            surface = 'outer_cont_surfaces_'+str(i-1)+'_'+str(j)+'[]'
        geom.add_raw_code(
'''start_v = newv;
ex'''+str(new_id)+'''[] = Extrude {0,0,'''+str(temp_length)+'''} {Surface{'''+surface+'''}; Layers{'''+str(temp_layers)+'''}; };
end_v = newv;
N = #ex'''+str(new_id)+'''[];
outer_cont_volumes_'''+str(i)+'''_'''+str(j)+''' = {};
outer_cont_surfaces_'''+str(i)+'''_'''+str(j)+''' = {};
For i In {start_v : end_v - 1}
    For j In {0 : N - 1}
        If(ex'''+str(new_id)+'''[j] == i)
            outer_cont_surfaces_'''+str(i)+'''_'''+str(j)+'''[#outer_cont_surfaces_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j-1];
            outer_cont_volumes_'''+str(i)+'''_'''+str(j)+'''[#outer_cont_volumes_'''+str(i)+'''_'''+str(j)+'''[]] = ex'''+str(new_id)+'''[j];
        EndIf
    EndFor
EndFor
'''
        )
    
    if i == 0:
        sal_extrude[i] = geom.extrude(
            sal,
            translation_axis=[0,0,temp_length],
            num_layers=temp_layers
        )
    else:
        sal_extrude[i] = geom.extrude(
            sal_extrude[i-1][0],
            translation_axis=[0,0,temp_length],
            num_layers=temp_layers
        )
    
# %% Tag Each Region
# contacts
num_contacts = len(contact_x)
all_contact_coords = [[x,y,z] for(x,y,z) in zip(contact_x,contact_y,contact_z)]
num_fasc = len(fasc_cont)
num_peri = len(peri_cont)
num_epi = len(epi_cont)
num_outer = len(outer_cont)
    
contact_matrix = []
raw_contact_matrix = []
for i in range(num_contacts):
    z_index = 1 + 2 * unique_z.index(contact_z[i])
    xy = [contact_x[i], contact_y[i]]
    if xy in fasc_cont_coords:
        index = [j for j,x in enumerate(fasc_cont_coords) if x == xy] #can be multiple (1 per fascicle)
        contact_matrix = contact_matrix + [fasc_cont_extrude[z_index][j][1] for j in index]
    if xy in peri_cont_coords:
        index = peri_cont_coords.index(xy)
        raw_contact_matrix = raw_contact_matrix + ['peri_cont_volumes_'+str(z_index)+'_'+str(index)+'[]']
    if xy in epi_cont_coords:
        index = epi_cont_coords.index(xy)
        raw_contact_matrix = raw_contact_matrix + ['epi_cont_volumes_'+str(z_index)+'_'+str(index)+'[]']
    if xy in outer_cont_coords:
        index = outer_cont_coords.index(xy)
        raw_contact_matrix = raw_contact_matrix + ['outer_cont_volumes_'+str(z_index)+'_'+str(index)+'[]']
contact_list = [get_id(i) for i in contact_matrix]
if contact_list and not raw_contact_matrix:
    geom.add_raw_code('Physical Volume("contacts") = {'+','.join(contact_list)+'};')
elif not contact_list and raw_contact_matrix:
    geom.add_raw_code('Physical Volume("contacts") = {'+','.join(raw_contact_matrix)+'};')
else:
    geom.add_raw_code('Physical Volume("contacts") = {'+','.join(contact_list)+','+','.join(raw_contact_matrix)+'};')
    

cuff_coords = []
for i in range(num_contacts):
    x = contact_x[i]
    y = [j for j in unique_y[0:-1] if j > contact_y[i]]
    z = contact_z[i]
    cuff_coords = cuff_coords + [[x,j,z] for j in y]
cuff_matrix = []
raw_cuff_matrix = []
for i in range(len(cuff_coords)):
    z_index = 1 + 2 * unique_z.index(cuff_coords[i][2])
    xy = [cuff_coords[i][0], cuff_coords[i][1]]
    if xy in fasc_cont_coords:
        index = [j for j,x in enumerate(fasc_cont_coords) if x == xy] #can be multiple (1 per fascicle)
        cuff_matrix = cuff_matrix + [fasc_cont_extrude[z_index][j][1] for j in index]
    if xy in peri_cont_coords:
        index = [j for j,x in enumerate(peri_cont_coords) if x == xy] #can be multiple (1 per fascicle)
        raw_cuff_matrix = raw_cuff_matrix + ['peri_cont_volumes_'+str(z_index)+'_'+str(j)+'[]' for j in index]
    if xy in epi_cont_coords:
        index = epi_cont_coords.index(xy)
        raw_cuff_matrix = raw_cuff_matrix + ['epi_cont_volumes_'+str(z_index)+'_'+str(index)+'[]']
    if xy in outer_cont_coords:
        index = outer_cont_coords.index(xy)
        raw_cuff_matrix = raw_cuff_matrix + ['outer_cont_volumes_'+str(z_index)+'_'+str(index)+'[]']
cuff_list = [get_id(i) for i in cuff_matrix]

# fascicle
for i in range(num_fascicles):
    fascicle_matrix = []
    raw_fascicle_matrix = []
    for j in range(segments):
        raw_fascicle_matrix = raw_fascicle_matrix + ['fasc_volumes_'+str(j)+'_'+str(i)+'[]']
    z_indices = [1 + 2 * j for j in range(len(unique_z))]
    not_z = [j for j in range(segments) if not j in z_indices]
    for ii in range(len(not_z)):
        for jj in range(num_fasc):
            if fasc_cont_index[jj] == i:
                fascicle_matrix = fascicle_matrix + [fasc_cont_extrude[not_z[ii]][jj][1]]
    for ii in range(len(z_indices)):
        for jj in range(num_fasc):
            coords = fasc_cont_coords[jj]+[unique_z[ii]]
            if coords not in all_contact_coords and coords not in cuff_coords and fasc_cont_index[jj] == i:
                fascicle_matrix = fascicle_matrix + [fasc_cont_extrude[z_indices[ii]][jj][1]]
    fasc_list = [get_id(i) for i in fascicle_matrix]
    if fasc_list:
        geom.add_raw_code('Physical Volume("fascicle'+str(i+1)+'") = {'+','.join(fasc_list)+','+','.join(raw_fascicle_matrix)+'};')
    else:
        geom.add_raw_code('Physical Volume("fascicle'+str(i+1)+'") = {'+','.join(raw_fascicle_matrix)+'};')

# perineurium
perineurium_matrix = []
for i in range(num_fascicles):
    for j in range(segments):
        perineurium_matrix = perineurium_matrix + ['peri_volumes_'+str(j)+'_'+str(i)+'[]']
for i in range(len(not_z)):
    for j in range(num_peri):
        perineurium_matrix = perineurium_matrix + ['peri_cont_volumes_'+str(not_z[i])+'_'+str(j)+'[]']
for ii in range(len(z_indices)):
    for jj in range(num_peri):
        coords = peri_cont_coords[jj]+[unique_z[ii]]
        if coords not in all_contact_coords and coords not in cuff_coords:
            perineurium_matrix = perineurium_matrix + ['peri_cont_volumes_'+str(z_indices[ii])+'_'+str(jj)+'[]']
geom.add_raw_code('Physical Volume("perineurium") = {'+','.join(perineurium_matrix)+'};')

# epineurium
epineurium_matrix = []
for i in range(segments):
    epineurium_matrix = epineurium_matrix + ['epi_volumes_'+str(i)+'[]']
for i in range(len(not_z)):
    for j in range(num_epi):
        epineurium_matrix = epineurium_matrix + ['epi_cont_volumes_'+str(not_z[i])+'_'+str(j)+'[]']
for ii in range(len(z_indices)):
    for jj in range(num_epi):
        coords = epi_cont_coords[jj]+[unique_z[ii]]
        if coords not in all_contact_coords and coords not in cuff_coords:
            epineurium_matrix = epineurium_matrix + ['epi_cont_volumes_'+str(z_indices[ii])+'_'+str(jj)+'[]']
geom.add_raw_code('Physical Volume("epineurium") = {'+','.join(epineurium_matrix)+'};')
    
if cuff_list and not raw_cuff_matrix:
    geom.add_raw_code('Physical Volume("cuff") = {'+','.join(cuff_list)+'};')
elif not cuff_list and raw_cuff_matrix:
    geom.add_raw_code('Physical Volume("cuff") = {'+','.join(raw_cuff_matrix)+'};')
else:
    geom.add_raw_code('Physical Volume("cuff") = {'+','.join(cuff_list)+','+','.join(raw_cuff_matrix)+'};')

# saline
saline_matrix = []
raw_saline_matrix = []
for i in range(segments):
    saline_matrix = saline_matrix + [sal_extrude[i][1]]
for i in range(len(not_z)):
    for j in range(num_outer):
        raw_saline_matrix = raw_saline_matrix + ['outer_cont_volumes_'+str(not_z[i])+'_'+str(j)+'[]']
for ii in range(len(z_indices)):
    for jj in range(num_outer):
        coords = outer_cont_coords[jj]+[unique_z[ii]]
        if coords not in all_contact_coords and coords not in cuff_coords:
            raw_saline_matrix = raw_saline_matrix + ['outer_cont_volumes_'+str(z_indices[ii])+'_'+str(jj)+'[]']
saline_list = [get_id(i) for i in saline_matrix]
geom.add_raw_code('Physical Volume("saline") = {'+','.join(saline_list)+','+','.join(raw_saline_matrix)+'};')

# %% Adaptive Mesh
geom.add_raw_code(
'''
model_length = '''+str(model_length)+''';
blip = '''+str(lc_nerve)+'''; // ensure field is outside of geometry where necessary
lc_outer = '''+str(lc_outer)+''';
peri_lc = '''+str(abs_min)+''';
lc_nerve = '''+str(lc_nerve)+''';
nerve_diameter = '''+str(nerve_diameter)+''';
'''
)

radius2 = nerve_diameter/2

# Create adaptive mesh fields in gmsh
geom.add_raw_code(
'''
Field[1] = Cylinder;
Field[1].XCenter = '''+str(0)+''';
Field[1].YCenter = '''+str(0)+''';
Field[1].ZCenter = model_length/2;
Field[1].XAxis = 0;
Field[1].YAxis = 0;
Field[1].ZAxis = model_length + blip;
Field[1].Radius = '''+str(radius2)+''';
Field[1].VIn = lc_nerve;
Field[1].VOut = lc_outer;
'''
)

field_string = '1'
for i in range(len(fascicle_diameter_x)):
    f_num = 2 + i
    f_rad = max(fascicle_diameter_x[i],fascicle_diameter_y[i])*1.2/2
    peri_thickness = np.mean([perineurium_diameter_x[i]-fascicle_diameter_x[i],perineurium_diameter_y[i]-fascicle_diameter_y[i]])
    f_lc_min = min(3*peri_thickness,lc_nerve)
    field_string = field_string + ',' + str(f_num)
    geom.add_raw_code('''
Field['''+str(f_num)+'''] = Cylinder;
Field['''+str(f_num)+'''].XCenter = '''+str(fascicle_position_x[i])+''';
Field['''+str(f_num)+'''].YCenter = '''+str(fascicle_position_y[i])+''';
Field['''+str(f_num)+'''].ZCenter = model_length/2;
Field['''+str(f_num)+'''].XAxis = 0;
Field['''+str(f_num)+'''].YAxis = 0;
Field['''+str(f_num)+'''].ZAxis = model_length + blip;
Field['''+str(f_num)+'''].Radius = '''+str(f_rad)+''';
Field['''+str(f_num)+'''].VIn = '''+str(f_lc_min)+''';
Field['''+str(f_num)+'''].VOut = lc_outer;
''')

final_num = 1 + len(fascicle_diameter_x) + 1
geom.add_raw_code('''
Field['''+str(final_num)+'''] = Min;
Field['''+str(final_num)+'''].FieldsList = {'''+field_string+'''};
Background Field = '''+str(final_num)+''';
'''
)

# %% Save
f = open('mesh.geo','w')
f.write(geom.get_code());
f.close();

os.system('sudo gmsh/bin/gmsh -3 mesh.geo') #Create mesh using (sudo path/to/gmsh -dimension filename)
msh = meshio.read('mesh.msh') #Read mesh into meshio
meshio.write('mesh.xdmf', meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]})) #Convert to xdmf
meshio.write('mesh_physical_region.xdmf', meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]},
                 cell_data={"tetra": {"name_to_read": msh.cell_data["tetra"]["gmsh:physical"]}}))