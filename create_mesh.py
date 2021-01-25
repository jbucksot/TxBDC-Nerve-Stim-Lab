import pygmsh
import numpy as np
from math import radians
import os
import meshio

# %% Parameters
nerve_diameter = 0.8
nerve_position_x = 0
nerve_position_y = 0
fascicle_diameter_x = [0.62]
fascicle_diameter_y = [0.62]
fascicle_position_x = [0]
fascicle_position_y = [0]
fascicle_rotation = [0]
cuff_inner_diameter = 1
cuff_outer_diameter = 2
cuff_angle = 355
cuff_length = 3
contact_width = 0.01
contact_height = 0.01
contact_z_positions = [-0.5, 0.5]
contact_angles = [270]
contact_rotation = [0]

fascicle_rotation = np.radians(fascicle_rotation) #Convert to radians
old_cr = contact_rotation
contact_rotation = np.radians(contact_rotation) #Convert to radians
fascicle_position_x = np.add(fascicle_position_x, nerve_position_x) #Shift fascicles along with the nerve
fascicle_position_y = np.add(fascicle_position_y, nerve_position_y)
perineurium_diameter_x = np.multiply(fascicle_diameter_x, 1.06) #Calculate perineurium outer diameter (thickness = 3% fascicle diameter)
perineurium_diameter_y = np.multiply(fascicle_diameter_y, 1.06)
model_length = max(20, cuff_length * 3) #Total length of model
model_diameter = cuff_outer_diameter * 2 + 2 #Total radius of model

#Create list of segment positions (going in z-direction)
contact_z_positions = np.sort(contact_z_positions)
segs = [0.0]
segs = segs + [model_length/2 - cuff_length/2]
for i in range(len(contact_z_positions)):
    segs = segs + [model_length/2 + contact_z_positions[i] - contact_width/2]
    segs = segs + [model_length/2 + contact_z_positions[i] + contact_width/2]
segs = segs + [model_length/2 + cuff_length/2]
segs = segs + [model_length]
segs.sort()

seg_starts = segs[0:-1] #Where each segment starts
seg_lengths = np.diff(segs) #Length of each segment
segments = len(seg_starts) #Number of segments

# %% Define Mesh Parameters
lc_inner = round(cuff_inner_diameter/15, 8) #Mesh resolution inside cuff
lc_outer = round(model_diameter/15, 8) #Mesh resolution outside cuff
lc_nerve = round(nerve_diameter/15, 8) #Mesh resolution inside nerve 
    
min_x = min(np.subtract(perineurium_diameter_x,fascicle_diameter_x)) #Calculate perineurium thickness
min_y = min(np.subtract(perineurium_diameter_y,fascicle_diameter_y))
min_x = round(min_x, 8)
min_y = round(min_y, 8)
abs_min = min(min_x,min_y,contact_width,contact_height,lc_nerve) #Set absolute minimum size of mesh elements

#Create geometry object
geom = pygmsh.opencascade.Geometry(
  characteristic_length_min=abs_min,
  characteristic_length_max=lc_outer,
)

# %% Create Geometry
#Nerve
epineurium_disk = geom.add_disk([nerve_position_x,nerve_position_y,0],nerve_diameter/2,nerve_diameter/2)    
fasc = [0 for i in range(len(fascicle_diameter_x))]    
peri = [0 for i in range(len(fascicle_diameter_x))]
for i in range(len(fascicle_diameter_x)):
    if fascicle_diameter_x[i] < fascicle_diameter_y[i]:
        peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_y[i]/2,perineurium_diameter_x[i]/2)
        fasc[i] = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_y[i]/2,fascicle_diameter_x[i]/2)
        peri[i] = geom.boolean_difference([peri_disk],[fasc[i]],delete_other=False)
        geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(peri[i].id)+'};}')
        geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},Pi/2} {Surface{'+str(fasc[i].id)+'};}')
    else:
        peri_disk = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],perineurium_diameter_x[i]/2,perineurium_diameter_y[i]/2)
        fasc[i] = geom.add_disk([fascicle_position_x[i],fascicle_position_y[i],0],fascicle_diameter_x[i]/2,fascicle_diameter_y[i]/2)
        peri[i] = geom.boolean_difference([peri_disk],[fasc[i]],delete_other=False)
    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(peri[i].id)+'};}')
    geom.add_raw_code('Rotate {{0,0,1},{'+str(fascicle_position_x[i])+','+str(fascicle_position_y[i])+',0},'+str(fascicle_rotation[i])+'} {Surface{'+str(fasc[i].id)+'};}')
epi = geom.boolean_difference([epineurium_disk],peri+fasc,delete_first=False,delete_other=False)

#Contacts
cont = [0 for i in range(len(contact_angles))]
for i in range(len(contact_angles)):
    alpha1 = 360-contact_angles[i]
    beta1 = 90-alpha1/2
    px = (cuff_inner_diameter/2)*np.cos(radians(beta1))
    py = (cuff_inner_diameter/2)*np.sin(radians(beta1)) 
    if contact_angles[i] < 360:
        cont_outer_disk_whole = geom.add_disk([0,0,0],cuff_inner_diameter/2,cuff_inner_diameter/2)
        cont_angle_tool = geom.add_polygon([[0,0,0],[px,py,0],[model_diameter/2,py,0],[model_diameter/2,model_diameter/2,0],[-model_diameter/2,model_diameter/2,0],[-model_diameter/2,py,0],[-px,py,0]])
        cont_outer_disk = geom.boolean_difference([cont_outer_disk_whole],[cont_angle_tool],delete_first=True,delete_other=True)
    else:
        cont_outer_disk = geom.add_disk([0,0,0],cuff_inner_diameter/2,cuff_inner_diameter/2)
    if round(nerve_diameter, 5) == round((cuff_inner_diameter - 2 * contact_height), 5):
        cont[i] = geom.boolean_difference([cont_outer_disk],[epineurium_disk],delete_first=True,delete_other=False)
    else:
        cont_inner_disk = geom.add_disk([0,0,0],cuff_inner_diameter/2-contact_height,cuff_inner_diameter/2-contact_height)
        cont[i] = geom.boolean_difference([cont_outer_disk],[cont_inner_disk],delete_first=True,delete_other=True)
    geom.add_raw_code('Rotate {{0,0,1},{0,0,0},'+str(contact_rotation[i])+'} {Surface{'+str(cont[i].id)+'};}')

#Cuff
alpha2 = 360-cuff_angle
beta2 = 90-alpha2/2
pxc = (cuff_outer_diameter/2)*np.cos(radians(beta2))
pyc = (cuff_outer_diameter/2)*np.sin(radians(beta2))  
cuff_outer_disk = geom.add_disk([0,0,0],cuff_outer_diameter/2,cuff_outer_diameter/2)
cont_outer_disk_whole = geom.add_disk([0,0,0],cuff_inner_diameter/2,cuff_inner_diameter/2)
if cuff_angle < 360:
    cuff_angle_tool = geom.add_polygon([[0,0,0],[pxc,pyc,0],[model_diameter/2,pyc,0],[model_diameter/2,model_diameter/2,0],[-model_diameter/2,model_diameter/2,0],[-model_diameter/2,pyc,0],[-pxc,pyc,0]])
    cuff_cut = geom.boolean_difference([cuff_outer_disk],[cuff_angle_tool],delete_other=True)
    cuff = geom.boolean_difference([cuff_cut],[cont_outer_disk_whole],delete_other=True)
else:
    if round(nerve_diameter, 5) == round((cuff_inner_diameter - 2 * contact_height), 5):
        cuff = geom.boolean_difference([cuff_outer_disk],[cont_outer_disk_whole],delete_first=False,delete_other=True)
    else:
        cuff = geom.boolean_difference([cuff_outer_disk],[cont_outer_disk_whole],delete_first=False,delete_other=False)
    
#Saline
if cuff_angle < 360:
    if round(nerve_diameter, 5) == round((cuff_inner_diameter - 2 * contact_height), 5):
        old_cr = np.mod(np.array(old_cr) + 360, 360)
        sorted_contact_angles = np.array([x for _,x in sorted(zip(old_cr, contact_angles))])
        sorted_contact_pos = np.array(sorted(old_cr))
        contact_start_angles = sorted_contact_pos - sorted_contact_angles / 2
        contact_end_angles = sorted_contact_pos + sorted_contact_angles / 2
        section_start_angles = contact_end_angles
        section_end_angles = contact_start_angles[[x for x in range(1, len(contact_start_angles))] + [0]]
        section_lengths = np.mod((section_end_angles - section_start_angles) + 360, 360)
        cuff_start_angle = -cuff_angle / 2
        cuff_end_angle = cuff_angle / 2
        section_start_angles[section_start_angles > 180] = section_start_angles[section_start_angles > 180] - 360
        section_end_angles[section_end_angles > 180] = section_end_angles[section_end_angles > 180] - 360
        standalone_sections = [] #If there are any sections of saline not connected to the rest of the saline, they need to be created independently
        for i in range(len(section_start_angles)):
            if section_start_angles[i] > cuff_start_angle and section_end_angles[i] < cuff_end_angle and section_lengths[i] > 0 and section_start_angles[i] < section_end_angles[i]:
                standalone_sections.append(i)
        standalone_sal = [0 for i in range(len(standalone_sections))]
        for i in range(len(standalone_sections)):
            alpha1 = 360 - section_lengths[standalone_sections[i]]
            beta1 = 90 - alpha1 / 2
            px = (cuff_inner_diameter/2)*np.cos(radians(beta1))
            py = (cuff_inner_diameter/2)*np.sin(radians(beta1))
            sal_outer_disk_whole = geom.add_disk([0,0,0],cuff_inner_diameter/2,cuff_inner_diameter/2)
            sal_angle_tool = geom.add_polygon([[0,0,0],[px,py,0],[model_diameter/2,py,0],[model_diameter/2,model_diameter/2,0],[-model_diameter/2,model_diameter/2,0],[-model_diameter/2,py,0],[-px,py,0]])
            sal_outer_disk = geom.boolean_difference([sal_outer_disk_whole],[sal_angle_tool],delete_first=True,delete_other=True)
            standalone_sal[i] = geom.boolean_difference([sal_outer_disk],[epineurium_disk],delete_other=False)
            temp_rotation = section_start_angles[standalone_sections[i]] + section_lengths[standalone_sections[i]] / 2
            temp_rotation = radians(temp_rotation)
            geom.add_raw_code('Rotate {{0,0,1},{0,0,0},'+str(temp_rotation)+'} {Surface{'+str(standalone_sal[i].id)+'};}')
        sal_disk = geom.add_disk([0,0,0],model_diameter/2,model_diameter/2)
        sal = geom.boolean_difference([sal_disk],[epineurium_disk],delete_other=True)
        sal = [geom.boolean_difference([sal],[cuff]+cont+standalone_sal,delete_other=False)]
        sal = sal+standalone_sal
    else:
        sal_disk = geom.add_disk([0,0,0],model_diameter/2,model_diameter/2)
        sal = geom.boolean_difference([sal_disk],[epineurium_disk],delete_other=True)
        sal = [geom.boolean_difference([sal],[cuff]+cont,delete_other=False)]
else:
    if round(nerve_diameter, 5) == round((cuff_inner_diameter - 2 * contact_height), 5):
        old_cr = np.mod(np.array(old_cr) + 360, 360)
        sorted_contact_angles = np.array([x for _,x in sorted(zip(old_cr, contact_angles))])
        sorted_contact_pos = np.array(sorted(old_cr))
        contact_start_angles = sorted_contact_pos - sorted_contact_angles / 2
        contact_end_angles = sorted_contact_pos + sorted_contact_angles / 2
        section_start_angles = contact_end_angles
        section_end_angles = contact_start_angles[[x for x in range(1, len(contact_start_angles))] + [0]]
        section_lengths = np.mod((section_end_angles - section_start_angles) + 360, 360)
        section_start_angles[section_start_angles > 180] = section_start_angles[section_start_angles > 180] - 360
        section_end_angles[section_end_angles > 180] = section_end_angles[section_end_angles > 180] - 360
        standalone_sections = []
        for i in range(len(section_start_angles)):
            if section_lengths[i] > 0:
                standalone_sections.append(i)
        standalone_sal = [0 for i in range(len(standalone_sections))]
        for i in range(len(standalone_sections)):
            alpha1 = 360 - section_lengths[standalone_sections[i]]
            beta1 = 90 - alpha1 / 2
            px = (cuff_inner_diameter/2)*np.cos(radians(beta1))
            py = (cuff_inner_diameter/2)*np.sin(radians(beta1))
            sal_outer_disk_whole = geom.add_disk([0,0,0],cuff_inner_diameter/2,cuff_inner_diameter/2)
            sal_angle_tool = geom.add_polygon([[0,0,0],[px,py,0],[model_diameter/2,py,0],[model_diameter/2,model_diameter/2,0],[-model_diameter/2,model_diameter/2,0],[-model_diameter/2,py,0],[-px,py,0]])
            sal_outer_disk = geom.boolean_difference([sal_outer_disk_whole],[sal_angle_tool],delete_first=True,delete_other=True)
            if i == len(standalone_sections) - 1:
                standalone_sal[i] = geom.boolean_difference([sal_outer_disk],[epineurium_disk],delete_other=True)
            else:
                standalone_sal[i] = geom.boolean_difference([sal_outer_disk],[epineurium_disk],delete_other=False)                    
            temp_rotation = section_start_angles[standalone_sections[i]] + section_lengths[standalone_sections[i]] / 2
            temp_rotation = radians(temp_rotation)
            geom.add_raw_code('Rotate {{0,0,1},{0,0,0},'+str(temp_rotation)+'} {Surface{'+str(standalone_sal[i].id)+'};}')
        sal_outer_disk = geom.add_disk([0,0,0],model_diameter/2,model_diameter/2)
        sal_outer = geom.boolean_difference([sal_outer_disk],[cuff_outer_disk],delete_first=True,delete_other=True)
        sal = [sal_outer] + standalone_sal
    else:
        sal_outer_disk = geom.add_disk([0,0,0],model_diameter/2,model_diameter/2)
        sal_outer = geom.boolean_difference([sal_outer_disk],[cuff_outer_disk],delete_first=True,delete_other=True)
        sal_inner_disk = cont_outer_disk_whole
        sal_inner = geom.boolean_difference([sal_inner_disk],[epineurium_disk],delete_other=True)
        sal_inner = geom.boolean_difference([sal_inner],cont,delete_other=False)
        sal = [sal_outer, sal_inner]
    
geom.boolean_fragments(sal,[cuff]+cont+[epi]) #Run boolean fragments to clean up any tiny issues
geom.boolean_fragments(cont,[epi,cuff])
geom.boolean_fragments(cont, sal)

# %% Extrude
s = cont+fasc+peri+[epi,cuff]+sal #List of all the regions

v=[[0 for i in range(len(s))] for j in range(segments)] #Grand list of all extruded pieces
for j in range(len(s)):   
    char_len = lc_outer * 2
    v[0][j] = geom.extrude(
        s[j],
        translation_axis=[0,0,seg_lengths[0]],
        rotation_axis=None,
        point_on_axis=None,
        angle=None,
        num_layers=int(round(seg_lengths[0]/char_len)),
        recombine=False    
    )    
for i in range(segments-1):
    if i == segments-2:
        char_len = lc_outer * 2             
    elif (i-1)%2 == 0:
        char_len = min(contact_width,lc_inner)
    else:
        char_len = lc_inner * 2
    if seg_lengths[i+1] < char_len:
        char_len = seg_lengths[i+1]
    for j in range(len(s)): # i = segment, j = surface   
        v[i+1][j] = geom.extrude(
            v[i][j][0],
            translation_axis=[0,0,seg_lengths[i+1]],
            rotation_axis=None,
            point_on_axis=None,
            angle=None,
            num_layers=int(round(seg_lengths[i+1]/char_len)),
            recombine=False    
        )
    
# %% Tag Each Region
num_fascicles = len(fascicle_diameter_x)
num_contacts = len(contact_angles)
# contacts
contact_matrix = [0 for i in range(len(contact_z_positions)*len(contact_angles))]
for i in range(len(contact_z_positions)):
    for j in range(num_contacts):
        contact_matrix[i*num_contacts+j] = v[2+2*i][j][1]
geom.add_physical(contact_matrix, label='contacts')

# fascicle
for i in range(num_fascicles):
    fascicle_matrix = [0 for i in range(segments)]
    for j in range(segments):
        fascicle_matrix[j]=v[j][num_contacts+i][1]
    geom.add_physical(fascicle_matrix, label='fascicle'+str(i+1))

# perineurium
perineurium_matrix = [0 for i in range(segments*num_fascicles)]
for i in range(num_fascicles):
    for j in range(segments):
        perineurium_matrix[segments*i+j]=v[j][num_contacts+num_fascicles+i][1]
geom.add_physical(perineurium_matrix, label='perineurium')

# epineurium
epineurium_matrix = [0 for i in range(segments)]
for i in range(segments):
    epineurium_matrix[i] = v[i][num_contacts+num_fascicles*2][1]
geom.add_physical(epineurium_matrix, label='epineurium')
    
# cuff
cuff_matrix = [0 for i in range(segments - 2)]
for i in range(segments - 2):
    cuff_matrix[i] = v[i+1][num_contacts+num_fascicles*2+1][1]
geom.add_physical(cuff_matrix, label='cuff')

# saline
saline_matrix = [0 for i in range(len(sal) * segments)]
for i in range(len(sal)):
    for j in range(segments):
        saline_matrix[segments*i+j]=v[j][num_contacts+num_fascicles*2+2+i][1]
for i in range(num_contacts):
    saline_matrix = saline_matrix + [v[0][i][1]]                                #First electrode segment
    saline_matrix = saline_matrix + [v[1][i][1]]                                #Second electrode segment
for i in range(len(contact_z_positions)-1):
    for j in range(num_contacts):
        saline_matrix = saline_matrix + [v[2+2*i+1][j][1]]                      #Segments in-between electrodes
for i in range(num_contacts):
    saline_matrix = saline_matrix + [v[-2][i][1]]                               #Second-to-last electrode segment
    saline_matrix = saline_matrix + [v[-1][i][1]]                               #Last electrode segment
saline_matrix = saline_matrix + [v[0][num_fascicles*2+1+num_contacts][1]]       #First cuff segment
saline_matrix = saline_matrix + [v[-1][num_fascicles*2+1+num_contacts][1]]      #Last cuff segment
geom.add_physical(saline_matrix, label='saline')

#Contact Surfaces
for i in range(len(contact_z_positions)):
    for j in range(num_contacts):
        ind = len(s) * 2 * (i + 1) + 1 + j                                      #Skip first two segments, then take 2nd extrusion from desired segment
        num = i * num_contacts + j + 1
        if contact_angles[j] == 360:
            geom.add_raw_code(
            '''Physical Surface("contact'''+str(num)+'''") = {ex'''+str(ind)+'''[0], ex'''+str(ind)+'''[2], ex'''+str(ind)+'''[3], ex'''+str(ind-len(s))+'''[0]};'''
            )
        elif contact_angles[j] > 180: 
            geom.add_raw_code(
            '''Physical Surface("contact'''+str(num)+'''") = {ex'''+str(ind)+'''[0], ex'''+str(ind)+'''[2], ex'''+str(ind)+'''[3], ex'''+str(ind)+'''[4], ex'''+str(ind)+'''[5], ex'''+str(ind)+'''[6], ex'''+str(ind)+'''[7], ex'''+str(ind-len(s))+'''[0]};'''
            )
        else:
            geom.add_raw_code(
            '''Physical Surface("contact'''+str(num)+'''") = {ex'''+str(ind)+'''[0], ex'''+str(ind)+'''[2], ex'''+str(ind)+'''[3], ex'''+str(ind)+'''[4], ex'''+str(ind)+'''[5], ex'''+str(ind-len(s))+'''[0]};'''
            )

# %% Adaptive Mesh
# Feed necessary parameters to gmsh
geom.add_raw_code(
'''
cuff_outer_diameter = '''+str(cuff_outer_diameter)+''';
cuff_inner_diameter = '''+str(cuff_inner_diameter)+''';
model_length = '''+str(model_length)+''';
blip = '''+str(lc_inner)+'''; // ensure field is outside of geometry where necessary
lc_inner = '''+str(lc_inner)+''';
lc_outer = '''+str(lc_outer)+''';
peri_lc = '''+str(abs_min)+''';
lc_nerve = '''+str(lc_nerve)+''';
nerve_diameter = '''+str(nerve_diameter)+''';
'''
)

if nerve_diameter == (cuff_inner_diameter - 2 * contact_height):
    radius2 = cuff_inner_diameter/2
else:
    radius2 = nerve_diameter/2

# Create adaptive mesh fields in gmsh
geom.add_raw_code(
'''
Field[1] = Cylinder;
Field[1].XCenter = 0;
Field[1].YCenter = 0;
Field[1].ZCenter = model_length/2;
Field[1].XAxis = 0;
Field[1].YAxis = 0;
Field[1].ZAxis = model_length + blip;
Field[1].Radius = cuff_outer_diameter/2 + blip;
Field[1].VIn = lc_inner;
Field[1].VOut = lc_outer;

Field[2] = Cylinder;
Field[2].XCenter = '''+str(nerve_position_x)+''';
Field[2].YCenter = '''+str(nerve_position_y)+''';
Field[2].ZCenter = model_length/2;
Field[2].XAxis = 0;
Field[2].YAxis = 0;
Field[2].ZAxis = model_length + blip;
Field[2].Radius = '''+str(radius2)+''';
Field[2].VIn = lc_nerve;
Field[2].VOut = lc_outer;
'''
)

field_string = '1,2'
for i in range(len(fascicle_diameter_x)):
    f_num = 3 + 2 * i
    peri_thickness = np.min([perineurium_diameter_x[i]-fascicle_diameter_x[i],perineurium_diameter_y[i]-fascicle_diameter_y[i]])
    peri_thickness = round(peri_thickness, 8)
    f_lc_min = min(2*peri_thickness,lc_nerve)
    field_string = field_string + ',' + str(f_num+1)
    peri_face_list = [peri[i].id[0:peri[i].id.index('[')]] + [x[num_contacts+num_fascicles+i][0].id for x in v]
    geom.add_raw_code('''
Field['''+str(f_num)+'''] = Distance;
Field['''+str(f_num)+'''].FacesList = {'''+','.join(peri_face_list)+'''};

Field['''+str(f_num+1)+'''] = Threshold;
Field['''+str(f_num+1)+'''].IField = '''+str(f_num)+''';
Field['''+str(f_num+1)+'''].LcMin = '''+str(f_lc_min)+''';
Field['''+str(f_num+1)+'''].LcMax = '''+str(lc_nerve)+''';
Field['''+str(f_num+1)+'''].DistMin = 0;
Field['''+str(f_num+1)+'''].DistMax = '''+str(peri_thickness * 2)+''';
Field['''+str(f_num+1)+'''].StopAtDistMax = 1;

''')
    
for i in range(len(contact_angles)):
    f_num = 3 + 2 * len(fascicle_diameter_x) + 2 * i
    f_lc_min = min(lc_nerve, contact_height)
    field_string = field_string + ',' + str(f_num + 1)
    cont_face_list = [cont[i].id[0:cont[i].id.index('[')]] + [x[i][0].id for x in v]
    geom.add_raw_code('''
Field['''+str(f_num)+'''] = Distance;
Field['''+str(f_num)+'''].FacesList = {'''+','.join(cont_face_list)+'''};

Field['''+str(f_num+1)+'''] = Threshold;
Field['''+str(f_num+1)+'''].IField = '''+str(f_num)+''';
Field['''+str(f_num+1)+'''].LcMin = '''+str(f_lc_min)+''';
Field['''+str(f_num+1)+'''].LcMax = '''+str(lc_nerve)+''';
Field['''+str(f_num+1)+'''].DistMin = 0;
Field['''+str(f_num+1)+'''].DistMax = '''+str(contact_height)+''';
Field['''+str(f_num+1)+'''].StopAtDistMax = 1;

''')
    
if round(nerve_diameter, 5) == round((cuff_inner_diameter - 2 * contact_height), 5) and len(sal) > 1:
    num_isolated = len(sal) - 1
    for i in range(num_isolated):
        f_num = 3 + 2 * len(fascicle_diameter_x) + 2 * len(contact_angles) + 2 * i
        f_lc_min = contact_height
        field_string = field_string + ',' + str(f_num + 1)
        sal_face_list = [sal[i+1].id[0:sal[i+1].id.index('[')]] + [x[num_contacts+num_fascicles*2+2+i+1][0].id for x in v]
        geom.add_raw_code('''
Field['''+str(f_num)+'''] = Distance;
Field['''+str(f_num)+'''].FacesList = {'''+','.join(sal_face_list)+'''};

Field['''+str(f_num+1)+'''] = Threshold;
Field['''+str(f_num+1)+'''].IField = '''+str(f_num)+''';
Field['''+str(f_num+1)+'''].LcMin = '''+str(f_lc_min)+''';
Field['''+str(f_num+1)+'''].LcMax = '''+str(lc_nerve)+''';
Field['''+str(f_num+1)+'''].DistMin = 0;
Field['''+str(f_num+1)+'''].DistMax = '''+str(contact_height)+''';
Field['''+str(f_num+1)+'''].StopAtDistMax = 1;

''')
        
final_num = 2 + 2 * len(fascicle_diameter_x) + 2 * len(contact_angles) + 1
if round(nerve_diameter, 5) == round((cuff_inner_diameter - 2 * contact_height), 5) and len(sal) > 1:
    final_num = 2 + 2 * len(fascicle_diameter_x) + 2 * len(contact_angles) + 2 * num_isolated + 1
geom.add_raw_code('''
Field['''+str(final_num)+'''] = Min;
Field['''+str(final_num)+'''].FieldsList = {'''+field_string+'''};
Background Field = '''+str(final_num)+''';
'''
)

# %% Save
#Write Mesh File
f = open('mesh.geo','w')
f.write(geom.get_code());
f.close();

os.system('sudo gmsh/bin/gmsh -3 mesh.geo') #Create mesh using (sudo path/to/gmsh -dimension filename)
msh = meshio.read('mesh.msh') #Read mesh into meshio
meshio.write('mesh.xdmf', meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]})) #Convert to xdmf
meshio.write('mesh_physical_region.xdmf', meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]},
                 cell_data={"tetra": {"name_to_read": msh.cell_data["tetra"]["gmsh:physical"]}}))
meshio.write('mesh_facet_region.xdmf', meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                 cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))