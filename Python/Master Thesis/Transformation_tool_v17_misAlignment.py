# -*- coding: utf-8 -*-
"""

@author: Alexander Jorud, Carl Jacobsson
"""


"""
-----------------------FUTURE IMPROVEMENTS-----------------------

* Implement all visualization to plotly
    -3Dquiver does not yet exist in plotly but there are alternative methods which are easy to implement
* Use an event based gui such as pyqt5/qtdesigner
* Implement into Django
* Improved progress bar/ also a progress bar for when writing to file.
* Include more decimals regarding the local_coord_sys.nas file from Ansa....
    A small error is created due to this, although it is in the range of 5e-05 both for
    magnitude,x,y,z when comparing in Meta... 
* Implement Class structure/divide into different .py files
* Plotly plotel visualization can be greatly improved and generalized if wanted.
    -E.g. use the Open_Local_coord_sys method to open a plotel file. Only small adjustments to the code is needed.
 
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import csv
from easygui import *
import sys
import plotly
import plotly.graph_objs as go
from progressbar import *

"""
------------Manual of how to create a .exe file------------

Use the following command in the terminal in PyCharm to create an .exe file:

    pyinstaller myscriptsname.py

In the "dist" folder located in the pycharm workspace, the executable file can be located.
If a single .exe file is wanted instead of numerous of files inside a folder, the following command can be used:

    pyinstaller --onefile myscriptsname.py



!!!!! The plotly visualization module will not work with --onefile !!!!!



To enable plotly with pyinstaller:

    Don't use --onefile when creating the .exe file.
    Completely copy the plotly package from:  Lib\site-packages\plotly 
    directory into the /dist/{exe name}/ directory

"""

# Method to clear the print window after every run
def cls():
    os.system('CLS')


# Method which reads a pch file and extracts the displacements and timesteps which
# can be used during the main method. Note only translational displacements are extracted.
# The rotational displacements can be extracted with minor alterations to the code.
def read_pch(filename):

    pch_file = open(filename,'r')

    fo = pch_file.readlines()
    pch_file.close()
    
    Timestep_list = []
    X_tr_list = []
    Y_tr_list = []
    Z_tr_list = []
    X_rot_list = []
    Y_rot_list = []
    Z_rot_list = []


    msg_Globalcoord = "Reading from pch file to ensure validity of the coordinate transformation program. \n " \
                "Choose which global coordinate displacements to transform to local displacements from pch-file."
    title = "Coordinate Transformation"
    choices_Globalcoord = ["101001", "102001", "107001", "108001", "109001", "110001",
                           "115001","116001","119001","120001","121001","122001",
                           "123001","124001"]
    choice_Globalcoord_tmp = choicebox(msg_Globalcoord, title, choices_Globalcoord)
    choice_Globalcoord = str(int(choice_Globalcoord_tmp)+1000)

    counter_tmp = 0
    for line in fo:

        words = line.split("\n")        
                 
        if len(words) == 2:
            
            del words[1]

        for element in words:
            str_part = element.split(' ')

            if str_part[0] == '$POINT':
                pointID = str_part[8]

            if "$" in str_part[0]:
                break

            while ("") in str_part: str_part.remove("")


            if str_part[0] == '-CONT-':  # Storing translation displacements in lists
                X_rot_list.append(str_part[1])
                Y_rot_list.append(str_part[2])
                Z_rot_list.append(str_part[3])
            else:                        # Storing rotational displacements in lists
                Timestep_list.append(str_part[0])
                X_tr_list.append(str_part[2])
                Y_tr_list.append(str_part[3])
                Z_tr_list.append(str_part[4])

        if words[0].find('$POINT') == 0:

            if pointID >= choice_Globalcoord and counter_tmp == 0:

                counter_tmp = 1

                Timestep_matrix = np.array([Timestep_list]).astype(np.float32)

                Translation_matrix = np.matrix([np.asarray(X_tr_list).astype(np.float32),
                                                np.asarray(Y_tr_list).astype(np.float32),
                                                np.asarray(Z_tr_list).astype(np.float32)])

                # Rotation_matrix_xyz = np.matrix([np.array(X_rot_list).astype(np.float32),
                #                                  np.array(Y_rot_list).astype(np.float32),
                #                                  np.array(Z_rot_list).astype(np.float32)])

            Timestep_list = []
            X_tr_list = []
            Y_tr_list = []
            Z_tr_list = []
            X_rot_list = []
            Y_rot_list = []
            Z_rot_list = []

    return {'Tra_mat':Translation_matrix,
            'Time_mat':Timestep_matrix,'choice_Globalcoord':choice_Globalcoord_tmp}


# Reads a .csv-file and searches for numbers. The .csv-file-values needs to
# be delimited by a ","-sign
# Returns "data[Accelerometer,timeStepNumber,variable(x,y,z as 0,1,2)]" as a numpy array
def readCSV(filename):
    All = []
    timeX = []
    timeY = []
    timeZ = []
    file = open(filename, "r")
    reader = csv.reader(file, delimiter=',')
    for line in reader:
        if not line[0] == "":
            if not line[0][0] in ["#", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p",
                                  "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "A", "B", "C", "D",
                                  "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                                  "V", "W", "X", "Y", "Z"]:
                counter = 2 # column to start reading in csv
                length = len(line)
                accX = []
                accY = []
                accZ = []
                while (counter < length - 2):
                    accX.append(line[counter])
                    accY.append(line[counter + 1])
                    accZ.append(line[counter + 2])
                    counter = counter + 3

                timeX.append(accX)
                timeY.append(accY)
                timeZ.append(accZ)

    All.append(timeX)
    All.append(timeY)
    All.append(timeZ)

    var, step, acc = len(All), len(All[0]), len(All[0][0])
    data = [[[0.0 for x in range(var)] for y in range(step)] for z in range(acc)]
    
    # Rearranging
    for variable in range(0, len(All)):
        for timeStep in range(0, len(All[variable])):
            for accelerometer in range(0, len(All[variable][timeStep])):
                data[accelerometer][timeStep][variable] = All[variable][timeStep][accelerometer]
   
    if ynbox('Fix misalignment of roof accelerometers?', 'Coordinate transformation', ('Yes', 'No')):
        data=fixAlignment(data) #Fixes the misalignment of the roof accelerometers

        print 'Roof accelerometer "local" displacements transformed to "global" displacements.'

    reader_time = csv.reader(open(filename, "rb"), delimiter=",") #To get time data from the csv-file
    tmp = list(reader_time)
    timestep_list = []

    for p in range(8,len(tmp)):
        timestep_list.append(tmp[p][1])

    return {'disp_mat': np.asarray(data).astype(np.float32), 'timestep_mat': timestep_list}


# Fixes the misalignment of the roof accelerometers.
def fixAlignment(data):

    angle_x_front=np.deg2rad(-9.42)   # Anti-clockwise rotation around front x-axis
    angle_y_front=np.deg2rad(9.64)    # Anti-clockwise rotation around front y-axis
    angle_z_front=np.deg2rad(1.37)    # Clockwise rotation around front z-axis
     
    angle_x_rear=np.deg2rad(5.47)   # Rotation around rear x-axis
    angle_y_rear=np.deg2rad(-4.25) # Rotation around rear y-axis
    angle_z_rear=np.deg2rad(-6.9)   # Rotation around rear z-axis
     
    frontPos=9 # location of the used front roof accelerometer in the data list
    rearPos=8 # location of the used rear roof accelerometer in the data list
    
    def invrotation(angle_x, angle_y, angle_z):
         #Roe Convention - 3-D rotaton matrix / Euler angles
         # http://mathworld.wolfram.com/EulerAngles.html 
         # https: // www.cs.utexas.edu / ~theshark / courses / cs354 / lectures / cs354 - 14.
         # pdf (Followed this convention)
         X = np.matrix([[1, 0, 0],
                      [0, np.cos(angle_x), -np.sin(angle_x)],
                       [0, np.sin(angle_x), np.cos(angle_x)]])

         Y = np.matrix([[np.cos(angle_y), 0, np.sin(angle_y)],
                      [0, 1, 0],
                       [-np.sin(angle_y), 0, np.cos(angle_y)]])

         Z = np.matrix([[np.cos(angle_z), -np.sin(angle_z), 0],
                      [np.sin(angle_z), np.cos(angle_z), 0],
                       [0, 0, 1]])

         R = Z.dot(Y).dot(X)
         return np.linalg.inv(R)
    
    invRfront=invrotation(angle_x_front, angle_y_front, angle_z_front)
    invRrear=invrotation(angle_x_rear, angle_y_rear, angle_z_rear)
     
    for i in range (0, len(data[frontPos])):         
        globalFront=invRfront.dot(np.transpose(np.array([float(data[frontPos][i][0]), float(data[frontPos][i][1]), float(data[frontPos][i][2])])))
        globalRear=invRrear.dot(np.transpose(np.array([float(data[rearPos][i][0]), float(data[rearPos][i][1]), float(data[rearPos][i][2])])))
        for j in range (0,3):
            data[frontPos][i][j]=globalFront[0,j]
            data[rearPos][i][j]=globalRear[0,j]
            
    plt.close('all') #Close all previous open figures.
    figRear = plt.figure()
    axRear = figRear.gca(projection='3d')
    plt.title(r'Acc. misalignment - CORD_SYS:'+' 219001', fontsize=14)
    figFront = plt.figure()
    axFront = figFront.gca(projection='3d')
    plt.title(r'Acc. misalignment - CORD_SYS:'+' 220001', fontsize=14)

    #Create Reference orthogonoal basis in x,y,z=0,0,0
    refOrthBasis = np.matrix([[1,0,0], [0,1,0],[0,0,1]])
    Ref_orth_Basis_pos = [0,0,0]
    Grid_origin = [0,0,0]
    x = 100
    
    axRear.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], refOrthBasis[0,:], refOrthBasis[1,:], refOrthBasis[2,:],length=100, pivot='tail')
    axRear.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], invRrear.dot(refOrthBasis)[0,:], invRrear.dot(refOrthBasis)[1,:], invRrear.dot(refOrthBasis)[2,:],length=100, pivot='tail', color='red')
    axRear.scatter(Ref_orth_Basis_pos[0], Ref_orth_Basis_pos[1], Ref_orth_Basis_pos[2] , color='red')
    axFront.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], refOrthBasis[0,:], refOrthBasis[1,:], refOrthBasis[2,:],length=100, pivot='tail')
    axFront.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], invRfront.dot(refOrthBasis)[0,:], invRfront.dot(refOrthBasis)[1,:], invRfront.dot(refOrthBasis)[2,:],length=100, pivot='tail', color='red')
    axFront.scatter(Ref_orth_Basis_pos[0], Ref_orth_Basis_pos[1], Ref_orth_Basis_pos[2] , color='red')    
    
    axRear.set_xlim(Grid_origin[0]-x, Grid_origin[0]+x)
    axRear.set_ylim(Grid_origin[1]-x, Grid_origin[1]+x)
    axRear.set_zlim(Grid_origin[2]-x, Grid_origin[2]+x)
    axFront.set_xlim(Grid_origin[0]-x, Grid_origin[0]+x)
    axFront.set_ylim(Grid_origin[1]-x, Grid_origin[1]+x)
    axFront.set_zlim(Grid_origin[2]-x, Grid_origin[2]+x)

    axRear.set_xlabel('X axis')
    axRear.set_ylabel('Y axis')
    axRear.set_zlabel('Z axis')
    axFront.set_xlabel('X axis')
    axFront.set_ylabel('Y axis')
    axFront.set_zlabel('Z axis')

    plt.show()
            
    return data


# Transforming global displacements to local displacements.
def transformation(Local_coord, Global_coord):

    Transformation_matrix = Local_coord.dot((np.linalg.det(Global_coord))*np.transpose(Global_coord))
    check_if_orthogonal(Local_coord)

    return Transformation_matrix


# Visualizing method for the pch file
def visualization_pch(Local_coord,Global_coord,Grid_origin,Grid_point_selected):

    fig = plt.figure()
    #    Visualization of the local coordinate system and the global coordinate system

    ax = fig.gca(projection='3d')
    plt.title('Grid point:'+Grid_point_selected, fontsize=14)

    ax.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], Global_coord[:, 0],
              Global_coord[:, 1], Global_coord[:, 2], length=100, pivot='tail', color='black')

    ax.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], Local_coord[0, 0],
              Local_coord[0, 1], Local_coord[0, 2], length=100, pivot='tail', color='blue')

    ax.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], Local_coord[1, 0],
              Local_coord[1, 1], Local_coord[1, 2], length=100, pivot='tail', color='blue')

    ax.quiver(Grid_origin[0], Grid_origin[1], Grid_origin[2], Local_coord[2, 0],
              Local_coord[2, 1], Local_coord[2, 2], length=100, pivot='tail', color='red')

    print 'Legends for the visualization:'
    print '---------------------------------------'
    print 'Black axes: Global coordinate system.'
    print 'Red axis:   Local z-axis.'
    print 'Blue axes: Local x and y-axis.'
    print '---------------------------------------'

    x = 100
    ax.set_xlim(Grid_origin[0] - x, Grid_origin[0] + x)
    ax.set_ylim(Grid_origin[1] - x, Grid_origin[1] + x)
    ax.set_zlim(Grid_origin[2] - x, Grid_origin[2] + x)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


# Visualizing method for the csv file
def visualization_csv(Local_coord_list,Global_coord_list,Grid_origin_list):


    title = "Coordinate Transformation"

    # The following while statement and its contents enables the user to visualize multiple times.
    check_counter = 0
    while True:
        if check_counter > 0:
            check = ynbox('Visualize again?.', title, ('Yes', 'No'))
            if not check:
                break

        msg_Globalcoord = "Choose visualization option."
        choices_vis = ["Show all evaluated grid point locations", "Show specific grid point location",
                       "Continue without any visualization."]
        choice_vis = choicebox(msg_Globalcoord, title, choices_vis)


        print 'Legends for the visualization:'
        print '------------------------------'
        print 'Black axes: Global coordinate system.'
        print 'Red axis:   Local z-axis.'
        print 'Blue axes: Local x and y-axis.'

        # Visualization of all evaluated grid points, both the local coordinate system and global coordinate system is
        # included.
        if choice_vis == "Show all evaluated grid point locations":

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            plt.title('All grid point locations evaluated', fontsize=14)

            for j in range(0, len(Local_coord_list)):
                ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                          Global_coord_list[j][:, 0], Global_coord_list[j][:, 1], Global_coord_list[j][:, 2],
                          length=100, pivot='tail', color='black')

                ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                          Local_coord_list[j][0, 0], Local_coord_list[j][0, 1], Local_coord_list[j][0, 2],
                          length=100, pivot='tail', color='blue')

                ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                          Local_coord_list[j][1, 0], Local_coord_list[j][1, 1], Local_coord_list[j][1, 2],
                          length=100, pivot='tail', color='blue')

                ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                          Local_coord_list[j][2, 0], Local_coord_list[j][2, 1], Local_coord_list[j][2, 2],
                          length=100, pivot='tail', color='red')



            ax.set_xlim(2400, 5600)
            ax.set_ylim(-900, 900)
            ax.set_zlim(300, 1800)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            plt.show()


        # Visualization of a chosen local coordinate system, included is also a global coordinate system.
        elif choice_vis == "Show specific grid point location":

            msg_acc = "Choose which accelerometer to view."
            title = "Coordinate Transformation"
            choices_acc = ["201001", "202001", "207001", "208001", "209001", "210001", "215001", "216001",
                           "219001", "220001", "221001", "222001", "223001", "224001"]
            choice_acc = choicebox(msg_acc, title, choices_acc)

            accelerometer_ID = generate_grid_points()['accelerometers']

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            plt.title('Grid point: ' + choice_acc, fontsize=14)

            for z in range(0, len(accelerometer_ID)):

                if choice_acc == str(accelerometer_ID[z]):
                    j = z

            ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                      Global_coord_list[j][:, 0], Global_coord_list[j][:, 1], Global_coord_list[j][:, 2],
                      length=100, pivot='tail', color='black')

            ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                      Local_coord_list[j][0, 0], Local_coord_list[j][0, 1], Local_coord_list[j][0, 2],
                      length=100, pivot='tail', color='blue')

            ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                      Local_coord_list[j][1, 0], Local_coord_list[j][1, 1], Local_coord_list[j][1, 2],
                      length=100, pivot='tail', color='blue')

            ax.quiver(Grid_origin_list[j][0], Grid_origin_list[j][1], Grid_origin_list[j][2],
                      Local_coord_list[j][2, 0], Local_coord_list[j][2, 1], Local_coord_list[j][2, 2],
                      length=100, pivot='tail', color='red')

            x = 100
            ax.set_xlim(Grid_origin_list[j][0] - x, Grid_origin_list[j][0] + x)
            ax.set_ylim(Grid_origin_list[j][1] - x, Grid_origin_list[j][1] + x)
            ax.set_zlim(Grid_origin_list[j][2] - x, Grid_origin_list[j][2] + x)

            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')

            plt.show()

        elif choice_vis == "Continue without any visualization.":
            break

        check_counter += 1


# Method to check that the imported coordinate system are orthogonal.
def check_if_orthogonal(Local_coord):

    Temp_matrix = np.absolute(Local_coord.dot(np.transpose(Local_coord)))

    pos01 = Temp_matrix[0,1]
    pos02 = Temp_matrix[0,2]
    pos10 = Temp_matrix[1,0]
    pos12 = Temp_matrix[1,2]
    pos20 = Temp_matrix[2,0]
    pos21 = Temp_matrix[2,1]

    # tolerance
    TOL = 0.0001
    title = "Orthogonality basis check."

    if pos01 > TOL and pos02 > TOL and pos10 > TOL and pos12 > TOL and pos20 > TOL and pos21 > TOL:
        check = ynbox('Local coordinate system is not a orthonormal basis. \n '
                      'Shall I continue?', title, ('Yes', 'No'))

        if not check:
            sys.exit(0)


# Reads a file including the grid points which defines a coordinate system from Ansa.
def Open_Local_coord_sys():
    msgbox(msg="Choose a local_coord_sys.nas file"
               + "\n" +
               "\n" +
               "\n" +
               "-------------------------------------------------"
               + "\n" +
               "An example of a local_coord_sys.nas file:"
               + "\n"
               + "\n" +
               "GRID      101001        2397.692-750.801519.1727"
               + "\n" +
               "GRID      201001        2397.709-750.801519.2727"
               + "\n" +
               "GRID      202001        3342.936-594.0451565.543"
               + "\n" +
               "CORD1R    201001  201001  202001  101001"
               + "\n" +
               "-------------------------------------------------",
           title="Coordinate Transformation")
    filename = fileopenbox()
    localcoordsys = open(filename, 'r')
    fo = localcoordsys.readlines()
    localcoordsys.close()


    CORD1R_list = []
    G1_list = []
    G2_list = []
    G3_list = []
    GRID_x = []
    GRID_y = []
    GRID_z = []
    t = 8
    GRID_ID = []

    for line in fo:
        words = line.split("\n")

        if len(words) == 2:
            del words[1]


        for element in words:


            str_part = element.split(' ')
            str_part_tmp = map(''.join, zip(*[iter(element)] * 8))


            if "$" in str_part[0]:
                break


            while ("") in str_part: str_part.remove("")



            if str_part[0] == 'GRID':

                gX = float(str_part_tmp[3])
                gY = float(str_part_tmp[4])
                gZ = float(str_part_tmp[5])

                GRID_ID.append(str_part[1])
                GRID_x.append(gX)
                GRID_y.append(gY)
                GRID_z.append(gZ)


            if str_part[0] == 'CORD1R':

                CORD1R_list.append(str_part[1])
                G1_list.append(str_part[2])
                G2_list.append(str_part[3])
                G3_list.append(str_part[4])




    CORD1R_list_VIS = [CORD1R_list, G1_list, G2_list, G3_list]
    GRID_list = [GRID_ID, GRID_x, GRID_y, GRID_z]


    Temp_list = CORD1R_list_VIS


    for h in range(0,len(GRID_x)):
        for k in range(0,len(G1_list)):
            for l in range(1,4):
                if GRID_list[0][h] == CORD1R_list_VIS[l][k]:
                    Temp_list[l][k] = [GRID_list[1][h], GRID_list[2][h], GRID_list[3][h]]

    msg_Globalcoord = "Choose visualization option. \nFor better performance open the .html files in Google Chrome."
    choices_vis = ["Visualize imported grid points from the " + filename + " file",
                   "Continue without any visualization."]
    choice_vis = choicebox(msg_Globalcoord, "Coordinate Transformation", choices_vis)

    if choice_vis == "Visualize imported grid points from the " + filename + " file":
        trace1 = go.Scatter3d(
            x=GRID_x,
            y=GRID_y,
            z=GRID_z,
            mode='markers',
            text=GRID_ID,
            name='Grid ID',
            textposition='top',
            hoverinfo='x,y,z',
            marker=dict(
                size=12,
                symbol='circle',
                line=dict(
                    color='rgba(217, 217, 217, 0.14)',
                    width=0.5
                ),
                opacity=0.8
            )
        )

        # E-points visualization
        trace2 = go.Scatter3d(
            x=[Temp_list[1][0][0], Temp_list[2][0][0]],
            y=[Temp_list[1][0][1], Temp_list[2][0][1]],
            z=[Temp_list[1][0][2], Temp_list[2][0][2]],
            mode='lines',
            name="E-Line: " + Temp_list[0][0],
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(255, 0, 0, 0.8)',
            )
        )

        trace3 = go.Scatter3d(
            x=[Temp_list[1][1][0], Temp_list[2][1][0]],
            y=[Temp_list[1][1][1], Temp_list[2][1][1]],
            z=[Temp_list[1][1][2], Temp_list[2][1][2]],
            mode='lines',
            name="E-Line: " + Temp_list[0][1],
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(255, 0, 0, 0.8)',
            )
        )

        trace4 = go.Scatter3d(
            x=[Temp_list[1][2][0], Temp_list[2][2][0]],
            y=[Temp_list[1][2][1], Temp_list[2][2][1]],
            z=[Temp_list[1][2][2], Temp_list[2][2][2]],
            mode='lines',
            name="E-Line: " + Temp_list[0][2],
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(255, 0, 0, 0.8)',
            )
        )

        trace5 = go.Scatter3d(
            x=[Temp_list[1][3][0], Temp_list[2][3][0]],
            y=[Temp_list[1][3][1], Temp_list[2][3][1]],
            z=[Temp_list[1][3][2], Temp_list[2][3][2]],
            mode='lines',
            name="E-Line: " + Temp_list[0][3],
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(255, 0, 0, 0.8)',
            )
        )

        trace6 = go.Scatter3d(
            x=[Temp_list[1][4][0], Temp_list[2][4][0]],
            y=[Temp_list[1][4][1], Temp_list[2][4][1]],
            z=[Temp_list[1][4][2], Temp_list[2][4][2]],
            mode='lines',
            name="E-Line: " + Temp_list[0][4],
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(255, 0, 0, 0.8)',
            )
        )

        trace7 = go.Scatter3d(
            x=[Temp_list[1][5][0], Temp_list[2][5][0]],
            y=[Temp_list[1][5][1], Temp_list[2][5][1]],
            z=[Temp_list[1][5][2], Temp_list[2][5][2]],
            mode='lines',
            name="E-Line: " + Temp_list[0][5],
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(255, 0, 0, 0.8)',
            )
        )

        trace8 = go.Scatter3d(
            x=[Temp_list[1][6][0], Temp_list[2][6][0]],
            y=[Temp_list[1][6][1], Temp_list[2][6][1]],
            z=[Temp_list[1][6][2], Temp_list[2][6][2]],
            mode='lines',
            name="E-Line: " + Temp_list[0][6],
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(255, 0, 0, 0.8)',
            )
        )

        # Plotel visualization

        GRID_10481445 = [1124.818, 467.751, 902.188]
        GRID_10522797 = [5029.582, 816.5411, 901.7781]
        GRID_10526596 = [2314.639, 815.4728, 427.7071]
        GRID_10527935 = [3189.085, 812.6772, 410.7405]
        GRID_10529217 = [3911.622, 816.2675, 431.1044]
        GRID_10533256 = [4982.592, 712.9239, 1270.576]
        GRID_10537597 = [2187.324, 783.8357, 1187.428]
        GRID_10541677 = [2848.693, 637.2427, 1605.562]
        GRID_10542343 = [3385.78, 612.0434, 1666.176]
        GRID_10544297 = [4517.183, 606.3771, 1618.336]
        GRID_10679597 = [1124.091, -467.751, 903.6284]
        GRID_10692276 = [5028.411, -817.155, 901.0937]
        GRID_10696517 = [4984.58, -713.495, 1269.758]
        GRID_10705864 = [2314.057, -815.155, 425.738]
        GRID_10707422 = [3190.829, -813.302, 414.2796]
        GRID_10708919 = [3909.237, -817.076, 435.6888]
        GRID_10715201 = [2190.708, -783.564, 1189.418]
        GRID_10726112 = [2850.041, -636.197, 1606.592]
        GRID_10726706 = [3389.365, -612.29, 1665.938]
        GRID_10728359 = [4519.979, -606.945, 1617.436]
        GRID_10827365 = [1108.742, -557.05, 593.5827]
        GRID_10828857 = [1114.465, 567.05, 602.5099]

        trace9 = go.Scatter3d(
            x=[GRID_10827365[0], GRID_10705864[0]
                , GRID_10705864[0], GRID_10707422[0]
                , GRID_10707422[0], GRID_10708919[0]
                , GRID_10708919[0], GRID_10692276[0]
                , GRID_10692276[0], GRID_10696517[0]
                , GRID_10696517[0], GRID_10728359[0]
                , GRID_10728359[0], GRID_10726706[0]
                , GRID_10726706[0], GRID_10726112[0]
                , GRID_10726112[0], GRID_10715201[0]
                , GRID_10715201[0], GRID_10679597[0]
                , GRID_10679597[0], GRID_10827365[0]],

            y=[GRID_10827365[1], GRID_10705864[1]
                , GRID_10705864[1], GRID_10707422[1]
                , GRID_10707422[1], GRID_10708919[1]
                , GRID_10708919[1], GRID_10692276[1]
                , GRID_10692276[1], GRID_10696517[1]
                , GRID_10696517[1], GRID_10728359[1]
                , GRID_10728359[1], GRID_10726706[1]
                , GRID_10726706[1], GRID_10726112[1]
                , GRID_10726112[1], GRID_10715201[1]
                , GRID_10715201[1], GRID_10679597[1]
                , GRID_10679597[1], GRID_10827365[1]],

            z=[GRID_10827365[2], GRID_10705864[2]
                , GRID_10705864[2], GRID_10707422[2]
                , GRID_10707422[2], GRID_10708919[2]
                , GRID_10708919[2], GRID_10692276[2]
                , GRID_10692276[2], GRID_10696517[2]
                , GRID_10696517[2], GRID_10728359[2]
                , GRID_10728359[2], GRID_10726706[2]
                , GRID_10726706[2], GRID_10726112[2]
                , GRID_10726112[2], GRID_10715201[2]
                , GRID_10715201[2], GRID_10679597[2]
                , GRID_10679597[2], GRID_10827365[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace10 = go.Scatter3d(
            x=[GRID_10481445[0], GRID_10537597[0]
                , GRID_10537597[0], GRID_10541677[0]
                , GRID_10541677[0], GRID_10542343[0]
                , GRID_10542343[0], GRID_10544297[0]
                , GRID_10544297[0], GRID_10533256[0]
                , GRID_10533256[0], GRID_10522797[0]
                , GRID_10522797[0], GRID_10529217[0]
                , GRID_10529217[0], GRID_10527935[0]
                , GRID_10527935[0], GRID_10526596[0]
                , GRID_10526596[0], GRID_10828857[0]
                , GRID_10828857[0], GRID_10481445[0]],

            y=[GRID_10481445[1], GRID_10537597[1]
                , GRID_10537597[1], GRID_10541677[1]
                , GRID_10541677[1], GRID_10542343[1]
                , GRID_10542343[1], GRID_10544297[1]
                , GRID_10544297[1], GRID_10533256[1]
                , GRID_10533256[1], GRID_10522797[1]
                , GRID_10522797[1], GRID_10529217[1]
                , GRID_10529217[1], GRID_10527935[1]
                , GRID_10527935[1], GRID_10526596[1]
                , GRID_10526596[1], GRID_10828857[1]
                , GRID_10828857[1], GRID_10481445[1]],

            z=[GRID_10481445[2], GRID_10537597[2]
                , GRID_10537597[2], GRID_10541677[2]
                , GRID_10541677[2], GRID_10542343[2]
                , GRID_10542343[2], GRID_10544297[2]
                , GRID_10544297[2], GRID_10533256[2]
                , GRID_10533256[2], GRID_10522797[2]
                , GRID_10522797[2], GRID_10529217[2]
                , GRID_10529217[2], GRID_10527935[2]
                , GRID_10527935[2], GRID_10526596[2]
                , GRID_10526596[2], GRID_10828857[2]
                , GRID_10828857[2], GRID_10481445[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace11 = go.Scatter3d(
            x=[GRID_10481445[0], GRID_10679597[0]],

            y=[GRID_10481445[1], GRID_10679597[1]],

            z=[GRID_10481445[2], GRID_10679597[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace12 = go.Scatter3d(
            x=[GRID_10828857[0], GRID_10827365[0]],

            y=[GRID_10828857[1], GRID_10827365[1]],

            z=[GRID_10828857[2], GRID_10827365[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace13 = go.Scatter3d(
            x=[GRID_10537597[0], GRID_10715201[0]],

            y=[GRID_10537597[1], GRID_10715201[1]],

            z=[GRID_10537597[2], GRID_10715201[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace14 = go.Scatter3d(
            x=[GRID_10522797[0], GRID_10692276[0]],

            y=[GRID_10522797[1], GRID_10692276[1]],

            z=[GRID_10522797[2], GRID_10692276[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace15 = go.Scatter3d(
            x=[GRID_10828857[0], GRID_10827365[0]],

            y=[GRID_10828857[1], GRID_10827365[1]],

            z=[GRID_10828857[2], GRID_10827365[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace16 = go.Scatter3d(
            x=[GRID_10696517[0], GRID_10533256[0]],

            y=[GRID_10696517[1], GRID_10533256[1]],

            z=[GRID_10696517[2], GRID_10533256[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace17 = go.Scatter3d(
            x=[GRID_10728359[0], GRID_10544297[0]],

            y=[GRID_10728359[1], GRID_10544297[1]],

            z=[GRID_10728359[2], GRID_10544297[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace18 = go.Scatter3d(
            x=[GRID_10705864[0], GRID_10526596[0]],

            y=[GRID_10705864[1], GRID_10526596[1]],

            z=[GRID_10705864[2], GRID_10526596[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace19 = go.Scatter3d(
            x=[GRID_10708919[0], GRID_10529217[0]],

            y=[GRID_10708919[1], GRID_10529217[1]],

            z=[GRID_10708919[2], GRID_10529217[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        trace20 = go.Scatter3d(
            x=[GRID_10541677[0], GRID_10726112[0]],

            y=[GRID_10541677[1], GRID_10726112[1]],

            z=[GRID_10541677[2], GRID_10726112[2]],
            mode='lines',
            name="Plotel",
            hoverinfo='name',
            line=dict(
                width=5,
                color='rgba(0, 0, 0, 0.6)',
            )
        )

        data = [trace1, trace2, trace3, trace4, trace5, trace6, trace7, trace8, trace9, trace10, trace11, trace12,
                trace13,
                trace14, trace15, trace16, trace17, trace18, trace19, trace20]
        layout = go.Layout(
            title='Imported grid points and E-points from Ansa. Distances in [mm]',
            scene=dict(
                camera=dict(
                    center=dict(
                        x=0,
                        y=0,
                        z=0
                    ),
                    eye=dict(
                        x=1.96903462608,
                        y=-1.09022831971,
                        z=0.405345349304
                    ),
                    up=dict(
                        x=0,
                        y=0,
                        z=1
                    )
                ),
                dragmode="turntable",
                xaxis=dict(
                    title="x",
                    type="x",
                    titlefont=dict(
                        family='Courier New, monospace',
                        size=18,
                        color='#7f7f7f'
                    )
                ),
                yaxis=dict(
                    title="y",
                    type="y",
                    titlefont=dict(
                        family='Courier New, monospace',
                        size=18,
                        color='#7f7f7f'
                    )
                ),
                zaxis=dict(
                    title="z",
                    type="z",
                    titlefont=dict(
                        family='Courier New, monospace',
                        size=18,
                        color='#7f7f7f'
                    )
                ),
                annotations=[
                    dict(
                        showarrow=False,
                        x=GRID_x[0],
                        y=GRID_y[0],
                        z=GRID_z[0],
                        text=GRID_ID[0],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[1],
                        y=GRID_y[1],
                        z=GRID_z[1],
                        text=GRID_ID[1],
                        xanchor="left",
                        xshift=10,
                        yshift=-15,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[2],
                        y=GRID_y[2],
                        z=GRID_z[2],
                        text=GRID_ID[2],
                        xanchor="left",
                        xshift=10,
                        yshift=-15,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[3],
                        y=GRID_y[3],
                        z=GRID_z[3],
                        text=GRID_ID[3],
                        xanchor="left",
                        xshift=10,
                        yshift=-15,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[4],
                        y=GRID_y[4],
                        z=GRID_z[4],
                        text=GRID_ID[4],
                        xanchor="left",
                        xshift=10,
                        yshift=-15,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[5],
                        y=GRID_y[5],
                        z=GRID_z[5],
                        text=GRID_ID[5],
                        xanchor="left",
                        xshift=10,
                        yshift=-15,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[6],
                        y=GRID_y[6],
                        z=GRID_z[6],
                        text=GRID_ID[6],
                        xanchor="left",
                        xshift=10,
                        yshift=-15,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[7],
                        y=GRID_y[7],
                        z=GRID_z[7],
                        text=GRID_ID[7],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[8],
                        y=GRID_y[8],
                        z=GRID_z[8],
                        text=GRID_ID[8],
                        xanchor="left",
                        xshift=-55,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[9],
                        y=GRID_y[9],
                        z=GRID_z[9],
                        text=GRID_ID[9],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[10],
                        y=GRID_y[10],
                        z=GRID_z[10],
                        text=GRID_ID[10],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[11],
                        y=GRID_y[11],
                        z=GRID_z[11],
                        text=GRID_ID[11],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[12],
                        y=GRID_y[12],
                        z=GRID_z[12],
                        text=GRID_ID[12],
                        xanchor="left",
                        xshift=-55,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[13],
                        y=GRID_y[13],
                        z=GRID_z[13],
                        text=GRID_ID[13],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[14],
                        y=GRID_y[14],
                        z=GRID_z[14],
                        text=GRID_ID[14],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[15],
                        y=GRID_y[15],
                        z=GRID_z[15],
                        text=GRID_ID[15],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[16],
                        y=GRID_y[16],
                        z=GRID_z[16],
                        text=GRID_ID[16],
                        xanchor="left",
                        xshift=-55,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[17],
                        y=GRID_y[17],
                        z=GRID_z[17],
                        text=GRID_ID[17],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[18],
                        y=GRID_y[18],
                        z=GRID_z[18],
                        text=GRID_ID[18],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[19],
                        y=GRID_y[19],
                        z=GRID_z[19],
                        text=GRID_ID[19],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    ),
                    dict(
                        showarrow=False,
                        x=GRID_x[20],
                        y=GRID_y[20],
                        z=GRID_z[20],
                        text=GRID_ID[20],
                        xanchor="left",
                        xshift=10,
                        opacity=0.7
                    )]
            ),
        )
        fig = go.Figure(data=data, layout=layout)
        plotly.offline.plot(fig, filename='imported_ansa_coord.html')

    print 'imported_ansa_coord.html saved in working directory.'

    # Returns a list with the following content
    # Temp_list[0:end] different "grid points" or accelerometers
    #
    #                         ID     X,Y,Z   X,Y,Z   X,Y,Z
    #
    # Where Temp_list=     [201001  201001  202001  203001]
    #                      [203001  203001  204001  202001]
    #                      [205001  205001  206001  207001]
    #                      [207001  207001  208001  206001]

    # Where ID column is one list
    # Second column has a list [X,Y,Z] for the grid point 201001 for example

    return {'CORD1R_mat':Temp_list}


# Creates a local coordinate system from three grid points in space.
# Also creates a global coordinate system.
def createLocalcoord(G1,G2,G3):


    V1 = (G2 - G1) / np.linalg.norm(G2 - G1)
    V2 = (G3 - G1) / np.linalg.norm(G3 - G1)
    n = np.cross(V1, V2) / np.linalg.norm(np.cross(V1, V2))
    x_local = np.cross(n, V1) / np.linalg.norm(np.cross(n, V1))
    Local_coord = np.matrix([x_local, n, V1])


    Global_coord = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    return {'Local_coord':Local_coord, 'Global_coord':Global_coord, 'G1':G1,'G2':G2}


# Generates accelerometer grid points.
def generate_grid_points():

    accelerometers = [201001,202001,207001,208001,209001,210001,215001,216001,219001,220001,221001,222001,223001,224001]
    amount_of_accelerometers = len(accelerometers)


    return {'accelerometers': accelerometers, 'amount_of_accelerometers': amount_of_accelerometers}


def main():

    widgets = ['Progressbar: ', Percentage(), ' ', Bar(marker='#', left='[', right=']'),
               ' ', ETA()]

    pbar = ProgressBar(widgets=widgets, maxval=4)
    pbar.start()
    prbarcounter = 0

    # print 'Progress bar: [*            ]'
    prbarcounter += 1
    pbar.update(prbarcounter)


    title = "Coordinate transformation"
    msgbox(msg="Choose pch-file or csv-file.",title=title)
    filename = fileopenbox()
    print 'Open file:',filename


    #     READING PCH FILE SECTION
    if filename.endswith('.pch'):

        msgbox(msg=".pch file chosen.",
               title=title)

        File_choice = read_pch(filename)
        Disp_matrix = File_choice['Tra_mat']
        Time_step_matrix = File_choice['Time_mat']
        Grid_point_selected = File_choice['choice_Globalcoord']

        # print 'Progress bar: [***          ]'
        prbarcounter += 1
        pbar.update(prbarcounter)  # adds a symbol at each iteration

        CORD1R = Open_Local_coord_sys()['CORD1R_mat']

        # Extracts the correct grid points for the selected coordinate system.
        for k in range(0,len(CORD1R[0])):
            if CORD1R[0][k] == str((int(Grid_point_selected)+100000)):
                G1 = np.asarray(CORD1R[1][k]).astype(np.float32)
                G2 = np.asarray(CORD1R[2][k]).astype(np.float32)
                G3 = np.asarray(CORD1R[3][k]).astype(np.float32)
                CORD1R_string = CORD1R[0][k]


        try:
            G1
        except NameError:

            acc_list = generate_grid_points()['accelerometers']
            del acc_list[::2]

            for h in range(0,len(CORD1R[0])):

                if (int(Grid_point_selected) + 100000) == acc_list[h]:
                    G1 = np.asarray(CORD1R[1][h]).astype(np.float32)
                    G2 = np.asarray(CORD1R[2][h]).astype(np.float32)
                    G3 = np.asarray(CORD1R[3][h]).astype(np.float32)
                    CORD1R_string = CORD1R[0][h]

        # Creates the local coordinate system from the extracted grid points G1,G2 and G3.
        Local_coord = createLocalcoord(G1,G2,G3)['Local_coord']
        Global_coord = createLocalcoord(G1,G2,G3)['Global_coord']
        Grid_origin = createLocalcoord(G1,G2,G3)['G1']


        msgbox(msg="Transforming global displacements to local displacements. \n"
                   "Grid point ID:" + Grid_point_selected + " \n"
                    "Local coordinate system:" + CORD1R_string,
               title=title)



        # Creates a transformation matrix.
        Transformation_matrix = transformation(Local_coord, Global_coord)
        print 'Performing transformation'
        visualization_pch(Local_coord, Global_coord, Grid_origin,Grid_point_selected)
        # print 'Progress bar: [**********   ]'
        prbarcounter += 1
        pbar.update(prbarcounter)  # adds a symbol at each iteration
        Transformed_disp_matrix = []

        # Multiplies the transformation matrix with the "global" displacements for each timestep to
        # obtain the new "local" displacements
        for i in range(0, Disp_matrix[0].size):
            Transformed_disp_matrix.append(Transformation_matrix.dot(Disp_matrix[:, i]))

        Grid_point_selected = str(int(Grid_point_selected) + 100000)

        # Save to file
        filesave_path = filesavebox()
        with open(filesave_path, 'w') as file:
            file.write("$TITLE   = Test_data" + "\n" + "$SUBTITLE= " + "\n" + "$LABEL   =   " + "\n" + "$DISPLACEMENTS " + "\n" + "$REAL OUTPUT " + "\n")
            file.write("$SUBCASE ID =       Event " + "\n")
            file.write("$POINT ID =      " + Grid_point_selected + "  IDENTIFIED BY TIME " + "\n")
            for i in range(0, Time_step_matrix.size):
                L = "    " + str("%.6E" % Time_step_matrix[0, i]) + " G      " + str(
                    "%.6E" % Transformed_disp_matrix[i][0]) + "      " + str(
                    "%.6E" % Transformed_disp_matrix[i][1]) + "      " + str(
                    "%.6E" % Transformed_disp_matrix[i][2]) + "\n"
                newstr = L.replace(" -", "-")
                file.write(newstr)
                file.write("-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00" + "\n")
        # print 'Progress bar: [*************]'
        prbarcounter += 1
        pbar.update(prbarcounter)  # adds a symbol at each iteration
        print 'Filed saved at:', filesave_path

        pbar.finish()
        print

    #     READING CSV FILE SECTION
    elif filename.endswith('.csv'):

        msgbox(msg=".csv file chosen.",
               title=title)

        csvOutput = readCSV(filename)
        File_choice = csvOutput['disp_mat']
        Time_step_matrix = csvOutput['timestep_mat']
        accelerometer_ID = generate_grid_points()['accelerometers']
        nbr_of_accelerometers = generate_grid_points()['amount_of_accelerometers']
        CORD1R = Open_Local_coord_sys()['CORD1R_mat']
        # print 'Progress bar: [***          ]'
        prbarcounter += 1
        pbar.update(prbarcounter)  # adds a symbol at each iteration

        Transformed_disp_matrix_Storage = []

        # Creating lists to store values for plotting
        Local_coord_list = []
        Global_coord_list = []
        Grid_origin_list = []


        # Loop over all accelerometers
        for acc in range(0,nbr_of_accelerometers):
            Disp_matrix_tmp = File_choice[acc,:]
            Disp_matrix = np.matrix([Disp_matrix_tmp[:,0],
                                    Disp_matrix_tmp[:,1],
                                    Disp_matrix_tmp[:,2]])


            print 'Performing transformation for accelerometer ID:',accelerometer_ID[acc]

            # Extracts the correct grid points for the current coordinate system in the for-loop.
            for k in range(0,len(CORD1R[0])):
                if CORD1R[0][k] == str(accelerometer_ID[acc]):
                    G1 = np.asarray(CORD1R[1][k]).astype(np.float32)
                    G2 = np.asarray(CORD1R[2][k]).astype(np.float32)
                    G3 = np.asarray(CORD1R[3][k]).astype(np.float32)
                    Grid_origin_list.append(G1)
                elif CORD1R[0][k] == str(accelerometer_ID[acc]-1000):
                    G1 = np.asarray(CORD1R[1][k]).astype(np.float32)
                    G2 = np.asarray(CORD1R[2][k]).astype(np.float32)
                    G3 = np.asarray(CORD1R[3][k]).astype(np.float32)
                    Grid_origin_list.append(G2)


            Local_coord = createLocalcoord(G1,G2,G3)['Local_coord']
            Global_coord = createLocalcoord(G1,G2,G3)['Global_coord']


            # Storing values for plotting
            Local_coord_list.append(Local_coord)
            Global_coord_list.append(Global_coord)


            Transformation_matrix = transformation(Local_coord, Global_coord)

            Transformed_disp_matrix = []

            # Multiplies the transformation matrix with the "global" displacements for each timestep to
            # obtain the new "local" displacements
            for i in range(0, Disp_matrix[0].size):
                Transformed_disp_matrix.append(Transformation_matrix.dot(Disp_matrix[:,i]))

            Transformed_disp_matrix_Storage.append(Transformed_disp_matrix)

        # print 'Progress bar: [**********   ]'
        prbarcounter += 1
        pbar.update(prbarcounter)  # adds a symbol at each iteration


        # Visualization
        visualization_csv(Local_coord_list, Global_coord_list, Grid_origin_list)


        # Save to file
        filesave_path = filesavebox()
        with open(filesave_path, 'w') as file:
            for k in range(0,nbr_of_accelerometers):
                file.write("$TITLE   = Test_data" + "\n" + "$SUBTITLE= " + "\n" + "$LABEL   =   " + "\n" + "$DISPLACEMENTS " + "\n" + "$REAL OUTPUT " + "\n")
                file.write("$SUBCASE ID =       Event " + "\n")
                file.write("$POINT ID =      " + str(accelerometer_ID[k]) + "  IDENTIFIED BY TIME " + "\n")
                for i in range(0, len(Time_step_matrix)):
                    L = "    " + str("%.6E" % float(Time_step_matrix[i])) + " G      " + str(
                        "%.6E" % Transformed_disp_matrix_Storage[k][i][0]) + "      " + str(
                        "%.6E" % Transformed_disp_matrix_Storage[k][i][1]) + "      " + str(
                        "%.6E" % Transformed_disp_matrix_Storage[k][i][2]) + "\n"
                    newstr = L.replace(" -", "-")
                    file.write(newstr)
                    file.write("-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00" + "\n")

        # print 'Progress bar: [*************]'
        prbarcounter += 1
        pbar.update(prbarcounter)  # adds a symbol at each iteration
        print 'Filed saved at:', filesave_path

        pbar.finish()
        print
    else:
        msgbox(msg="Neither .pch or .csv was chosen. \n Coordinate Transformation terminated.",
               title=title)
        sys.exit(0)



if __name__ == "__main__":

    cls()
    while ynbox('Perform coordinate transformation?', 'Coordinate transformation', ('Yes', 'No')):
        main()

    
    
    
    
    