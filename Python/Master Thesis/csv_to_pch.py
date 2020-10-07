# -*- coding: utf-8 -*-
"""

@author: Alexander Jorud, Carl Jacobsson
"""

"""
-----------------------IMPROVEMENTS-----------------------
Parts of the code is not necessary to perform a
transformation from csv file format to 
pch file format. Improvements should be carried out
to remove these parts of the code.

"""

import numpy as np
import os
import csv
from easygui import *
import sys


# Method to clear the print window after every run
def cls():
    os.system('CLS')

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

    reader_time = csv.reader(open(filename, "rb"), delimiter=",")
    tmp = list(reader_time)
    timestep_list = []

    for p in range(8,len(tmp)):
        timestep_list.append(tmp[p][1])

    return {'disp_mat': np.asarray(data).astype(np.float32), 'timestep_mat': timestep_list}

# Generates accelerometer grid points.
def generate_grid_points():

    accelerometers = [201001,202001,207001,208001,209001,210001,215001,216001,219001,220001,221001,222001,223001,224001]
    amount_of_accelerometers = len(accelerometers)


    return {'accelerometers': accelerometers, 'amount_of_accelerometers': amount_of_accelerometers}


def main():

    title = "csv to pch"
    msgbox(msg="Choose a csv-file.", title=title)
    filename = fileopenbox()
    print 'Open file:', filename

    #     READING CSV FILE SECTION
    if filename.endswith('.csv'):

        msgbox(msg=".csv file chosen.",
               title=title)

        File_choice = readCSV(filename)['disp_mat']
        Time_step_matrix = readCSV(filename)['timestep_mat']
        accelerometer_ID = generate_grid_points()['accelerometers']
        nbr_of_accelerometers = generate_grid_points()['amount_of_accelerometers']


        Transformed_disp_matrix_Storage = []

        # Loop over all accelerometers
        for acc in range(0, nbr_of_accelerometers):
            print 'Reading accelerometer ID:', accelerometer_ID[acc]

            Disp_matrix_tmp = File_choice[acc, :]

            Disp_matrix = np.matrix([Disp_matrix_tmp[:, 0],
                                     Disp_matrix_tmp[:, 1],
                                     Disp_matrix_tmp[:, 2]])


            Transformed_disp_matrix = []
            for i in range(0, Disp_matrix[0].size):
                Transformed_disp_matrix.append((Disp_matrix[:, i]))

            Transformed_disp_matrix_Storage.append(Transformed_disp_matrix)


        filesave_path = filesavebox()
        with open(filesave_path, 'w') as file:
            for k in range(0,nbr_of_accelerometers):
                file.write("$TITLE   = Test_data" + "\n" + "$SUBTITLE= " + "\n" + "$LABEL   =   " + "\n" + "$DISPLACEMENTS " + "\n" + "$REAL OUTPUT " + "\n")
                file.write("$SUBCASE ID =         event " + "\n")
                file.write("$POINT ID =      " + str(accelerometer_ID[k]) + "  IDENTIFIED BY TIME " + "\n")
                for i in range(0, len(Time_step_matrix)):
                    L = "    " + str("%.6E" % float(Time_step_matrix[i])) + " G      " + str(
                        "%.6E" % Transformed_disp_matrix_Storage[k][i][0]) + "      " + str(
                        "%.6E" % Transformed_disp_matrix_Storage[k][i][1]) + "      " + str(
                        "%.6E" % Transformed_disp_matrix_Storage[k][i][2]) + "\n"
                    newstr = L.replace(" -", "-")
                    file.write(newstr)
                    file.write("-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00" + "\n")


        print 'Filed saved at:', filesave_path



    else:
        msgbox(msg="A .csv file was not chosen. \n Program terminated.",
               title=title)
        sys.exit(0)


    
if __name__ == "__main__":

    cls()
    main()

    
    
    
    
    