# -*- coding: utf-8 -*-
"""

@author: Alexander Jorud, Carl Jacobsson
"""

"""

*** Improve visualization

*** Improve reading from pch method.

"""


import numpy as np
from easygui import *
import plotly
import plotly.graph_objs as go
from progressbar import *

# Generates accelerometer grid points.
def generate_grid_points():

    accelerometers = [201001,202001,207001,208001,209001,210001,215001,216001,219001,220001,221001,222001,223001,224001]
    amount_of_accelerometers = len(accelerometers)


    return {'accelerometers': accelerometers, 'amount_of_accelerometers': amount_of_accelerometers}


# Reads a file including the grid points which defines a coordinate system from Ansa.
def Open_Local_coord_sys():

    title = "Diagonal_method"
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
           title = title)
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

    msg_Globalcoord = "Choose visualization option. For better performance open the .html files in Google Chrome"
    choices_vis = ["Visualize imported grid points from the "+ filename +" file",
                   "Continue without any visualization."]
    choice_vis = choicebox(msg_Globalcoord, title, choices_vis)


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


# Method which reads a pch file and extracts the displacements and timesteps which
# can be used during the main method. Note only translational displacements are extracted.
# The rotational displacements can be extracted with minor alterations to the code.
def read_pch(filename):
    pch_file = open(filename, 'r')

    fo = pch_file.readlines()
    pch_file.close()

    Timestep_list = []
    X_tr_list = []
    Y_tr_list = []
    Z_tr_list = []
    X_rot_list = []
    Y_rot_list = []
    Z_rot_list = []

    title = "Diagonal method"
    msgbox(msg="Reading all E-points from pch...",title=title)
    choices_Globalcoord = ["201001", "202001", "207001", "208001", "209001", "210001",
                           "215001", "216001", "219001", "220001", "221001", "222001",
                           "223001", "224001"]
    Disp_storage_matrix = []
    widgets = ['Progressbar: ', Percentage(), ' ', Bar(marker='#', left='[', right=']'),
               ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=len(choices_Globalcoord))
    pbar.start()



    for pp in range(0,len(choices_Globalcoord)):
        pbar.update(pp)  # Adds a symbol at each iteration
        choice_Globalcoord = str(int(choices_Globalcoord[pp]) + 1000)
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

                if str_part[0] == '-CONT-':
                    X_rot_list.append(str_part[1])
                    Y_rot_list.append(str_part[2])
                    Z_rot_list.append(str_part[3])
                else:
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
                Timestep_list = []
                X_tr_list = []
                Y_tr_list = []
                Z_tr_list = []
                X_rot_list = []
                Y_rot_list = []
                Z_rot_list = []
        else:
            # No more lines to be read from file
            if choice_Globalcoord == str(int(choices_Globalcoord[-1])+1000):

                Timestep_matrix = np.array([Timestep_list]).astype(np.float32)

                Translation_matrix = np.matrix([np.asarray(X_tr_list).astype(np.float32),
                                                np.asarray(Y_tr_list).astype(np.float32),
                                                np.asarray(Z_tr_list).astype(np.float32)])
        Disp_storage_matrix.append(Translation_matrix)

    pbar.finish()
    print

    return {'Tra_mat': Disp_storage_matrix,
            'Time_mat': Timestep_matrix}


def main():

    title = "Diagonal Method"
    msgbox(msg="Choose pch-file.", title=title)
    filename1 = fileopenbox()
    print 'Open file:', filename1

    File_choice = read_pch(filename1)
    Disp_storage_matrix = File_choice['Tra_mat']
    Time_step_matrix = File_choice['Time_mat']
    CORD1R = Open_Local_coord_sys()['CORD1R_mat']
    accelerometer_ID = generate_grid_points()['accelerometers']
    nbr_of_accelerometers = generate_grid_points()['amount_of_accelerometers']
    Time_step_matrix = Time_step_matrix[0]


    Diag_storage_list = []

    cc = -1
    for acc in range(0,nbr_of_accelerometers):
        cc += 1
        if cc % 2 == 0:
            Disp_matrix1 = Disp_storage_matrix[acc]
            Disp_matrix2 = Disp_storage_matrix[acc+1]

            for k in range(0,len(CORD1R[0])):
                if CORD1R[0][k] == str(accelerometer_ID[acc]):
                    G1 = np.asarray(CORD1R[1][k]).astype(np.float32)
                    G2 = np.asarray(CORD1R[2][k]).astype(np.float32)
                elif CORD1R[0][k] == str(accelerometer_ID[acc]-1000):
                    G1 = np.asarray(CORD1R[1][k]).astype(np.float32)
                    G2 = np.asarray(CORD1R[2][k]).astype(np.float32)
            Diag_list = []

            for z in range(0,Disp_matrix1[0].size):

                G1_x = Disp_matrix1[0,z] + G1[0]
                G1_y = Disp_matrix1[1,z] + G1[1]
                G1_z = Disp_matrix1[2,z] + G1[2]

                G2_x = Disp_matrix2[0,z] + G2[0]
                G2_y = Disp_matrix2[1,z] + G2[1]
                G2_z = Disp_matrix2[2,z] + G2[2]


                # Compute the length between two points in space
                D = np.sqrt(((G1_x-G2_x)**2)+((G1_y-G2_y)**2)+((G1_z-G2_z)**2))

                # Store each diagonal to a list.
                Diag_list.append(D)
            Diag_plot = []

            for t in range(0,len(Diag_list)):

                if t == 0:

                    Diag_plot.append(0)

                else:
                    Diag_plot.append(Diag_list[t]-Diag_list[0])



            # Plot Dn-D0
            trace = go.Scatter(
                x=Time_step_matrix,
                y=Diag_plot,
                mode='lines',
                name="E-point: " + str(accelerometer_ID[acc]),
            )

            data = [trace]
            layout = go.Layout(
                title="E-point: " + str(accelerometer_ID[acc]),
                xaxis=dict(
                    title='Time [s]',
                    titlefont=dict(
                        family='Courier New, monospace',
                        size=26,
                        color='#7f7f7f'
                    )
                ),
                yaxis=dict(
                    title='[mm]',
                    titlefont=dict(
                        family='Courier New, monospace',
                        size=26,
                        color='#7f7f7f'
                    )
                )
            )
            fig = go.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename='Diagonal_method_plot_E-Point_' +
                                              str(accelerometer_ID[acc]) + '.html')

            Diag_storage_list.append(Diag_plot)



    # Save to file
    Diag_cc = -1
    counter = -1
    filesave_path = filesavebox()
    with open(filesave_path, 'w') as file:
        for k in range(0, len(accelerometer_ID)):
            file.write(
                "$TITLE   = Test_data" + "\n" + "$SUBTITLE= " + "\n" + "$LABEL   =   " + "\n" + "$DISPLACEMENTS " + "\n" + "$REAL OUTPUT " + "\n")
            file.write("$SUBCASE ID =         event " + "\n")
            file.write("$POINT ID =      " + str(accelerometer_ID[k]) + "  IDENTIFIED BY TIME " + "\n")
            counter += 1
            if counter % 2 == 0:
                Diag_cc += 1
                for i in range(0, len(Time_step_matrix)):
                    L = "    " + str("%.6E" % float(Time_step_matrix[i])) + " G      " + str("%.6E" % float(Diag_storage_list[Diag_cc][i])) +\
                        "      " + str("%.6E" % 0) + "      " + str("%.6E" % 0) + "\n"
                    newstr = L.replace(" -", "-")
                    file.write(newstr)
                    file.write("-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00" + "\n")
            else:
                for i in range(0, len(Time_step_matrix)):
                    L = "    " + str("%.6E" % float(Time_step_matrix[i])) + " G      " + str("%.6E" % 0) +\
                        "      " + str("%.6E" % 0) + "      " + str("%.6E" % 0) + "\n"
                    newstr = L.replace(" -", "-")
                    file.write(newstr)
                    file.write("-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00" + "\n")

    print 'Filed saved at:', filesave_path


if __name__ == "__main__":


    main()

    
    
    
    
    