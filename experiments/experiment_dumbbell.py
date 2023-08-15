from utilities.graphs import ConnectionGraph
from utilities.connection_resistance import ConnectionResistance
import numpy as np
from pathlib import Path
import os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from utilities.plotstyle_maps import get_color_map, get_dash_map, get_marker_map

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from definitions import CONFIG_MATPLOTSTYLE_PATH
plt.style.use(CONFIG_MATPLOTSTYLE_PATH)

from definitions import CONFIG_DUMBBELL_GRAPH
from config.yaml_functions import yaml_loader

if __name__=='__main__':
    experiment_name = 'dumbbell_graph'
    indices = [(0,1), (0,2), (0,3), (0,7)]
    indices_markfan_exp = [(0,1), (1, 2), (3, 4), (6, 7)]
    conf_cg = yaml_loader(CONFIG_DUMBBELL_GRAPH)
    tag = "_deg90"
    #tag = ""
    rotations = np.arange(0, 365, 5)
    rotations = list(rotations)
    if tag == "_deg90":
        rotations.remove(270)
        special_rotations = [270]
    else:
        rotations.remove(0)
        rotations.remove(360)
        special_rotations = [0, 360]
    rotations = np.array(rotations)
    rotations = np.sort(rotations)

    connectionGraph = ConnectionGraph(conf_cg.copy())
    connectionResistance = ConnectionResistance(conf_cg.copy(), connectionGraph)

    # Calculate standard effective resistance
    er_standard = connectionResistance.calc_effective_resistance()

    cr_chung2014 = []
    cr_trace = []
    for rotation in rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        edge_rotations[0][0] = rotation

        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        # Calculate chung 2014 connection resistance
        cr_chung2014.append(connectionResistance.calc_connection_resistance_chung2014())

        # Calculate trace connection resistance
        cr_trace.append(connectionResistance.calc_connection_resistance())

    cr_chung2014_special_rot = []
    for rotation in special_rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        edge_rotations[0][0] = rotation

        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        cr_chung2014_special_rot.append(connectionResistance.calc_connection_resistance_chung2014())


    cr_trace = np.array(cr_trace)
    cr_chung2014 = np.array(cr_chung2014)
    cr_chung2014_special_rot = np.array(cr_chung2014_special_rot)

    # Formatting
    color_map = get_color_map()
    dash_map = get_dash_map()
    marker_map = get_marker_map()

    n = er_standard.shape[0]
    #xticks = np.arange(0, 380, 120)
    xticks = [0, 90, 180, 270, 360]
    xtickslabels = [r'$0$', r'$\pi/2$', r'$\pi$', r'$2\pi/3$', r'$2\pi$']
    yticks = [-3, 0, 3]
    edges = connectionResistance.get_edges()

    fig, axes = plt.subplots(1, len(indices), figsize=(17, 3))
    for k, (i,j) in enumerate(indices_markfan_exp):
        axes[k].plot([0, rotations[-1]], [er_standard[i, j], er_standard[i, j]],
                     color=color_map['color_darkgray'],
                     label='ER')
        if (i,j) in edges:
            axes[k].plot(rotations, cr_chung2014[:, i, j],
                        dashes=dash_map['dash3'],
                        color=color_map['color_cyan'],
                        label='CR Chung et al (2014)')
            for p, rot in enumerate(special_rotations):
                axes[k].plot(rot, cr_chung2014_special_rot[p, i, j],
                             marker='x',
                             markersize='12',
                             color=color_map['color_cyan'])
        axes[k].plot(rotations, cr_trace[:, i, j],
                     dashes=dash_map['dash3'],
                     color=color_map['color_darkgreen'],
                     label='CR Def. 6.7')
        axes[k].set_xlim([0, 360])
        axes[k].set_ylim([0.5, 1.5])
        axes[k].set_xticks(xticks)
        axes[k].set_xticklabels(xtickslabels)

        if k==0:
            axes[k].legend()
            axes[k].set_ylabel('Effective resistance')
        axes[k].set_title(f'Vertices {(i+1, j+1)}')
        axes[k].set_xlabel('Rotation')
    plt.tight_layout()
    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'dumbbell_graph_comparisons_markfan'+ tag + '.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'dumbbell_graph_comparisons_markfan'+ tag + '.png'), format='png')
