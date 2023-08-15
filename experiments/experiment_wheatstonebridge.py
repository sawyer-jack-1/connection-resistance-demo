from utilities.graphs import ConnectionGraph
from utilities.connection_resistance import ConnectionResistance
import numpy as np
from pathlib import Path
import os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from definitions import CONFIG_MATPLOTSTYLE_PATH
plt.style.use(CONFIG_MATPLOTSTYLE_PATH)

from utilities.plotstyle_maps import get_color_map, get_dash_map, get_marker_map
from definitions import CONFIG_WHEATSTONEBRIDGE_GRAPH
from config.yaml_functions import yaml_loader

if __name__=='__main__':
    experiment_name = 'wheatstone_bridge_graph'
    indices = [(0, 1), (1, 2), (2, 3), (3, 0), (1, 3), (0, 2)]
    indices_markfan_exp = [(1, 0), (1, 2), (1, 3), (0, 2)]
    conf_cg = yaml_loader(CONFIG_WHEATSTONEBRIDGE_GRAPH)

    rotations = np.arange(-180, 180, 5)
    rotations = list(rotations)
    rotations.remove(0)
    rotations = np.array(rotations)
    rotations = np.sort(rotations)
    special_rotations = [0]

    connectionGraph = ConnectionGraph(conf_cg.copy())
    connectionResistance = ConnectionResistance(conf_cg.copy(), connectionGraph)

    # Calculate standard effective resistance
    er_standard = connectionResistance.calc_effective_resistance()

    cr_chung2014 = []
    cr_trace = []
    for rotation in rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        edge_rotations[2][0] = rotation

        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        # Calculate mark and fan connection resistance
        cr_chung2014.append(connectionResistance.calc_connection_resistance_chung2014())

        # Calculate trace connection resistance
        cr_trace.append(connectionResistance.calc_connection_resistance())

    cr_chung2014_special_rot = []
    for rotation in special_rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        edge_rotations[2][0] = rotation

        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        cr_chung2014_special_rot.append(connectionResistance.calc_connection_resistance_chung2014())

    cr_chung2014 = np.array(cr_chung2014)
    cr_chung2014_special_rot = np.array(cr_chung2014_special_rot)
    cr_trace = np.array(cr_trace)

    ### Make plots
    color_map = get_color_map()
    dash_map = get_dash_map()
    marker_map = get_marker_map()

    n = er_standard.shape[0]
    xticks = [-180, -90, 0, 90, 180]
    xtickslabels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
    edges = connectionResistance.get_edges()

    # Figure for main
    fig, axes = plt.subplots(2, 3, figsize=(17, 6))
    l=0
    for k, (i,j) in enumerate(indices):
        if k > 2:
            l = 1
        axes[l, k % 3].plot([-180, rotations[-1]], [er_standard[i, j], er_standard[i, j]],
                     color=color_map['color_darkgray'],
                     label='ER')
        axes[l, k % 3].plot(rotations, cr_chung2014[:, i, j],
                         dashes=dash_map['dash3'],
                         color=color_map['color_cyan'],
                         label='CR Chung et al (2014)')
        axes[l, k % 3].plot(special_rotations[0], cr_chung2014_special_rot[0, i, j],
                         marker='x',
                         markersize='12',
                         color=color_map['color_cyan'])
        axes[l, k % 3].plot(rotations, cr_trace[:, i, j],
                     dashes=dash_map['dash4'],
                     color=color_map['color_darkgreen'],
                     label='CR Def. 6.7')
        axes[l, k % 3].set_xlim([-180, 180])
        axes[l, k % 3].set_ylim([0.13, 0.38])
        axes[l, k % 3].set_xticks(xticks)
        axes[l, k % 3].set_xticklabels(xtickslabels)

        if k==0:
            axes[l, k % 3].legend()
            axes[l, k % 3].set_ylabel('Magnitude')
        axes[l, k % 3].set_title(f'Vertices {(i+1, j+1)}')
        axes[l, k % 3].set_xlabel('Rotations')
    plt.tight_layout()
    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'wheatstonebridge_graph_comparisons.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'wheatstonebridge_graph_comparisons.png'), format='png')


