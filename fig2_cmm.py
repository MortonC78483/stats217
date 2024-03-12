# -*- coding: utf-8 -*-
# author: Yu Fu, Louise A. Dennis
# Update version used for large network(>99/use Networkx)

from DTMCSN import Infection_Model, generate_FCN_graph
import sys
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

nodes_list = [10] # nodes in graph
for nodes in nodes_list:
    multiples = 30 # how many messages to run graphs for
    runs = 100 # how many times to run each graph
    prob1 = 0.5
    prob2 = 0.5

    seed_list = [608, 22399, 30799, 45619, 86914, 88199, 95479, 97124, 102228, 110381, 123191, 128066, 136198, 138622, 144029, 145665, 150140, 155098, 160442, 167664, 174155, 178254, 182472, 205539, 207874, 211275, 242730, 253936, 273295, 277832, 283560, 283621, 287753, 293983, 342667, 368661, 388498, 398682, 399457, 401114, 405783, 411091, 414681, 419865, 425440, 436392, 437673, 445677, 463650, 466866, 466979, 480885, 490791, 493598, 512922, 515091, 515872, 520311, 523216, 523537, 535565, 540857, 540884, 543581, 558932, 560267, 570716, 582510, 614492, 629770, 632541, 672255, 689666, 706116, 708028, 712723, 741788, 746230, 810533, 812966, 813561, 816066, 825575, 831753, 832583, 849195, 861630, 866382, 888709, 888810, 917493, 925522, 936628, 938122, 955047, 972166, 974518, 976014, 991753, 993612]

    pn_diff = []
    pn_sim = []
    for i in range(multiples):
        pn_diff.append([])
        pn_sim.append([])

    edges = 3
    res = np.zeros((len(seed_list), multiples))

    for i in np.arange(len(seed_list)):
        seed = seed_list[i]
        G = generate_FCN_graph(nodes, edges, seed)
        random.seed(seed)
        m = 0
        print("\nWorking on Graph " + str(i+1), end = " ")
        while m < multiples:
            experiment = 0
            # running Infection_Model outputs
            v  = Infection_Model(nodes, G, prob1, prob2, (m), runs, experiment)
            res[i,m] = v

            m = m + 1

    mean = np.mean(res, axis = 0)*nodes/100
    high = np.quantile(res, .9, axis = 0)*nodes/100
    low = np.quantile(res, .1, axis = 0)*nodes/100

    df_out = pd.DataFrame(np.vstack((mean, high, low)))
    df_out.index = ["mean", "high", "low"]
    df_out.columns +=1
    df_out.to_csv("results/cmm/cmm_results_fig2_.csv")


    plt.figure()
    plt.scatter(np.arange(len(mean))+1, mean, color = "blue")
    plt.fill_between(np.arange(len(mean))+1, low, high, alpha=0.2, color = "blue")
    plt.xlabel("Number of Messages")
    plt.ylabel("Expected Infection")
    plt.savefig("results/cmm/cmm_results_fig2.png")
