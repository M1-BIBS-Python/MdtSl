#!/usr/bin/env python
# -*- coding : utf8 -*-

from RMSD import *

import numpy as np
import matplotlib.pyplot as plt
import math

def distDico(res1, res2, mode):
    """
    Calcule la distance entre 2 residus.
    :param res1: Dictionnaire correspondant au premier residu.
    :param res2: Dictionnaire correspondant au second residu.
    :param mode: Mode de calcul de la distance (distance entre les 2 atomes les plus proches, ou entre les centres de masse).
    :return: La distance entre les 2 residus.
    """
    if mode == "CM":
        CM_res1 = centerOfMass(res1)
        CM_res2 = centerOfMass(res2)
        minval = math.sqrt(distanceCarree(CM_res1, CM_res2))

    elif mode == "atom":
        minval = 1000000
        for atom1 in res1['atomlist']:
            for atom2 in res2['atomlist']:
                dist = math.sqrt(distanceCarree(res1[atom1], res2[atom2]))
                if dist < minval:
                    minval = dist

    return minval


def distMatrix(dico, dom1, dom2, mode, **keyword_param):
    """
    Calcule la matrice des distances entre 2 domaines proteiques ou entre 1 domaine et l'ARN.
    :param dico: Dictionnaire correspondant au complexe proteine-ARN.
    :param dom1: Nom du premier domaine a utiliser.
    :param dom2: Nom du second domaine a utiliser.
    :param mode: Mode de calcul de la distance entre les domaines: par rapport au centre de masse ('CM') ou entre les 2 atomes les plus proches ('atom').
    :param keyword_param: Parametre optionnel pour la representation graphique de la matrice.
    :return: 
    """
    list1 = dico[dom1]['reslist']
    list2 = dico[dom2]['reslist']
    nb_row = len(list1)
    nb_col = len(list2)
    distmat = []
    dic1 = dico[dom1]
    dic2 = dico[dom2]

    for i in range(nb_row):
        listi = []
        for j in range(nb_col):
            dij = distDico(dic1[list1[i]], dic2[list2[j]], mode)
            listi.append(dij)
        distmat.append(listi)

    if 'heatmap' in keyword_param:
        if keyword_param['heatmap'] == "yes":   # si l'appel de la fonction comprend l'argument heatmap="yes"
            data = np.array(distmat).reshape(nb_row, nb_col)
            plt.pcolor(data, cmap="RdBu_r")
            plt.xticks(np.arange(nb_col) + 0.5, list2, fontsize=5)
            plt.yticks(np.arange(nb_row) + 0.5, list1, fontsize=5)
            plt.xlabel(dom2)
            plt.ylabel(dom1)
            cb = plt.colorbar()
            cb.set_label('Distance between residues (in Angstrom)')
            plt.show()

    return distmat





