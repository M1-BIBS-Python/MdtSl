#!/usr/bin/env python
# -*- coding : utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

def distanceCarree(p1, p2):
    """
    Calcule le carre de la distance entre 2 points dans l'espace.
    :param p1: Premier point.
    :param p2: Second point.
    :return: Le carre de la distance (nombre reel) entre les 2 points.
    """
    return ((p1['x'] - p2['x'])**2 + (p1['y'] - p2['y'])**2 + (p1['z'] - p2['z'])**2)

def centerOfMass(dico):
    """
    Calcule le centre de masse d'une molecule.
    :param dico: Dictionnaire contenant les coordonnees des atomes de la molecule.
    :return: Dictionnaire contenant les coordonnees du centre de masse de la molecule.
    """
    x = 0
    y = 0
    z = 0
    nbAtomes = 0
    CM_res = {}
    for atom in dico['atomlist']:
        x += dico[atom]['x']
        y += dico[atom]['y']
        z += dico[atom]['z']
        nbAtomes += 1
    CM_res['x'] = x / nbAtomes
    CM_res['y'] = y / nbAtomes
    CM_res['z'] = z / nbAtomes

    return CM_res


#-------------------------------------------
def RMSD_prot(dico1, dico2, mode):
    """
    Calcule le RMSD entre 2 proteines entieres (par exemple des structures superposees).
    :param dico1: Dictionnaire contenant les donnees parsees depuis le fichier pdb de la premiere proteine.
    :param dico2: Dictionnaire contenant les donnees parsees depuis le fichier pdb de la seconde proteine.
    :param mode: Mode de calcul du rmsd (par rapport aux carbones alphe, 'CA', ou au centre de masse de la molecule, 'CM').
    :return: La valeur du rmsd entre les 2 structures.
    """

    nb_pairs = 0
    somme = 0
    for chain in dico1.keys():
        for res in dico1[chain]['reslist']:					# pour chaque residu de la proteine
            if mode == 'CA':								# calcul du rmsd par rapport au carbone alpha du residu
                for atom in dico1[chain][res].keys():
                    if atom == mode:
                        d = distanceCarree(dico1[chain][res][atom], dico2[chain][res][atom])
                        somme += d
                        nb_pairs += 1
            elif mode == 'CM':								# calcul du rmsd par rapport au centre de masse du residu
                CM_res1 = centerOfMass(dico1[chain][res])	# calcul des coordonnes du centre de masse du residu de la premiere structure
                CM_res2 = centerOfMass(dico2[chain][res])	# idem pour le residu correspondant dans la seconde structure
                d = distanceCarree(CM_res1, CM_res2)
                somme += d
                nb_pairs += 1

    rmsd = sqrt(somme / nb_pairs)
    return rmsd



def RMSD_domain(dico1, dico2, mode):
    """
    Calcule le RMSD entre 2 domaines de proteine.
    :param dico1: Dictionnaire contenant les donnees parsees pour le domaine de la proteine 1.
    :param dico2: Dictionnaire contenant les donnees parsees pour le domaine de la proteine 2.
    :param mode: Nom des atomes a partir desquels le RMSD sera calcule.
    :return: La valeur du rmsd entre les 2 domaines.
    """
    nb_pairs = 0
    somme = 0
    for res in dico1['reslist']:
        if mode == 'CA':
            for atom in dico1[res]['atomlist']:
                if atom == mode:
                    d = distanceCarree(dico1[res][atom], dico2[res][atom])
                    somme += d
                    nb_pairs += 1
        elif mode == 'CM':
            CM_res1 = centerOfMass(dico1[res])
            CM_res2 = centerOfMass(dico2[res])
            d = distanceCarree(CM_res1, CM_res2)
            somme += d
            nb_pairs += 1
    rmsd = sqrt(somme / nb_pairs)
    return rmsd


# -------------------------------------------------------------
def computeRMSD(ref, frames, list_dom_prot, rmsd_mode, output):
    """
    Calcule le RMSD entre la structure de reference et chacune des conformations de la dynamique.
    :param ref: Dictionnaire correspondant a la structure de reference.
    :param frames: Dictionnaire correspondant aux differentes conformations.
    :param list_dom_prot: Liste des domaines proteiques.
    :param rmsd_mode: Mode de calcul du RMSD.
    :param output: Fichier de sortie contenant pour chaque conformation, le RMSD global et celui des domaines.
    :return: Graphes RMSD en fonction de conformations.
    """
    f = open(output, "w")

    x_plot = []
    y_global = []
    y_dom = dict()  # cle = nom du domaine, valeur = liste des RMSD du domaine

    for model in sorted(map(int, frames.keys())):   # pour chaque conformation on calcule le RMSD global et le RMSD de chaque domaine

        rmsd_global = RMSD_prot(ref, frames[str(model)], rmsd_mode)

        x_plot.append(model)
        y_global.append(rmsd_global)
        f.write("Model " + str(model) + "\t" + str(rmsd_global) + "\n")

        for dom in list_dom_prot:
            if dom not in y_dom.keys():
                y_dom[dom] = []

            rmsd_dom = RMSD_domain(ref[dom], frames[str(model)][dom], rmsd_mode)

            y_dom[dom].append(rmsd_dom)

            f.write("\t" + str(dom) + "\t" + str(rmsd_dom) + "\n")
    f.close()

    # Graphes:
    plt.plot(x_plot, y_global)
    plt.xlabel('Frame')
    plt.ylabel('RMSD')
    plt.title('Global RMSD of the protein')
    plt.show()

    # tracer les courbes des domaines sur le meme graphe:
    colors = ['b', 'r', 'g', 'y', 'c', 'm', 'k']
    for i in range(len(list_dom_prot)):
        plt.plot(x_plot, y_dom[list_dom_prot[i]], colors[i], label=list_dom_prot[i])

    plt.xlim(0, 5000)
    plt.xlabel('Frame')
    plt.ylabel('RMSD')
    plt.title('RMSD for each domain')
    plt.legend(loc=4)
    plt.show()