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
    if mode == "CM":  # distance entre les centres de masse des residus
        CM_res1 = centerOfMass(res1)
        CM_res2 = centerOfMass(res2)
        minval = math.sqrt(distanceCarree(CM_res1, CM_res2))

    elif mode == "atom":  # distance entre les atomes les plus proches des 2 residus
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
    :return: La matrice des distances.
    """
    list1 = dico[dom1]['reslist']  # liste des residus du premier domaine
    list2 = dico[dom2]['reslist']  # liste des residus du second domaine
    nb_row = len(list1)  # nombre de lignes de la matrice = nombre de residus dans le premier domaine
    nb_col = len(list2)  # nombre de colonnes de la matrice = nombre de residus dans le second domaine
    distmat = []
    dic1 = dico[dom1]  # pour simplifier l'ecriture par la suite
    dic2 = dico[dom2]

    for i in range(nb_row):  # pour chaque residu  i du premier domaine ...
        listi = []  # chaque ligne de la matrice est une liste
        for j in range(nb_col):  # pour chaque residu j du second domaine
            dij = distDico(dic1[list1[i]], dic2[list2[j]],
                           mode)  # ... on calcule sa distance avec le j-eme residu du second domaine
            listi.append(dij)  # et on ajoute cette distance a la liste representant la i-eme ligne de la matrice
        distmat.append(listi)  # la ligne est ajoutee a la matrice : distmat est une liste de liste

    # Optionnel: representer la matrice sous forme de heatmap
    if 'heatmap' in keyword_param:
        if keyword_param['heatmap'] == "yes":  # si l'appel de la fonction comprend l'argument heatmap="yes"
            data = np.array(distmat).reshape(nb_row, nb_col)  # la matrice est transformee en array
            plt.pcolor(data, cmap="RdBu_r")  # choix du systeme de coloration
            plt.xticks(np.arange(nb_col) + 0.5, list2, fontsize=5)  # positionnement des etiquettes sur l'abscisse
            plt.yticks(np.arange(nb_row) + 0.5, list1, fontsize=5)  # idem sur l'ordonnee
            plt.xlabel(dom2)  # le domaine 2 est represente par les colonnes (abscisse) de la matrice
            plt.ylabel(dom1)
            cb = plt.colorbar()  # ajout d'une echelle des couleurs
            cb.set_label('Distance between residues (in Angstrom)')  # titre de la legende
            plt.show()

    return distmat


def resInterface(dico, prot_domains, rna_dom, threshold, mode, **keyword_param):
    """
    Calcule les frequences d'appartenance a l'interface proteine/ARN pour chacun des residus de la proteine, et les retourne dans un fichier texte de sortie.
    :param dico: Dictionnaire representant le complexe proteine/ARN.
    :param prot_domains: Liste contenant les identifiants des domaines proteiques a considerer.
    :param rna_dom: Identifiant du domaine correspondant a l'ARN.
    :param threshold: Seuil en Angstrom, correspond a la distance ARN-residu en-dessous de laquelle le residu appartient a l'interface.
    :param mode: Mode de calcul de la distance entre les domaines.
    :param keyword_param: Parametre optionnel: si present, le dictionnaire d'entree est modifie en ajoutant pour chaque residu le B-factor (0.00 si le residu n'appartient pas a l'interface, 1.00 sinon).
    :return: Dictionnaire correspondant aux residus dont la frequence d'appartenance a l'interface est non nulle.
    """
    inInterface = dict()
    for dom in prot_domains:  # pour chaque domaine proteique
        inInterface[dom] = dict()  # on cree un dictionnaire de dictionnaire
        inInterface[dom]['reslist'] = []  # contiendra la liste des residus du domaine
        for model in dico.keys():  # pour chaque frame
            for res1 in dico[model][dom]['reslist']:  # pour chaque residu
                if res1 not in inInterface[dom]['reslist']:
                    inInterface[dom][
                        res1] = 0  # on initialise a 0 le nombre de fois que le residu se trouve dans l'interface
                    inInterface[dom]['reslist'].append(res1)  # on ajoute le residu dans la liste
                min = 60  # initialisation de la distance minimale residu-ARN
                for res2 in dico[model][rna_dom]['reslist']:
                    d = distDico(dico[model][dom][res1], dico[model][rna_dom][res2],
                                 mode)  # on calcule la distance entre le residu et tous les nucleotides de l'ARN
                    dico[model][rna_dom][res2]['bfactor'] = 0.00
                    if d < min:  # si la distance entre le residu et le nucleotide est plus petite que la valeur minimale
                        min = d  # min sera donc la distance entre le residu et le nucleotide le plus proche
                if min <= threshold:
                    inInterface[dom][res1] = inInterface[dom][res1] + 1
                    dico[model][dom][res1]['bfactor'] = 1.00
                else:
                    dico[model][dom][res1]['bfactor'] = 0.00

    # Optionnel: remet le dictionnaire modifie au format pdb, pour une visualisation dans PyMol
    if 'writePDB' in keyword_param:
        writePDBframes(dico, keyword_param['writePDB'])

    # Retourne les frequences (si non nulles) dans un fichier texte:
    f = open("interfaceFreq.txt", "w")
    res_interface = dict()
    for dom in inInterface.keys():
        res_interface[dom] = dict()
        f.write("Domain " + dom + "\n")
        for res in inInterface[dom]['reslist']:
            freq = (inInterface[dom][res] / len(dico.keys())) * 100
            if freq != 0:
                f.write("\tResidue " + res + " : " + str(freq) + " %\n")
                res_interface[dom][res] = freq
    f.close()

    return res_interface


def writePDBframes(dico, output):
    """
    Utilise un dictionnaire contenant plusieurs conformations pour ecrire un fichier pdb visualisable sous PyMol.
    :param dico: Dictionnaire contenant les differentes conformations.
    :param output: Nom du fichier de sortie.
    :return: Un fichier pdb au format ATOM.
    """
    fout = open(output, "w")

    for model in dico.keys():
        fout.write("MODEL\t" + str(model) + "\n")
        for dom in ['A1', 'A2', 'A3', 'A4']:
            for res in dico[model][dom]['reslist']:
                for atom in dico[model][dom][res]['atomlist']:
                    fout.write(
                        "{:6s}{:5s} {:^4s}{:1s}{:3s} {:1s}{:4s}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}       {:^4s}\n".format(
                            "ATOM", dico[model][dom][res][atom]['id'], atom, '', dico[model][dom][res]['resname'],
                            '', res, '', dico[model][dom][res][atom]['x'], dico[model][dom][res][atom]['y'],
                            dico[model][dom][res][atom]['z'],
                            1.00, dico[model][dom][res]['bfactor'], dom))

        for nucl in dico[model]['B']['reslist']:
            for atom in dico[model]['B'][nucl]['atomlist']:
                fout.write(
                    "{:6s}{:5s} {:^4s}{:1s}{:3s} {:1s}{:4s}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}       {:^4s}\n".format(
                        "ATOM", dico[model]['B'][nucl][atom]['id'], atom, '', dico[model]['B'][nucl]['resname'],
                        '', nucl, '', dico[model]['B'][nucl][atom]['x'], dico[model]['B'][nucl][atom]['y'],
                        dico[model]['B'][nucl][atom]['z'],
                        1.00, 0.00, 'B'))


def contactTime(pairs, frames_dico, threshold, duration, mode):
    """
    Calcule le temps de contact entre differentes paires de residus.
    :param pairs: 
    :param pairs: Dictionnaire contenant les paires de residus.
    :param frames_dico: Dictionnaire contenant toutes les conformations.
    :param threshold: Seuil (en Angstrom) pour definir le contact.
    :param duration: Duree de la dynamique.
    :param mode: Mode de calcul des distances.
    :return: Le dictionnaire pairs auquel ont ete ajoutes les temps de contact.
    """
    for res in pairs.keys():
        pairs[res]['contact time'] = 0
        for model in frames_dico.keys():
            d = distDico(frames_dico[model][pairs[res]['dom']][res], frames_dico[model][pairs[res]['rna']][pairs[res]['nucl']], mode)
            if d <= threshold:
                pairs[res]['contact time'] += 1

        pairs[res]['contact time'] = pairs[res]['contact time']*duration/len(frames_dico.keys())

    return pairs


