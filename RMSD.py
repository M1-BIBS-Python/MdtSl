#!/usr/bin/env python
# -*- coding : utf8 -*-

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
    for atom in dico.keys():
        if atom != 'resname' and atom != 'atomlist':
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
    :return: La valeur du rmsd (reel).
    """

    nb_pairs = 0
    somme = 0
    for chain in dico1.keys():
        for res in dico1[chain].keys():
            if res != 'reslist':
                if mode == 'CA':
                    for atom in dico1[chain][res].keys():
                        if atom == mode:
                            d = distanceCarree(dico1[chain][res][atom], dico2[chain][res][atom])
                            somme += d
                            nb_pairs += 1
                elif mode == 'CM':
                    CM_res1 = centerOfMass(dico1[chain][res])
                    CM_res2 = centerOfMass(dico2[chain][res])
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
    :param calc_atom: Nom des atomes a partir desquels le RMSD sera calcule.
    :return: La valeur du rmsd (reel).
    """
    nb_pairs = 0
    somme = 0
    for res in dico1.keys():
        if res != 'reslist':
            if mode == 'CA':
                for atom in dico1[res].keys():
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

