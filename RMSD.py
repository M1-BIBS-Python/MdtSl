#!/usr/bin/env python
# -*- coding : utf8 -*-

"""
Description: Functions that calculate the RMSD between 2 proteins or for a domain
"""

import ParserPDB

from math import sqrt
import re
import matplotlib.pyplot as plt
import numpy as np
import os


def Distance(p1, p2):
    """
    Calcule la distance entre 2 points dans l'espace.
    :param p1: Premier point.
    :param p2: Second point.
    :return: La distance (nombre reel) entre les 2 points.
    """
    return sqrt(
        (p1['x'] - p2['x']) * (p1['x'] - p2['x']) + (p1['y'] - p2['y']) * (p1['y'] - p2['y']) + (p1['z'] - p2['z']) * (p1['z'] - p2['z']))  # calcul du RMSD entre deux proteines entieres



def RMSD_prot(dico1, dico2, calc_atom):
    """
    Calcule le RMSD entre 2 proteines entieres (par exemple des structures superposees).
    :param dico1: Dictionnaire contenant les donnees parsees depuis le fichier pdb de la premiere proteine.
    :param dico2: Dictionnaire contenant les donnees parsees depuis le fichier pdb de la seconde proteine.
    :param calc_atom: Nom des atomes a partir desquels le RMSD sera calcule.
    :return: La valeur du rmsd (reel).
    """
    nb_pairs = 0
    somme = 0
    for chain in dico1.keys():
        for res in dico1[chain].keys():
            if res != 'reslist':
                for atom in dico1[chain][res].keys():
                    if re.match(calc_atom, atom):
                        d = Distance(dico1[chain][res][atom], dico2[chain][res][atom])
                        somme += d * d
                        nb_pairs += 1

    rmsd = sqrt(somme / nb_pairs)
    return rmsd



def RMSD_domain(dico1, dico2, calc_atom):
    """
    Calcule le RMSD entre 2 domaines de proteine.
    :param dico1: Dictionnaire contenant les donnees parsees pour le domaine de la proteine 1.
    :param dico2: Dictionnaire contenant les donnees parsees pour le domaine de la proteine 2.
    :param calc_atom: Nom des atomes a partir desquels le RMSD sera calcule.
    :return: La valeur du rmsd (reel).
    """
    for res in dico1.keys():
        nb_pairs = 0
        somme = 0
        if res != 'reslist':
            for atom in dico1[res].keys():
                if re.match(calc_atom, atom):
                    d = Distance(dico1[res][atom], dico2[res][atom])
                    somme += d * d
                    nb_pairs += 1
    rmsd = sqrt(somme / nb_pairs)
    return rmsd
                
