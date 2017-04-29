#!/usr/bin/env python
# -*- coding : utf8 -*-

from ParserPDB import *
from computeInterface import *
import sys, os

def usage():
    print ("""
    
    obligatory:
    ===========
    
    -ref    -> pdb file containing the reference structure of the protein/RNA complex
    
    -conf   -> pdb file containing the different conformations of the complex
    
    
    optional:
    =========
    
    -th     -> threshold to define a contact,in Ansgtrom (default = 9.0)
    
    -rmsd   -> if rmsd = 'CA', computes the RMSD between alpha carbons of two residues
               and returns the smallest distance.
               if rmsd = 'CM', commputes the RMSD between the centers of mass of the two
               residues and returns it. (default = 'CM')
                   
    -mode   -> if mode = 'atom', computes the distance between all the atoms of two residues
               or of a residue and a nucleotide and returns the smallest distance.
               if mode = 'CM', computes the distance between the centers of mass and returns it.
               (default = 'CM)
    """)


def exists_file(f):
    """
    Teste si un fichier de sortie existe deja dans le repertoire courant.
    :param f: Nom du fichier de sortie.
    :return: True si le fichier existe deja, False sinon
    """
    if os.path.exists(f):
        return True
    return False

def overwrite_file(f):
    """
    Si un fichier de sortie existe deja, demande a l'utilisateur s'il veut l'ecraser.
    :param f: Nom du fichier de sortie.
    :return: Nom du fichier de sortie.
    """
    overwrite = "No"
    while ((overwrite == "No" or overwrite == "no") and exists_file(f)):
        overwrite = input("This file already exists ! Do you want to overwrite this file ? (yes / no) ")
        if overwrite == "No" or overwrite == "no":
            f = input("Please, enter the name of the output file: ")
    return f


# Get arguments
# =============

try:
    ref_file = sys.argv[sys.argv.index("-ref")+1]
except:
    usage()
    print("ERROR: please, enter the name of the reference pdb input")
    sys.exit()

try:
    conf_file = sys.argv[sys.argv.index("-conf")+1]
except:
    usage()
    print("ERROR: please, enter the name of the conformations pdb input")
    sys.exit()

try:
    f = open(ref_file)
    f.close()
except:
    print("ERROR: this file does not exist: ", ref_file)
    sys.exit()

try:
    f = open(conf_file)
    f.close()
except:
    print("ERROR: this file does not exist: ", conf_file)
    sys.exit()

try:
    threshold = float(sys.argv[sys.argv.index("-th")+1])
    #print("The threshold used to define a contact is ", threshold)
except:
    threshold = 9.0

try:
    rmsd_mode = sys.argv[sys.argv.index("-rmsd")+1]
    #print("The mode of calculation for the RMSD is ", rmsd_mode)
except:
    rmsd_mode = "CM"

try:
    dist_mode = sys.argv[sys.argv.index("-mode")+1]
    #print("The mode of calculation for the distance is ", dist_mode)
except:
    dist_mode = "CM"


list_dom_prot = input("Please, enter the list of proteic domains identifiers (example: A1,A2,A3,A4) :").split(sep=",")
dom_rna = input("Please, enter the RNA domain identifier (example: B) :")
parsing_list = list(list_dom_prot)
parsing_list.append(dom_rna)

# ----------------------------------------
# Parsing de la conformation de reference
# ----------------------------------------

ref = PDBparser(ref_file, parsing_list)


# ------------------------------
# Parsing des 500 conformations
# ------------------------------

frames = PDBparserMulti(conf_file, parsing_list)


# --------------------------------------------------------------------------------------------
# Calcul du RMSD entre la structure de reference et chacune des conformations de la dynamique
# --------------------------------------------------------------------------------------------

rmsd_output = input("Please, enter the name of the output file to store the RMSDs\n(example: rmsd.txt): ")
# Verifie si le fichier existe deja:
rmsd_output = overwrite_file(rmsd_output)

computeRMSD(ref, frames, list_dom_prot, rmsd_mode, rmsd_output)


# ---------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------
# Calcul de la matrice des distances entre chaque domaine proteique et l'ARN pour la structure de reference
# ----------------------------------------------------------------------------------------------------------

distmat_ref = input("Do you want to show the distance matrices (between the proteic domains and RNA)\nfor the reference structure as heatmaps ? (yes / no) ")
if distmat_ref == "yes" or distmat_ref == "Yes":
    for dom in list_dom_prot:
        distMatrix(ref, dom, dom_rna, dist_mode)


# -------------------------------------------------------------------------------------------------
# Calcul de la frequence d'appartenance a l'interface avec l'ARN pour chaque residu de la proteine
# -------------------------------------------------------------------------------------------------

freq_output = input("Please, enter the name of the output file to store the frequences\nof belonging to the interface (example: freq.txt): ")
freq_output = overwrite_file(freq_output)

writing = input("Do you want to get a pdb file with B-factors different for\nresidues belonging to the interface ? (yes / no)")
if writing == "yes" or writing == "Yes":
    bfactor_output = input("Please enter the name of the output pdb file (example: bfactor.pdb):")
    bfactor_output = overwrite_file(bfactor_output)
    freq_interface = resInterface(frames, list_dom_prot, dom_rna, threshold, dist_mode, freq_output, writePDB=bfactor_output)
else:
    freq_interface = resInterface(frames, list_dom_prot, dom_rna, threshold, dist_mode, freq_output)


# ------------------------------------------------------------------------------------
# Calcul des temps de contact entre paires de residus choisis a partir de la Figure 2
# ------------------------------------------------------------------------------------

pairs = {'41': {'dom':'A4', 'nucl':'32', 'rna':'B'},
         '100':{'dom':'A4', 'nucl':'31', 'rna':'B'},
         '46': {'dom':'A4', 'nucl':'25', 'rna':'B'},
         '34': {'dom':'A3', 'nucl':'33', 'rna':'B'}}

contact_output = input("Please, enter the name of the output file to store the durations\n of contact: ")
contact_output = overwrite_file(contact_output)

duration = int(input("Please, enter the duration of the dynamic (in ns):"))

contacts = contactTime(pairs, frames, threshold, duration, dist_mode, contact_output)

