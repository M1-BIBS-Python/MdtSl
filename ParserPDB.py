#!/usr/bin/env python
# -*- coding : utf8 -*-

def PDBparser(pdbFile):
    """
    Parse un fichier pdb au format ATOM en un dictionnaire utilisable par Python.
    :param pdbFile: Fichier pdb (format ATOM) contenant les coordonnees des atomes d'une proteine.
    :return: Un dictionnaire utilisable par Python.
    """
    #~ infile = open(pdbFile, "r")
    #~ lines = infile.readlines
    
    with open(pdbFile) as f:
	    molecule = {}  # dictionnaire le plus externe
	    chainList = []
	
	    cptAlt = False
	    
	    for line in f:
	        if line[:4:] == 'ATOM':  # si la ligne commence par 'ATOM'
	            if cptAlt == False:
	                alt = line[16]
	                cptAlt = True
	
	            if line[16] == alt:
	                chaine = line[72:76].strip()
	                
	                if chaine not in ["BWAT", "POT", "CL"]:
	    
		                if chaine not in chainList:
		                    chainList.append(chaine)
		                    molecule[chaine] = {}
		                    resList = []
		    
		                curres = line[22:26].strip()
		    
		                if curres not in resList:
		                    resList.append(curres)
		                    molecule[chaine][curres] = {}
		                    atomList = []
		                    molecule[chaine][curres]['resname'] = line[17:20].strip()
		    
		                atom = line[12:16].strip()
		                if atom not in atomList:
		                    atomList.append(atom)
		                    molecule[chaine][curres][atom] = {}
		    
		                molecule[chaine]['reslist'] = resList
		                molecule[chaine][curres]['atomlist'] = atomList
		    
		                molecule[chaine][curres][atom]['x'] = float(line[30:38])
		                molecule[chaine][curres][atom]['y'] = float(line[38:46])
		                molecule[chaine][curres][atom]['z'] = float(line[46:54])
		                molecule[chaine][curres][atom]['id'] = line[6:11].strip()

    return molecule


def PDBparserConf(list):
    """
    Parse une liste contenant les donnees du fichier pdb correspondant a une conformation d'une proteine.
    :param list: Liste contenant les donnees du fichier pdb pour une conformation de proteine.
    :return: Un dictionnaire utilisable par Python.
    """
    molecule = {}  # dictionnaire le plus externe
    chainList = []

    cptAlt = False

    for line in list:
        if line[:4:] == 'ATOM':  # si la ligne commence par 'ATOM'
            if cptAlt == False:
                alt = line[16]
                cptAlt = True

            if line[16] == alt:
                chaine = line[72:76]
                
                if chaine not in ["BWAT", "POT", "CL"]:

	                if chaine not in chainList:
	                    chainList.append(chaine)
	                    molecule[chaine] = {}
	                    resList = []
	
	                curres = line[22:26].strip()
	
	                if curres not in resList:
	                    resList.append(curres)
	                    molecule[chaine][curres] = {}
	                    atomList = []
	                    molecule[chaine][curres]['resname'] = line[17:20].strip()
	
	                atom = line[12:16].strip()
	                if atom not in atomList:
	                    atomList.append(atom)
	                    molecule[chaine][curres][atom] = {}
	
	                molecule[chaine]['reslist'] = resList
	                molecule[chaine][curres]['atomlist'] = atomList
	
	                molecule[chaine][curres][atom]['x'] = float(line[30:38])
	                molecule[chaine][curres][atom]['y'] = float(line[38:46])
	                molecule[chaine][curres][atom]['z'] = float(line[46:54])
	                molecule[chaine][curres][atom]['id'] = line[6:11].strip()

    return molecule



def PDBparserMulti(pdbFile):
    """
    Parse un fichier pdb au format ATOM te contenant plusieurs proteines en un dictionnaire. 
    :param pdbFile: Fichier pbd (format ATOM) contenant plusieurs proteines.
    :return: Un dictionnaire utilisable par Python.
    """
    with open(pdbFile) as f:

	    frames = {}     # dictionnaire contenant toutes les conformations
	    nb_frames = 0   # nombre de conformations
	
	    conf = []       # donnees du pdb correspondant a la conformation
	    model = ""
	
	    for line in f:
	        if "MODEL" in line:
	            if model != "":
	                frames[model] = PDBparserConf(conf)
	
	            conf = []
	            model = line[10:14].strip()
	            nb_frames += 1
	
	        else:
	            conf.append(line)
	
	    frames[model] = PDBparserConf(conf)     # on parse la derniere conformation
	
    return frames
