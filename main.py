#!/usr/bin/env python
# -*- coding : utf8 -*-

from ParserPDB import *
from RMSD import *

"""
Parsing de la conformation de reference
"""
ref = PDBparser("pab21_structure_de_ref.pdb")

f = open("ref.txt", "w")
f.write(str(ref))
f.close()


"""
Parsing des 500 conformations
"""
frames = PDBparserMulti("pab21_prod_solute_500frames.pdb")

f = open("frames.txt", "w")
f.write(str(frames))
f.close()

"""
Calcul du RMSD entre la structure de reference et chacune des conformations de la dynamique
"""
ref = PDBparser("pab21_structure_de_ref.pdb")
frames = PDBparserMulti("pab21_prod_solute_500frames.pdb")

f = open("rmsd.txt", "w")
for model in frames.keys():
	rmsd_global = RMSD_prot(ref, frames[model], "CA")
	f.write("Model "+str(model)+"\t"+str(rmsd_global)+"\n")
	for dom in ['A1', 'A2', 'A3', 'A4']:
		rmsd_dom = RMSD_domain(ref[dom], frames[model][dom], "CA")
		f.write("\t"+str(dom)+"\t"+str(rmsd_dom)+"\n")
f.close()