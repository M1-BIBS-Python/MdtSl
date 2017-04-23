#!/usr/bin/env python
# -*- coding : utf8 -*-

from ParserPDB import *
from RMSD import *
from computeInterface import *


# ----------------------------------------
# Parsing de la conformation de reference
# ----------------------------------------

ref = PDBparser("pab21_structure_de_ref.pdb")
#writePDBparsed(ref, "ref.txt")


# ------------------------------
# Parsing des 500 conformations
# ------------------------------

frames = PDBparserMulti("pab21_prod_solute_500frames.pdb")
#writePDBparsed(frames, "frames.txt")


# --------------------------------------------------------------------------------------------
# Calcul du RMSD entre la structure de reference et chacune des conformations de la dynamique
# --------------------------------------------------------------------------------------------
"""
f = open("rmsd.txt", "w")

absc_plot = []
ord_global = []
ord_A1 = []
ord_A2 = []
ord_A3 = []
ord_A4 = []
for model in sorted(map(int, frames.keys())):
	absc_plot.append(model)
	rmsd_global = RMSD_prot(ref, frames[str(model)], "CM")
	ord_global.append(rmsd_global)
	f.write("Model "+str(model)+"\t"+str(rmsd_global)+"\n")
	for dom in ['A1', 'A2', 'A3', 'A4']:
		rmsd_dom = RMSD_domain(ref[dom], frames[str(model)][dom], "CM")
		if dom == 'A1':
			ord_A1.append(rmsd_dom)
		elif dom == 'A2':
			ord_A2.append(rmsd_dom)
		elif dom == 'A3':
			ord_A3.append(rmsd_dom)
		else:
			ord_A4.append(rmsd_dom)
		f.write("\t"+str(dom)+"\t"+str(rmsd_dom)+"\n")
f.close()


plt.plot(absc_plot,ord_global)
plt.xlabel('Frame')
plt.ylabel('RMSD')
plt.title('Global RMSD of the protein')
plt.show()

a1, = plt.plot(absc_plot, ord_A1, 'b', label='A1')
a2, = plt.plot(absc_plot, ord_A2, 'r', label='A2')
a3, = plt.plot(absc_plot, ord_A3, 'g', label='A3')
a4, = plt.plot(absc_plot, ord_A4, 'y', label='A4')
plt.xlim(0,5000)
plt.xlabel('Frame')
plt.ylabel('RMSD')
plt.title('RMSD for each domain')
plt.legend(handles = [a1, a2, a3, a4], fontsize = 8)
plt.show()
"""

# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Calcul de la matrice des distances entre chaque domaine proteique et l'ARN
# ---------------------------------------------------------------------------

distMatrix(ref, 'A1', 'B', "CM", heatmap="yes")
distMatrix(ref, 'A2', 'B', "CM", heatmap="yes")
distMatrix(ref, 'A3', 'B', "CM", heatmap="yes")
distMatrix(ref, 'A4', 'B', "CM", heatmap="yes")