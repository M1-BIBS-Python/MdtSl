#!/usr/bin/env python
# -*- coding : utf8 -*-

from ParserPDB import *

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