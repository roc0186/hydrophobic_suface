# color_h
# -------
 
# PyMOL command to color protein molecules according to the Eisenberg hydrophobicity scale
 
#
# Source: http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html
# Amino acid scale: Normalized consensus hydrophobicity scale
# Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
# Reference: J. Mol. Biol. 179:125-142 (1984)
#
# Amino acid scale values:
#
# Ala:  0.620
# Arg: -2.530
# Asn: -0.780
# Asp: -0.900
# Cys:  0.290
# Gln: -0.850
# Glu: -0.740
# Gly:  0.480
# His: -0.400
# Ile:  1.380
# Leu:  1.060
# Lys: -1.500
# Met:  0.640
# Phe:  1.190
# Pro:  0.120
# Ser: -0.180
# Thr: -0.050
# Trp:  0.810
# Tyr:  0.260
# Val:  1.080
#
# Usage:
# color_h (selection)
#
from pymol import cmd

def color_h(selection='all'):
        s = str(selection)
        print(s)
        cmd.set_color('color_ile',[0.996,0.062,0.062])
        cmd.set_color('color_phe',[0.996,0.109,0.109])
        cmd.set_color('color_val',[0.992,0.156,0.156])
        cmd.set_color('color_leu',[0.992,0.207,0.207])
        cmd.set_color('color_trp',[0.992,0.254,0.254])
        cmd.set_color('color_met',[0.988,0.301,0.301])
        cmd.set_color('color_ala',[0.988,0.348,0.348])
        cmd.set_color('color_gly',[0.984,0.394,0.394])
        cmd.set_color('color_cys',[0.984,0.445,0.445])
        cmd.set_color('color_tyr',[0.984,0.492,0.492])
        cmd.set_color('color_pro',[0.980,0.539,0.539])
        cmd.set_color('color_thr',[0.980,0.586,0.586])
        cmd.set_color('color_ser',[0.980,0.637,0.637])
        cmd.set_color('color_his',[0.977,0.684,0.684])
        cmd.set_color('color_glu',[0.977,0.730,0.730])
        cmd.set_color('color_asn',[0.973,0.777,0.777])
        cmd.set_color('color_gln',[0.973,0.824,0.824])
        cmd.set_color('color_asp',[0.973,0.875,0.875])
        cmd.set_color('color_lys',[0.899,0.922,0.922])
        cmd.set_color('color_arg',[0.899,0.969,0.969])
        cmd.color("color_ile","("+s+" and resn ile)")
        cmd.color("color_phe","("+s+" and resn phe)")
        cmd.color("color_val","("+s+" and resn val)")
        cmd.color("color_leu","("+s+" and resn leu)")
        cmd.color("color_trp","("+s+" and resn trp)")
        cmd.color("color_met","("+s+" and resn met)")
        cmd.color("color_ala","("+s+" and resn ala)")
        cmd.color("color_gly","("+s+" and resn gly)")
        cmd.color("color_cys","("+s+" and resn cys)")
        cmd.color("color_tyr","("+s+" and resn tyr)")
        cmd.color("color_pro","("+s+" and resn pro)")
        cmd.color("color_thr","("+s+" and resn thr)")
        cmd.color("color_ser","("+s+" and resn ser)")
        cmd.color("color_his","("+s+" and resn his)")
        cmd.color("color_glu","("+s+" and resn glu)")
        cmd.color("color_asn","("+s+" and resn asn)")
        cmd.color("color_gln","("+s+" and resn gln)")
        cmd.color("color_asp","("+s+" and resn asp)")
        cmd.color("color_lys","("+s+" and resn lys)")
        cmd.color("color_arg","("+s+" and resn arg)")
cmd.extend('color_h',color_h)

def color_h2(selection='all'):
        s = str(selection)
        print(s)
        cmd.set_color("color_ile2",[0.938,1,0.938])
        cmd.set_color("color_phe2",[0.891,1,0.891])
        cmd.set_color("color_val2",[0.844,1,0.844])
        cmd.set_color("color_leu2",[0.793,1,0.793])
        cmd.set_color("color_trp2",[0.746,1,0.746])
        cmd.set_color("color_met2",[0.699,1,0.699])
        cmd.set_color("color_ala2",[0.652,1,0.652])
        cmd.set_color("color_gly2",[0.606,1,0.606])
        cmd.set_color("color_cys2",[0.555,1,0.555])
        cmd.set_color("color_tyr2",[0.508,1,0.508])
        cmd.set_color("color_pro2",[0.461,1,0.461])
        cmd.set_color("color_thr2",[0.414,1,0.414])
        cmd.set_color("color_ser2",[0.363,1,0.363])
        cmd.set_color("color_his2",[0.316,1,0.316])
        cmd.set_color("color_glu2",[0.27,1,0.27])
        cmd.set_color("color_asn2",[0.223,1,0.223])
        cmd.set_color("color_gln2",[0.176,1,0.176])
        cmd.set_color("color_asp2",[0.125,1,0.125])
        cmd.set_color("color_lys2",[0.078,1,0.078])
        cmd.set_color("color_arg2",[0.031,1,0.031])
        cmd.color("color_ile2","("+s+" and resn ile)")
        cmd.color("color_phe2","("+s+" and resn phe)")
        cmd.color("color_val2","("+s+" and resn val)")
        cmd.color("color_leu2","("+s+" and resn leu)")
        cmd.color("color_trp2","("+s+" and resn trp)")
        cmd.color("color_met2","("+s+" and resn met)")
        cmd.color("color_ala2","("+s+" and resn ala)")
        cmd.color("color_gly2","("+s+" and resn gly)")
        cmd.color("color_cys2","("+s+" and resn cys)")
        cmd.color("color_tyr2","("+s+" and resn tyr)")
        cmd.color("color_pro2","("+s+" and resn pro)")
        cmd.color("color_thr2","("+s+" and resn thr)")
        cmd.color("color_ser2","("+s+" and resn ser)")
        cmd.color("color_his2","("+s+" and resn his)")
        cmd.color("color_glu2","("+s+" and resn glu)")
        cmd.color("color_asn2","("+s+" and resn asn)")
        cmd.color("color_gln2","("+s+" and resn gln)")
        cmd.color("color_asp2","("+s+" and resn asp)")
        cmd.color("color_lys2","("+s+" and resn lys)")
        cmd.color("color_arg2","("+s+" and resn arg)")
cmd.extend('color_h2',color_h2)



def color_h3(selection='all'):
        s = str(selection)
        print(s)
        cmd.set_color('color_Ala3',[1.00,0.38,0.38])
        cmd.set_color('color_Arg3',[0.00,0.00,1.00])
        cmd.set_color('color_Asn3',[0.89,0.89,1.00])
        cmd.set_color('color_Asp3',[0.83,0.83,1.00])
        cmd.set_color('color_Cys3',[1.00,0.56,0.56])
        cmd.set_color('color_Gln3',[0.85,0.85,1.00])
        cmd.set_color('color_Glu3',[0.92,0.92,1.00])
        cmd.set_color('color_Gly3',[1.00,0.45,0.45])
        cmd.set_color('color_His3',[1.00,0.91,0.91])
        cmd.set_color('color_Ile3',[1.00,0.00,0.00])
        cmd.set_color('color_Leu3',[1.00,0.16,0.16])
        cmd.set_color('color_Lys3',[0.53,0.53,1.00])
        cmd.set_color('color_Met3',[1.00,0.38,0.38])
        cmd.set_color('color_Phe3',[1.00,0.09,0.09])
        cmd.set_color('color_Pro3',[1.00,0.64,0.64])
        cmd.set_color('color_Ser3',[1.00,0.80,0.80])
        cmd.set_color('color_Thr3',[1.00,0.73,0.73])
        cmd.set_color('color_Trp3',[1.00,0.29,0.29])
        cmd.set_color('color_Tyr3',[1.00,0.57,0.57])
        cmd.set_color('color_Val3',[1.00,0.15,0.15])                                                                                                                                                                             
        cmd.color("color_ile3","("+s+" and resn ile)")
        cmd.color("color_phe3","("+s+" and resn phe)")
        cmd.color("color_val3","("+s+" and resn val)")
        cmd.color("color_leu3","("+s+" and resn leu)")
        cmd.color("color_trp3","("+s+" and resn trp)")
        cmd.color("color_met3","("+s+" and resn met)")
        cmd.color("color_ala3","("+s+" and resn ala)")
        cmd.color("color_gly3","("+s+" and resn gly)")
        cmd.color("color_cys3","("+s+" and resn cys)")
        cmd.color("color_tyr3","("+s+" and resn tyr)")
        cmd.color("color_pro3","("+s+" and resn pro)")
        cmd.color("color_thr3","("+s+" and resn thr)")
        cmd.color("color_ser3","("+s+" and resn ser)")
        cmd.color("color_his3","("+s+" and resn his)")
        cmd.color("color_glu3","("+s+" and resn glu)")
        cmd.color("color_asn3","("+s+" and resn asn)")
        cmd.color("color_gln3","("+s+" and resn gln)")
        cmd.color("color_asp3","("+s+" and resn asp)")
        cmd.color("color_lys3","("+s+" and resn lys)")
        cmd.color("color_arg3","("+s+" and resn arg)")
cmd.extend('color_h3',color_h3)

def color_h4(selection='all'):
        #Using the same hydrophobicity values as used in DS，to show the similar fig.
        #
        #The Hydrophobicity Plot is a plot of the amino acids in a protein against their hydrophobicity index. 
        #The hydrophobicity scale is taken from Kyte and Doolittle [Kyte and Doolittle, 1982].
        #The Kyte and Doolittle hydrophobicity values are as follow:
        # PHE  2.8               SER -0.8
        # MET  1.9               PRO -1.6
        # ILE  4.5               TYR -1.3
        # LEU  3.8               HIS -3.2
        # VAL  4.2               GLN -3.5
        # CYS  2.5               ASN -3.5
        # TRP -0.9               GLU -3.5
        # ALA  1.8               LYS -3.9
        # THR -0.7               ASP -3.5
        # GLY -0.4               ARG -4.5 

        s = str(selection)
        print(s)
        cmd.set_color('color_PHE4',[1.000,0.376,0.376])
        cmd.color("color_PHE4","("+s+" and resn PHE)")
        cmd.set_color('color_SER4',[0.824,0.824,1.000])
        cmd.color("color_SER4","("+s+" and resn SER)")
        cmd.set_color('color_MET4',[1.000,0.573,0.573])
        cmd.color("color_MET4","("+s+" and resn MET)")
        cmd.set_color('color_PRO4',[0.643,0.643,1.000])
        cmd.color("color_PRO4","("+s+" and resn PRO)")
        cmd.set_color('color_ILE4',[1.000,0.000,0.000])
        cmd.color("color_ILE4","("+s+" and resn ILE)")
        cmd.set_color('color_TYR4',[0.714,0.714,1.000])
        cmd.color("color_TYR4","("+s+" and resn TYR)")
        cmd.set_color('color_LEU4',[1.000,0.149,0.149])
        cmd.color("color_LEU4","("+s+" and resn LEU)")
        cmd.set_color('color_HIS4',[0.282,0.282,1.000])
        cmd.color("color_HIS4","("+s+" and resn HIS)")
        cmd.set_color('color_VAL4',[1.000,0.063,0.063])
        cmd.color("color_VAL4","("+s+" and resn VAL)")
        cmd.set_color('color_GLN4',[0.220,0.220,1.000])
        cmd.color("color_GLN4","("+s+" and resn GLN)")
        cmd.set_color('color_CYS4',[1.000,0.439,0.439])
        cmd.color("color_CYS4","("+s+" and resn CYS)")
        cmd.set_color('color_ASN4',[0.220,0.220,1.000])
        cmd.color("color_ASN4","("+s+" and resn ASN)")
        cmd.set_color('color_TRP4',[0.800,0.800,1.000])
        cmd.color("color_TRP4","("+s+" and resn TRP)")
        cmd.set_color('color_GLU4',[0.220,0.220,1.000])
        cmd.color("color_GLU4","("+s+" and resn GLU)")
        cmd.set_color('color_ALA4',[1.000,0.596,0.596])
        cmd.color("color_ALA4","("+s+" and resn ALA)")
        cmd.set_color('color_LYS4',[0.133,0.133,1.000])
        cmd.color("color_LYS4","("+s+" and resn LYS)")
        cmd.set_color('color_THR4',[0.847,0.847,1.000])
        cmd.color("color_THR4","("+s+" and resn THR)")
        cmd.set_color('color_ASP4',[0.220,0.220,1.000])
        cmd.color("color_ASP4","("+s+" and resn ASP)")
        cmd.set_color('color_GLY4',[0.910,0.910,1.000])
        cmd.color("color_GLY4","("+s+" and resn GLY)")
        cmd.set_color('color_ARG4',[0.000,0.000,1.000])
        cmd.color("color_ARG4","("+s+" and resn ARG)")
cmd.extend('color_h4',color_h4)