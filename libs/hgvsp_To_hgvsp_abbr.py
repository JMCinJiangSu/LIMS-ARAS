#-*- coding:gbk -*-
import re

'''
	输入三字母氨基酸，转化为三字母氨基酸 

'''
def splitAA(hgvs_p):
	transAA = {
		"A":"Ala",
		"C":"Cys",
		"D":"Asp",
		"E":"Glu",
		"F":"Phe",
		"G":"Gly",
		"H":"His",
		"I":"Ile",
		"K":"Lys",
		"L":"Leu",
		"M":"Met",
		"N":"Asn",
		"P":"Pro",
		"Q":"Gln",
		"R":"Arg",
		"S":"Ser",
		"T":"Thr",
		"V":"Val",
		"W":"Trp",
		"Y":"Tyr",
		"*":"Ter"
		}
#	AA_str = []
#	for i in list(hgvs_p):
#		if i in transAA.keys():
#			AA_str.append(transAA[i])
#		else:
#			AA_str.append(i)
#		hgvs_p_abbr = "".join(AA_str)	
#	return hgvs_p_abbr
	AA_str = []
	for i in list(hgvs_p):
		AA_str.append(transAA.get(i, i))
#	print (hgvs_p, "".join(AA_str))
	return "".join(AA_str)