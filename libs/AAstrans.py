#-*- coding:gbk -*-
import re

'''
	输入hgvs_p和突变类型
	返回变异描述：
		“该突变导致基因编码蛋白第2471位氨基酸由苏氨酸突变为组氨酸并于2474位发生提前终止” 

'''

def strans(hgvs_p,var_type):
	AA_dict = {"A" : "丙氨酸",
			   "C" : "半胱氨酸",
			   "D" : "天冬氨酸",
			   "E" : "谷氨酸",
			   "F" : "苯丙氨酸",
			   "G" : "甘氨酸",
			   "H" : "组氨酸",
			   "I" : "异亮氨酸",
			   "K" : "赖氨酸",
			   "L" : "亮氨酸",
			   "M" : "蛋氨酸",
			   "N" : "天冬酰胺",
			   "P" : "脯氨酸",
			   "Q" : "谷氨酰胺",
			   "R" : "精氨酸",
			   "S" : "丝氨酸",
			   "T" : "苏氨酸",
			   "V" : "缬氨酸",
			   "W" : "色氨酸",
			   "Y" : "酪氨酸",
			   "*" : "终止密码"}
	outstr = ""
	mat = re.compile(r"\d+")
	#错义突变和无义突变
	if (var_type == "nonSynonymous_Substitution" or var_type == "Nonsense_Mutation") and not re.search("del|ins|dup", hgvs_p):
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in list(hgvs_p) if i in AA_dict.keys()]
		outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位氨基酸由"+hgvs_p_stran[0]+"突变为"+hgvs_p_stran[-1]
	#移码突变
	if re.search("fs*", hgvs_p):
		for i in mat.findall(re.split("fs\*", hgvs_p)[0]):
			pos = i
		endpos = int(pos) + int(re.split("fs\*", hgvs_p)[-1].replace(")", "")) -1
		hgvs_p_stran = [AA_dict[i] for i in re.split("fs*", hgvs_p)[0] if i in AA_dict.keys()]
		outstr = "该突变导致基因编码蛋白第"+pos+"位氨基酸由"+hgvs_p_stran[0]+"突变为"+hgvs_p_stran[-1]+"并于"+str(endpos)+"位发生提前终止"
	#延伸突变
	if var_type == "Extension":
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in re.split("ext", hgvs_p)[0].replace("*", "") if i in AA_dict.keys()]
		outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位终止密码突变为"+hgvs_p_stran[0]+"之后形成新的终止密码"
	#非移码突变
	#缺失
	if re.search("del", hgvs_p) and not re.search("delins", hgvs_p):
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in hgvs_p if i in AA_dict.keys()]
		if len(pos_list) == 1:
			outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位氨基酸缺失"
		else:
			outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位到"+pos_list[-1]+"位氨基酸缺失"
	#插入
	if re.search("ins", hgvs_p) and not re.search("delins", hgvs_p):
		pos_list = [i for i in mat.findall(re.split("ins", hgvs_p)[0])]
		insAA = re.split("ins", hgvs_p)[-1].replace(")","")
		if re.search("^[0-9]*$", insAA):
			outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位到"+pos_list[1]+"位氨基酸之间插入"+insAA+"个氨基酸"
		else:
			hgvs_p_stran = [AA_dict[i] for i in insAA if i in AA_dict.keys()]
			outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位到"+pos_list[1]+"位氨基酸之间插入"+"、".join(hgvs_p_stran)
	#重复
	if re.search("dup", hgvs_p):
		pos_list = [i for i in mat.findall(hgvs_p)]
		hgvs_p_stran = [AA_dict[i] for i in hgvs_p if i in AA_dict.keys()]
		if len(pos_list) == 1:
			outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位氨基酸重复"
		else:
			outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位到"+pos_list[1]+"位氨基酸重复"
	#缺失并插入
	if re.search("delins", hgvs_p):
		pos_list = [i for i in mat.findall(re.split("delins", hgvs_p)[0])]
		delinsAA = re.split("delins", hgvs_p)[-1].replace(")","")
		if len(pos_list) > 1:
			if re.search("^[0-9]*$", delinsAA):
				outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位到"+pos_list[1]+"位氨基酸缺失并插入"+delinsAA+"个氨基酸"
			else:
				hgvs_p_stran = [AA_dict[i] for i in delinsAA if i in AA_dict.keys()]
				outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位到"+pos_list[1]+"位氨基酸缺失并插入"+"、".join(hgvs_p_stran)
		else:
			if re.search("^[0-9]*$", delinsAA):
				outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位氨基酸缺失并插入"+delinsAA+"个氨基酸"
			else:
				hgvs_p_stran = [AA_dict[i] for i in delinsAA if i in AA_dict.keys()]
				outstr = "该突变导致基因编码蛋白第"+pos_list[0]+"位氨基酸缺失并插入"+"、".join(hgvs_p_stran)
	#内含子和剪切突变
#	if var_type == "Intronic" or var_type == "Splicing":
#		outstr = "该突变可能造成异常剪接"

	return outstr
