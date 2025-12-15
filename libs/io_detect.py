#-*- coding:gbk -*-
#import os
#import xlrd
import re
from libs.getConfig import getconfigxlsx
from libs.rule import judge_var

'''
Discription
	
	免疫正负相关检测结果。 
	适用：
	1. v3临检通用MP IVD
	2. v3临检通用116
	3. v3浙肿Master 跟临检的格式不一样，为了防止出错，额外再处理一下
	
'''
def get_io_detect(var_data, config):
	io = {}
	io["result"], io["io_p_summary"], io["io_n_summary"] = io_detect(var_data)
	io["result_116"], io["io_p_summary_116"], io["io_n_summary_116"] = io_detect_for_116(var_data)
	io["io_p_summary_ZJZL"], io["io_n_summary_ZJZL"] = io_detect_for_ZJZL(var_data, config)
	io["result_new_cnvlist"], io["io_p_summary_new_cnvlist"], io["io_n_summary_new_cnvlist"] = io_detect_new_cnvlist(var_data)
	return io

def io_detect(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1","FGF3","FGF4","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	#io_p_list = []
	#io_n_list = []
	#for k, v in io_result.items():
	#	io_p_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_P])
	#	io_n_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_N])


	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF4/FGF19扩增")
		#print (", ".join(io_n_list))
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_new_cnvlist(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	#io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
	#			 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	#cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]
	# 删除FGF4扩增-2022307.20
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	#both_cnv_list = ["CCND1","FGF3","FGF4","FGF19"]
	# 删除FGF4扩增-2022307.20
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	#io_p_list = []
	#io_n_list = []
	#for k, v in io_result.items():
	#	io_p_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_P])
	#	io_n_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_N])


	# 处理CNV共突变
	#if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
	#	io_n_list.append("CCND1/FGF3/FGF4/FGF19扩增")
	# 删除FGF4扩增-2022307.20
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19扩增")
	#print (", ".join(io_n_list))
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_116(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","FANCA","MRE11","PALB2","MLH1","MSH2","MSH6",\
				 "PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","TERT","CDK12"]
	#io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN"]
	# 新增CCND1/FGF3/FGF19共扩增-2023.07.12
	io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN", "CCND1", "FGF3", "FGF19"]
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合")
			# 加一个孙逸仙的 Gene1-ALK(G8:A12)基因重排 -2023.08.02
			if "ALK_syx" not in io_result.keys():
				io_result.setdefault("ALK_syx", [])
			io_result["ALK_syx"].append("{0}-{1}({2}{3}:{4}{5})基因重排".format(var["five_prime_gene"], var["three_prime_gene"], var["five_prime_gene"][0], \
								var["five_prime_cds"].replace("exon", ""), var["three_prime_gene"][0], var["three_prime_cds"].replace("exon", "")))
			# 新增结束-2023.08.02
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1", "FGF3", "FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19扩增")
	
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_ZJZL(var_data, config):
	fdzs_dict, fjzl_database, Data = getconfigxlsx(config)
	# io_result用于填充IO表
	io_result = {}

	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]
	# 展示I/II类和肿瘤发生发展相关变异
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "cnv", 
					"var_info" : "扩增"
					}
				)
		# 仅展示融合的基因
		# 可能存在两个融合基因都在检测范围里的情况，系统json返回的gene_symbol为gene1,gene2，需额外做识别
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(
				{
					"var_type" : "sv", 
					"three_prime_gene" : var["three_prime_gene"], 
					"three_prime_cds" : var["three_prime_cds"], 
					"five_prime_gene" : var["five_prime_gene"], 
					"five_prime_cds" : var["five_prime_cds"]
					}
				)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "snvindel", 
					"hgvs_p" : var["hgvs_p"], 
					"hgvs_c" : var["hgvs_c"], 
					"var_origin" : var["var_origin"], 
					"gene_region" : var["gene_region"], 
					"transcript_primary" : var["transcript_primary"]
					}
				)
	# summary展示
	io_inter = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF4","FGF19"]:
			if k in Data.keys():
				io_inter.append({
					"gene_symbol" : k,
					"var_info" : v,
					"inter" : Data[k]
				})
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
		io_inter.append({
			"gene_symbol" : "CCND1/FGF3/FGF4/FGF19",
			"var_info" : io_result["CCND1"],
			"inter" : Data["CCND1"] 
		})
	
	io_p_list = [i for i in io_inter if i["gene_symbol"] in io_gene_P]
	io_n_list = [i for i in io_inter if i["gene_symbol"] in io_gene_N or i["gene_symbol"] == "CCND1/FGF3/FGF4/FGF19"]

	return io_p_list, io_n_list