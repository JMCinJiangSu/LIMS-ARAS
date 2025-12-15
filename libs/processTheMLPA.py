#-*- coding:gbk -*-
import re
from libs.getEvi import varRegimen
import copy
from libs.specialRequest import varInfo_FDZS

'''
Discription
	
	处理MLPA格式。 
	返回格式：
	mlpa : {
		"B1_LOSS" : [],
		"B1_Gain" : [],
		"B2_LOSS" : [],
		"B2_Gain" : []
	}

'''

def process_mlpa(jsonDict, config):
	mlpa = copy.deepcopy(jsonDict["mlpa"])
	result = {}
	for var in mlpa:
		# evi_sum转换
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var) if "evi_sum" in var.keys() and var["evi_sum"] else []
		# 复旦中山变异特殊需求
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], "", "", config)

	# 返回格式
	for gene in ["BRCA1", "BRCA2"]:
		for i in ["Loss", "Gain"]:
			result["B"+gene[-1]+"_"+i] = [var for var in mlpa if var["gene_symbol"] == gene and re.search(i, var["type"])]

	# 返回MLPA图-2022.11.15
	mlpa_image = []
	for var in mlpa:
		if re.search("Loss|Gain", var["type"]):
			if var["file_path"] and var["file_path"] not in mlpa_image:
				mlpa_image.append(var["file_path"])

	# 返回MLPA del图-2023.03.30
	mlpa_image_del = []
	for var in mlpa:
		if re.search("Loss", var["type"]):
			if var["file_path"] and var["file_path"] not in mlpa_image_del:
				mlpa_image_del.append(var["file_path"])
	
	return result, mlpa_image,  mlpa_image_del
		
