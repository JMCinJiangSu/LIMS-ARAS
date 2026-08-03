#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.getEvi import varRegimen
from libs.rule import S_function
import copy
from libs.specialRequest import varInfo_FJZL
from customize_filters import get_analysis_version, judge_version

'''
Discription
	
	处理hd格式。 (非胚系变异)

'''

def process_hd(jsonDict, config):
	hd = [var for var in copy.deepcopy(jsonDict["hd"])] if "hd" in jsonDict.keys() else []
	# 判断HD是否要展示
	json_batch_name = jsonDict["sample_info"].get("json_batch_name", "")
	judge_hd_inter = False
	for var in hd:
		if "evi_sum" in var.keys():
			judge_hd_inter = True
			break
	# 2026.06.24 XW4221 V2需要展示HD
	if judge_hd_inter and (jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）", "CP200（组织）"] and jsonDict["sample_info"]["product_name"] == "XW4221-Oncopro200"):
		for var in hd:
			var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
								  clinicalNumStran(config).get(var["function_classification"], 3)
			var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
								  functionNumStran(config).get(var["clinical_significance"], 3)
			var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
			var["clinic_num_s"], var["top_level"] = S_function(var)
			# 福建肿瘤：返回变异频率相关信息（来源配置表）和治疗方案汇总
			var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		return hd
	# 以上无需判断分析版本，都输出hd结果.CP 200新流程v0.1.4后展示HD,嵇梦晨，2026.05.18
	elif "ADXHS-OncoPro" in json_batch_name and judge_version([json_batch_name, "v0.1.4"]) and judge_hd_inter:
		for var in hd:
			var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
								  clinicalNumStran(config).get(var["function_classification"], 3)
			var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
								  functionNumStran(config).get(var["clinical_significance"], 3)
			var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
			var["clinic_num_s"], var["top_level"] = S_function(var)
			# 福建肿瘤：返回变异频率相关信息（来源配置表）和治疗方案汇总
			var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		return hd
	else:
		return []
	