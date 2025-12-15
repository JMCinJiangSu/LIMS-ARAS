#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.getEvi import varRegimen
from libs.specialRequest import varInfo_FJZL
from libs.rule import S_function
import copy
from libs.rule import decimal_float, decimal_percen
from libs.getConfig import get_gene_class

'''
Discription
	
	处理cnv格式。 

'''

def process_cnv(jsonDict, config):
	cnv = copy.deepcopy(jsonDict["cnv"])
	for var in cnv:
		#var["cn_mean"] = format(round(float(var["cn_mean"]), 2)+0.00, ".2f") if "cn_mean" in var.keys() and var["cn_mean"] else 0
		var["cn_mean"] = decimal_float(var["cn_mean"]) if "cn_mean" in var.keys() and var["cn_mean"] else 0
		
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		
		#AD3101 新增一个字段，只有致病性解读的返回致病性等级，其他情况按致癌性解读等级输出,2025年3月6日
		#var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 3
		#var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 3
		# 更新ad3101 嵇梦晨 2025年6月18日，修改默认值
		var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 0
		var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 0

		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		# 福建肿瘤：返回变异频率相关信息（来源配置表）和治疗方案汇总
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
	# 按cn_mean排序下
	cnv = sorted(cnv, key=lambda i:float(i["cn_mean"]), reverse=True)

	return cnv