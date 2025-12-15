#-*- coding:gbk -*-
import copy
from libs import listResultToDict
from libs.getEvi import varRegimen

'''
Discription
	
	该脚本用来获取HRD结果。
	提取内容：
	1. 结果
	2. 治疗方案
	3. 合并BRCA突变的治疗方案（不包含解释） 

'''

def getHRD(jsonDict, BRCA_data, config):
	hrd_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["hrd"]))
	if "evi_sum" in hrd_dict.keys() and hrd_dict["evi_sum"]:
		# hrd治疗方案转化
		hrd_dict["evi_sum"] = varRegimen(jsonDict, hrd_dict["evi_sum"], config, hrd_dict)
		# hrd治疗方案汇总
		regimen_sum = [{"regimen_name":i["regimen_name"], "evidence_type":i["evidence_type"], "clinical_significance_cn":i["clinical_significance_cn"], "evi_conclusion_simple":i["evi_conclusion_simple"],"regimen_name_py":i["regimen_name_py"]} for i in hrd_dict["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Diagnostic","Predictive","Prognostic"]]
		# BRCA治疗方案汇总
		if BRCA_data:
			for var in BRCA_data:
				regimen_sum += [{"regimen_name":i["regimen_name"], "evidence_type":i["evidence_type"], "clinical_significance_cn":i["clinical_significance_cn"], "evi_conclusion_simple":i["evi_conclusion_simple"],"regimen_name_py":i["regimen_name_py"]} for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Diagnostic","Predictive","Prognostic"]]
		# hrd + BRCA治疗方案去重、排序
		regimen_sum_redup = []
		for i in regimen_sum:
			if i not in regimen_sum_redup:
				regimen_sum_redup.append(i)	
		regimen_sum_redup = sorted(regimen_sum_redup, key=lambda i:(i["evi_conclusion_simple"], i["clinical_significance_cn"], i["regimen_name_py"]))
		# 等级判断（判断依据包含hrd和BRCA的结果）
		regimen_level = [i["evi_conclusion_simple"] for i in regimen_sum_redup]
		hrd_dict["level_num"] = 5 if set(["A", "B"]) & set(regimen_level) else 4 if set(["C", "D"]) & set(regimen_level) else 3
		# 治疗方案分类展示
		hrd_dict["regimen"] = {}
		for regimen in regimen_sum_redup:
			if regimen["evidence_type"] not in hrd_dict["regimen"]:
				hrd_dict["regimen"].setdefault(regimen["evidence_type"], [])
			hrd_dict["regimen"][regimen["evidence_type"]].append(regimen)

	return hrd_dict