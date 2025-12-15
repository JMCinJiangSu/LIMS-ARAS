#-*- coding:gbk -*-
import copy
import re

'''
Discription
	
	该脚本用来获取可能获益临床试验信息。
	规则：
	1. 中国临床试验3条+外国临床试验3条
	2. 优先选择正在招募的、临床分期靠后的 

'''
# 用于去除治疗方案中的"Drug: /Biological: /Procedure: "等字符
def stran_druginfo(druginfo):
	regimen_dict = {
		"Drug" : "_Drug",
		"Biological" : "_Biological",
		"Procedure" : "_Procedure",
		"Other" : "_Other",
		"Device" : "_Device",
		"Behavioral" : "_Behavioral",
		"Radiation" : "_Radiation",
		"Diagnostic Test" : "_Diagnostic Test",
		"Genetic" : "_Genetic"
	}
	for k, v in regimen_dict.items():
		if k in druginfo:
			druginfo = druginfo.replace(k, v)
	return druginfo

def getClinic(jsonDict):
	# json返回的是所有检测到的基因的相关临床试验信息，需进行筛选
	clinic_list = copy.deepcopy(jsonDict["clinic_trial"])
	gene_list = []
	result = []
	for i in clinic_list:
		if i["gene_symbol"] not in gene_list:
			gene_list.append(i["gene_symbol"])

	for gene in gene_list:
		clinic_cn = [i for i in clinic_list if re.search("CTR", i["clinicaltrial_number"]) and i["gene_symbol"] == gene]
		clinic_en = [i for i in clinic_list if re.search("NCT", i["clinicaltrial_number"]) and "Drug" in i["interventions"] and i["gene_symbol"] == gene]
		
		# 中国临床试验信息，排序，取临床分期靠后的三条
		if clinic_cn:
			clinic_cn = sorted(clinic_cn, key=lambda i:i["phase"])
			for i in clinic_cn:
				i["interventions"] = [i["interventions"]]
			result += clinic_cn if len(clinic_cn) <= 3 else clinic_cn[0:3]

		# 外国临床试验信息，排序，取临床分期靠后且有治疗药物的三条
		# 治疗方案格式化
		# 临床分期改为中文格式
		phase_dict = {
			"Phase 1" : "I期",
			"Phase 2" : "II期",
			"Phase 3" : "III期",
			"Phase 4" : "IV期",
			"Phase 1 Phase 2" : "I期/II期",
			"Phase 2 Phase 3" : "II期/III期",
			"Early Phase 1" : "Early I期",
			"Phase 1;Phase 2" : "I期/II期"
		}
		if clinic_en:
			clinic_en = sorted(clinic_en, key=lambda i:i["phase"], reverse = True)
			for i in clinic_en:
				druginfo = stran_druginfo(i["interventions"])
				druglist = []
				for j in re.split("_", druginfo):
					if re.search("Drug", j):
						j = j.replace("Drug: ","")
						druglist.append(j)
				i["interventions"] = druglist
				i["phase"] = phase_dict.get(i["phase"], i["phase"])
			result += clinic_en if len(clinic_en) <= 3 else clinic_en[0:3]
			
	return result