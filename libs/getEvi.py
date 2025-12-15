#-*- coding:gbk -*-
import re
from libs.getConfig import senseStran
from libs.getPinyin import topinyin
from libs.getInterRef import getRef_from_inter

'''
Discription	

	将变异对应的治疗方案、证据描述等转化为适合填充的格式。 

'''

# 相同证据描述的合并治疗方案展示
def merge_Predictive_evi(datainfo):
	tmp_dict = {}
	for evi in datainfo:
		tmp_dict.setdefault(evi["evi_interpretation"], [])
		tmp_dict[evi["evi_interpretation"]].append(
			{
				"regimen_name" : evi["regimen_name"],
				"evi_conclusion_simple" : evi["evi_conclusion_simple"],
				"clinical_significance_cn" : evi["clinical_significance_cn"],
				"regimen_name_py" : evi["regimen_name_py"],
				"evi_conclusion" : evi["evi_conclusion"]
			}
		)

	merge_result = []
	for k, v  in tmp_dict.items():
		merge_result.append(
			{
				"regimen_name" : "、".join([i["regimen_name"] for i in v]),
				"evi_conclusion_simple" : "/".join([i["evi_conclusion_simple"] for i in v]),
				"clinical_significance_cn" : "/".join([i["clinical_significance_cn"] for i in v]),
				"regimen_name_py" : "/".join([i["regimen_name_py"] for i in v]),
				"evi_interpretation" : k,
				"evi_conclusion" : "/".join([i["evi_conclusion"] for i in v])
			}
		)

	return merge_result

def varRegimen(jsonDict, evi_sum, config, var):
	data = {}
	data["refer_evi"] = []
	data["evi_split"] = {}
	data["refer_evi_risk"] = []

	# 1. 新增字段、证据描述基础处理、排序等
	for evi in evi_sum:
		## 新增字段-用于报告展示-1.临床意义转化为中文
		evi["clinical_significance_cn"] = senseStran(config).get(evi["clinical_significance"], evi["clinical_significance"])
		## 新增字段-用于报告展示-2.等级由A0、A1等转化为A
		evi["evi_conclusion_simple"] = evi["evi_conclusion"][0] if evi["evi_conclusion"] else ""
		## 新增字段-用于报告展示-3.适用浙肿，治疗方案获批机构，仅保留FDA、MNPA、EMA
		evi["regimen_refer_agency_ZJZL"] = "/".join(list(set(["FDA","NMPA","EMA"]) & set(re.split(",", evi["regimen_refer_agency"])))) \
										   if "regimen_refer_agency" in evi.keys() and evi["regimen_refer_agency"] \
										   else "-"
		
		## 新增字段-用于排序-1. 治疗方案转化为拼音，便于后续排序
		evi["regimen_name_py"] = topinyin(evi["regimen_name"]) if evi["regimen_name"] else "0"
		## 新增字段-用于排序-2. 用于对敏感和耐药进行排序，敏感前耐药后
		evi["sense_rule"] = "0" if re.search("Sensitive", evi["clinical_significance"]) else \
							"1" if re.search("Resistant", evi["clinical_significance"]) else \
							evi["clinical_significance"]
		
		# 其他功能-1.证据描述去掉末尾空格
		evi["evi_interpretation"] = evi["evi_interpretation"].strip() if evi["evi_interpretation"] else ""
		
	# 对证据进行排序，按治疗方案等级、敏感/耐药、治疗方案拼音（字母统一大写，便于排序）排序
	#evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	# RET 融合，用药排序更新-由之前的卡博替尼、普拉替尼、赛普替尼更新为普拉提尼、赛普替尼、卡博替尼-2023.05.11
	drug_list_rule = {
		"普拉替尼" : 0, 
		"塞普替尼" : 1,
		"卡博替尼" : 2
		}
	if "bio_category" in var.keys() and var["bio_category"] in ["Sv", "PSeqRnaSv"] and (var["five_prime_gene"] == "RET" or var["three_prime_gene"] == "RET"):
		evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], drug_list_rule.get(i["regimen_name"], 3) , i["regimen_name_py"].upper()))
	else:
		evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	# 排序更新完成-2023.05.11

	# 2. 对排序后的证据进行参考文献提取、证据分类、治疗证据合并等
	for evi in evi_sum:
		## 获取参考文献-证据描述
		data["refer_evi"].extend(getRef_from_inter(jsonDict, evi["evi_interpretation"]))

		## 获取参考文献-遗传风险
		if evi["evidence_type"] == "Predisposing":
			data["refer_evi_risk"].extend(getRef_from_inter(jsonDict, evi["evi_interpretation"]))
		
		## 创建evi_split，将证据进行分类，治疗、辅助诊断、预后、遗传风险
		if evi["evidence_type"] not in data["evi_split"].keys():
			data["evi_split"].setdefault(evi["evidence_type"], [])
		data["evi_split"][evi["evidence_type"]].append(evi)

		## 目前仅有Predictive时才要合并相同证据，额外设置“Predictive_merge”字段，根据报告需求进行选用
		if "Predictive" in data["evi_split"].keys():
			data["evi_split"]["Predictive_merge"] = merge_Predictive_evi(data["evi_split"]["Predictive"])

	# 一些特殊需求
	## 1. 汇总所有证据
	data["regimen_evi_sum"] = evi_sum
	## 2. 汇总A级敏感治疗方案
	data["regimen_FDA_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"]) and var["evi_conclusion_simple"] == "A"
		]
	## 3. 汇总非A级敏感治疗方案
	data["regimen_noFDA_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"]) and var["evi_conclusion_simple"] != "A"
		]

	## 4. 汇总敏感治疗方案-新增证据描述-2023.08.15
	data["regimen_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"],
			"evi_interpretation" : var["evi_interpretation"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"])
		]

	## 5. 汇总耐药治疗方案-新增证据描述-2023.08.15
	data["regimen_R"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"],
			"evi_interpretation" : var["evi_interpretation"]
			} 
		for var in evi_sum if re.search("Resistant",var["clinical_significance"])
		]

	## 6. 敏感治疗方案列表转化为字符串
	data["regimen_S_str"] = "、".join([i["regimen_name"] for i in data["regimen_S"]])

	## 7. 耐药治疗方案列表转化为字符串
	data["regimen_R_str"] = "、".join([i["regimen_name"] for i in data["regimen_R"]])

	# 孙逸仙：药物拆分
	# I类 FDA/NMPA/NCCN药物：A0/A1/A2/C3 
	# I类 临床试验药物：B1/B2/B3/C1/C2/C4/D1/D2/D3/D4/D5/D6
	# II类 FDA/NMPA批准在其他癌种的药物： C3
	# II类 临床试验药物： C1/C2/C4/D1/D2/D3/D4/D5/D6
	# 可总结为两个字段，获批和临床试验药物，I/II类可根据变异等级来判定，应该不会混淆
	SYX_regimen_appr = [evi for evi in evi_sum if evi["evi_conclusion"] in ["A0","A1","A2","C3"]]
	SYX_regimen_clinic = [evi for evi in evi_sum if evi["evi_conclusion"] in ["B1","B2","B3","C1","C2","C4","D1","D2","D3","D4","D5","D6"]]
	#根据治疗方案、预后、辅助诊断、风险等拆分
	data["SYX_regiman_appr"] = {}
	for evi in SYX_regimen_appr:
		if evi["evidence_type"] not in  data["SYX_regiman_appr"].keys():
			data["SYX_regiman_appr"].setdefault(evi["evidence_type"], [])
		data["SYX_regiman_appr"][evi["evidence_type"]].append(evi)
		
	data["SYX_regimen_clinic"] = {}
	for evi in SYX_regimen_clinic:
		if evi["evidence_type"] not in  data["SYX_regimen_clinic"].keys():
			data["SYX_regimen_clinic"].setdefault(evi["evidence_type"], [])
		data["SYX_regimen_clinic"][evi["evidence_type"]].append(evi)

	return data