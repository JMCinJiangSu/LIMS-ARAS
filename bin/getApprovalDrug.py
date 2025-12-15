#-*- coding:gbk -*-
import copy
import re
from libs import getDrugSortRule
from functools import reduce
import itertools
'''
Discription
	
	该脚本用来获取获批药物介绍信息。
	1. 需要进行排序（按药物出现先后顺序）
	2. 涉及药物介绍的部分: 变异、MSI、TMB、PD-L1、HRD
	* BRCA只涉及到变异(snvindel、MLPA) 

'''

def getDrug(jsonDict):
	rule = getDrugSortRule.sort_rule(jsonDict)
	drug_list = copy.deepcopy(jsonDict["drug"])
	drug_result = []
	# 删除TMB分子特征
	# 临检MP IVD表中不放TMB-H的药物介绍
	# 只有MP和448有TMB的内容，其他项目无影响
	# 若TMB-H在[var]中，则删除该分子标志物
	tmb_biomarker = {"biomarker_type" : "TMB-H"}
	for drug in drug_list:
		drug["name"] = drug["general_name_cn"] if drug["general_name_cn"] else drug["general_name_en"]
		#var 中的hgvs_p加括号
		if drug["var"]:
			for a in drug["var"]:
				if a and "hgvs_p" in a.keys() and a["hgvs_p"]:
					a["hgvs_p"] = a["hgvs_p"].replace("p.", "p.(")+")" if not re.search("=", a["hgvs_p"]) and a["hgvs_p"] != "p.?" and not re.search("\(", a["hgvs_p"]) else a["hgvs_p"]
				# hgvs_p_abbr添加括号-20220923
				if a and "hgvs_p_abbr" in a.keys() and a["hgvs_p_abbr"]:
					a["hgvs_p_abbr"] = a["hgvs_p_abbr"].replace("p.", "p.(")+")" if not re.search("=", a["hgvs_p_abbr"]) and a["hgvs_p_abbr"] != "p.?" else a["hgvs_p_abbr"]
				if a == tmb_biomarker:
					drug["var"].remove(a)
				# 将MLPA的var统一为gene exon del/gene exon dup(dup 暂时无药)-20220902
				if "biomarker_type" in a.keys() and a["biomarker_type"] and re.search("BRCA", a["biomarker_type"]) and re.search("Loss|Gain", a["biomarker_type"]):
					biomarker_type = re.split(":", a["biomarker_type"])
					a["biomarker_type"] = biomarker_type[1]+" "+biomarker_type[2]+" del" if "Loss" in biomarker_type else biomarker_type[1]+" "+biomarker_type[2]+" dup"

		# 将药物获批机构限定在FDA和NMPA中
		drug["approval_organization"] = list(set(drug["approval_organization"]) & set(["FDA", "NMPA"]))
		# 适应症去重-2022.08.10
		#drug["adaptation_disease_cn"] = reduce(lambda x, y:x if y in x else x + [y], [[],]+drug["adaptation_disease_cn"]) if drug["adaptation_disease_cn"] else []
		# 适应症根据“\n”进行拆分，并去重-2022.08.18
		adaptation = list(itertools.chain(*[re.split("\n", i.strip()) for i in drug["adaptation_disease_cn"]])) if drug["adaptation_disease_cn"] else []
		drug["adaptation_disease_cn"] = reduce(lambda x, y:x if y in x else x + [y], [[],]+adaptation)		

		if str(drug["name"]) in rule and drug["approval_organization"] and drug["var"]:
			drug_result.append(drug)

	# 这边不确定会不会有drug内容超出实际选择药物的范围
	drug_result = sorted(drug_result, key=lambda i:rule.index(i["name"]))

	return drug_result