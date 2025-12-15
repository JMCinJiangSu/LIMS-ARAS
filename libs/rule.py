#-*- coding:gbk -*-
import re
import copy
from decimal import Decimal

'''
	统一处理模板中的各种规则
	体细胞：
	# I类：治疗/辅助诊断/预后证据最高等级为A/B
	# II类：治疗/辅助诊断/预后证据最高等级为C/D
	# 肿瘤发生发展相关，治疗/辅助诊断/预后未匹配到证据，且变异为致病/疑似致病的
	# III类：治疗/辅助诊断/预后未匹配到证据，且为意义不明位点
'''

## 将常使用的判定条件分离出来
# 1. 返回变异来源，胚系/体细胞，未知来源的按体细胞处理
def origin(var):
	if var["var_origin"] == "germline":
		return "G"
	else:
		return "S"

# 2. 判断是否匹配到治疗/辅助诊断/预后的证据
def judgeRegimen(var):
	if var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()):
		return 1
	else:
		return 0

# 3. 变异基因是否在目标基因列表中（单个基因的可以转成列表格式进行处理）
def matchGene(var, genelist):
	return set(re.split(",", var["gene_symbol"])) & set(genelist)

# 4. 判定变异是否为体系[指定等级]或胚系[指定等级]
def judge_var(var, s_list, g_list):
	if (var["var_origin"] == "germline" and var["clinic_num_g"] in g_list) or (var["var_origin"] != "germline" and var["clinic_num_s"] in s_list):
		return 1
	else:
		return 0

# 5. 返回变异结果，如Snvindel返回hgvs_p或hgvs_c
def get_varsimpleinfo(var):
	if var["bio_category"] == "Snvindel":
		return var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
	elif var["bio_category"] == "Cnv":
		return "扩增"
	elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
		return "{0}-{1}融合".format(var["five_prime_gene"], var["three_prime_gene"])
## 常使用判定条件-结束

# 体细胞变异分为I/II/肿瘤发生发展相关/III类-不分基因 适用于大部分模板
def s_var_rule(var_data):
	result = {}
	result["level_I"] = [var for var in var_data if var["clinic_num_s"] in [5] and origin(var)=="S" and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if var["clinic_num_s"] in [4] and origin(var)=="S" and judgeRegimen(var)]
	result["level_onco_nodrug"] = [var for var in var_data if var["clinic_num_s"] in [5, 4] and origin(var)=="S" and not judgeRegimen(var)]
	result["level_III"] = [var for var in var_data if var["clinic_num_s"] == 3 and origin(var)=="S"]
	# IV类： 适用需要展示变异列表的HRR-20220915新增
	result["level_IV"] =  sorted([var for var in var_data if var["clinic_num_s"] in [1,2] and origin(var)=="S"], key=lambda i:i["clinic_num_s"], reverse=True)
	# III类不包括同义突变-适用福建省立CP40-20221012 "Synonymous_Substitution"
	result["level_III_without_Syn"] = list(filter(lambda i : "type" in i.keys() and i["type"] != "Synonymous_Substitution", result["level_III"]))
	return result

# 体细胞变异分为I/II/肿瘤发生发展相关/III类-分单个基因 适用于PTM、BPTM和部分定制模板
def s_var_rule_gene(var_data, gene):
	result = {}
	result["level_I"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] == 5 and origin(var)=="S" and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] == 4 and origin(var)=="S" and judgeRegimen(var)]
	result["level_onco_nodrug"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] in [5, 4] and origin(var)=="S" and not judgeRegimen(var)]
	result["level_III"] = [var for var in var_data if matchGene(var, [gene]) and var["clinic_num_s"] == 3 and origin(var)=="S"]
	return result

# 体细胞变异分为I/II/肿瘤发生发展相关/III类-分多个基因 适用于上海十院HRR和部分定制模板
def s_var_rule_genelist(var_data, genelist):
	result = {}
	result["level_I"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] == 5 and origin(var) == "S" and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] == 4 and origin(var) == "S" and judgeRegimen(var)]
	result["level_onco_nodrug"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] in [5, 4] and origin(var)=="S" and not judgeRegimen(var)]
	result["level_III"] = [var for var in var_data if matchGene(var, genelist) and var["clinic_num_s"] == 3 and origin(var)=="S"]
	return result

# 胚系分为5/4/3/2/1类
def g_var_rule(var_data):
	result = {}
	for i in range(1, 6):
		result["level_"+str(i)] = [var for var in var_data if var["clinic_num_g"] == i and origin(var) == "G"]
	return result

# 胚系分为5/4/3/2/1类，并匹配基因列表
def g_var_rule_genelist(var_data, genelist):
	result = {}
	for i in range(1, 6):
		result["level_"+str(i)] = [var for var in var_data if var["clinic_num_g"] == i and origin(var) == "G" and matchGene(var, genelist)]
	return result

# 胚系变异根据治疗/辅助诊断/预后证据提取I/II类
def g_var_regimen_rule(var_data):
	result = {}
	result["regimen_level_I"] = [var for var in var_data if var["clinic_num_s"] == 5 and origin(var)=="G" and judgeRegimen(var)]
	result["regimen_level_II"] = [var for var in var_data if var["clinic_num_s"] == 4 and origin(var)=="G" and judgeRegimen(var)]
	return result

# 胚系+体细胞变异根据治疗/辅助诊断/预后证据提取I/II类
def var_regimen_rule(var_data):
	result = {}
	result["level_I"] = [var for var in var_data if var["clinic_num_s"] == 5 and judgeRegimen(var)]
	result["level_II"] = [var for var in var_data if var["clinic_num_s"] == 4 and judgeRegimen(var)]
	return result

# 浙江肿瘤胚系+体细胞有用药建议/预后/辅助诊断变异整理（报告中仅展示用药，无用药且有预后/辅助诊断的写“暂无靶向用药提示”，模板代码中较难实现，故这边再加字段）
def process_result_regimen_ZJZL(var_list_zjzl):
	for var in var_list_zjzl:
		if "Predictive" not in var["evi_sum"]["evi_split"]:
			var["evi_sum"]["evi_split"]["Predictive"] = {"note" : "without_regimen"}
	return var_list_zjzl

# PTM/BPTM需要展示各个基因的情况，因为TP53、POLE、BRCA基因都有证据，不存在肿瘤发生发展相关变异，有的话，可以在这边改。
def var_bptm_rule(var_data, gene):
	result = {}
	# 格式1：I/II、III类
	result[gene+"_level12"] = [var for var in var_data if var["clinic_num_s"] in [5, 4] and origin(var)=="S" and judgeRegimen(var) and matchGene(var, [gene])]
	result[gene+"_level3"] = [var for var in var_data if var["clinic_num_s"] == 3 and origin(var)=="S" and matchGene(var, [gene])]
	# 肿瘤发生发展相关变异放到III类中-2023.09.12
	result[gene+"_level3"] = [var for var in var_data if (var["clinic_num_s"] == 3 or (var["clinic_num_s"] in [5, 4]) and not judgeRegimen(var)) and origin(var)=="S" and matchGene(var, [gene])]
	# 2023.09.12-更新完成

	# 格式2：与格式1区别在于如果未检出变异，需要标记
	result[gene+"_level12_withECtype"] = copy.deepcopy(result[gene+"_level12"]) if result[gene+"_level12"] else [{"gene_symbol" : gene, "result" : "nofound"}]
	result[gene+"_level3_withECtype"] = copy.deepcopy(result[gene+"_level3"]) if result[gene+"_level3"] else [{"gene_symbol" : gene, "result" : "nofound"}]
	# 格式3：I/II/III结果放一起
	level_data = s_var_rule_gene(var_data, gene)
	result[gene+"_withEC_type"] = level_data["level_I"] + level_data["level_II"] + level_data["level_onco_nodrug"] + level_data["level_III"]
	result[gene+"_withEC_type"] = result[gene+"_withEC_type"] if result[gene+"_withEC_type"] else [{"gene_symbol" : gene, "result" : "nofound"}]
	return result

# 胚系-林奇综合征：五个基因结果需要分开展示（EPCAM、MLH1、MSH2、MSH6、PMS2），展示上述基因的3、4、5类变异
def var_lyn5_rule(var_data, gene_list):
	result = {}
	for gene in gene_list:
		if gene not in result.keys():
			result.setdefault(gene, [])
		result[gene] = [var for var in var_data if var["clinic_num_g"] in [5, 4, 3] and origin(var)=="G" and matchGene(var, [gene])]
	return result

# 提取未检测到I/II类体细胞变异的基因
def nofoundPath_genelist(var_data, gene_list):
	return sorted(list(set(gene_list) - set([var["gene_symbol"] for var in var_data if var["clinic_num_s"] in [4,5] and origin(var)=="S" and judgeRegimen(var)])))

# 提取未检测到体细胞变异的基因(包含I/II/III/肿瘤发生发展相关)
def nofound_genelist(var_data, gene_list):
	return sorted(list(set(gene_list) - set([var["gene_symbol"] for var in var_data if var["clinic_num_s"] in [3,4,5] and origin(var) == "S"])))

# 致病/疑似致病但无用药的变异，在结果汇总表中归为III类
def S_level(var):
	s_level = var["clinic_num_s"] if judgeRegimen(var) else 3
	return s_level

# 体细胞等级，最高等级AB为I类，CD为II类，无证据则按解读的分类
# 获取证据最高等级
def S_function(var):
	regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive", "Prognostic", "Diagnostic"]] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else []

	clinic_num_s =  5 if set(["A", "B"]) & set(regimen_level) else \
					4 if set(["C", "D"]) & set(regimen_level) else \
					var["clinic_num_s"]
					
	top_level = "A" if "A" in regimen_level else \
				"B" if "B" in regimen_level else \
				"C" if "C" in regimen_level else \
				"D" if "D" in regimen_level else \
				"N"
	return clinic_num_s, top_level

def decimal_float(a):
	'''
	四舍五入，处理浮点数
	'''
	#return str(Decimal(str(a)).quantize(Decimal("0.01"), rounding="ROUND_HALF_UP"))

	# 更新为<5舍去，>=5进1，和excel数据一致，但有可能和系统展示的不同-2023.07.04
	return "{:.2f}".format(float(int(float(a) * 100 + 0.5) / 100))

def decimal_percen(a):
	'''
	四舍五入，处理百分数
	'''
	#return str(Decimal(str(float(a)*100)).quantize(Decimal("0.01"), rounding="ROUND_HALF_UP"))+"%"

	# 更新为<5舍去，>=5进1，和excel数据一致，但有可能和系统展示的不同-2023.07.04
	return "{:.2%}".format(float(int(float(a) * 10000 + 0.5) / 10000))

# 识别是否为数值
def is_number(i):
	try:
		float(i)
		return True
	except:
		pass 
	if i.isnumeric():
		return True
	return False