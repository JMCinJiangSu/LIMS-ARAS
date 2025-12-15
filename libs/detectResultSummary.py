#-*- coding:gbk -*-
import os
import re
from libs.getConfig import CP40Gene
from libs.getConfig import CP40_SplitGene
import copy
from libs.rule import S_level, s_var_rule
from libs.getConfig import PAN116_other_gene

'''
Discription
	
	检测结果汇总，格式多种。 
	
'''
def getSummary_detect(var_data, mlpa_data, config):
	summary_result = {}
	summary_result["format1_forCP"] = format1(var_data, config)
	summary_result["format2_forCP_splitGene"] = {}
	summary_result["format2_forCP_splitGene"]["gene10"], \
	summary_result["format2_forCP_splitGene"]["gene30"], \
	summary_result["format2_forCP_splitGene"]["gene_116_other"] = format2(var_data, config)
	summary_result["format3_forHRDC_QL"] = format3(var_data, mlpa_data)
	summary_result["format4_forHRR_QL"] = format4(var_data, mlpa_data)
	summary_result["format5_for116"] = {}
	summary_result["format5_forLC10"] = {}
	summary_result["format5_for116"]["CRC25"], \
	summary_result["format5_for116"]["TC21"], \
	summary_result["format5_for116"]["GA18"], \
	summary_result["format5_forLC10"]["NBBL"] = format5(var_data)
	summary_result["format6_forLC10_YNZL"] = format6(var_data)
	summary_result["format7_forHRR_YNZL"] = format7(var_data)
	summary_result["format8_forHRDC_YNZL"] = format8(var_data)
	summary_result["format9_forLC10_FJZL"] = format9(var_data)
	summary_result["format10_forYNZL"] = {}
	summary_result["format10_forYNZL"]["CP40"], \
	summary_result["format10_forYNZL"]["TC21"], \
	summary_result["format10_forYNZL"]["GA18"], \
	summary_result["format10_forYNZL"]["CRC12"] = format10(var_data)	
	summary_result["format11_forFDZS_CP40"], summary_result["format11_forFDZS_CP40_sum"] = format11(var_data)
	summary_result["format12_forGDRM_gHRR"] = format12(var_data, mlpa_data)
	summary_result["format13_forBRCA_YNZL"] = format13(var_data, mlpa_data)
	summary_result["format14_forBPTM_YNZL"] = format14(var_data)

	return summary_result

# 函数1：包含基因、检测内容和检测结果，变异类型包含Snvindel、CNV、SV，包含I/II/III类变异
# 主要适用CP40产品
def getSum_forCP40(var_data_copy,gene_dict):
	gene_list = [k for k in gene_dict.keys()]
	result = []
	var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [4,5]) or \
				(var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5])]
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_dict.keys():
				var_info = ""
				note = var["judge_mergeMET"] if "judge_mergeMET" in var.keys() else ""
				note_freq = var["freq_2"] if "judge_mergeMET" in var.keys() else ""
				freq = ""
				if var["bio_category"] == "Snvindel":
					freq = var["freq_str"]
				elif var["bio_category"] == "Cnv":
					freq = var["cn_mean"]
				else:
					if "freq_str" in var.keys() and var["freq_str"]:
						freq = var["freq_str"]
					elif "copies" in var.keys() and var["copies"]:
						freq = var["copies"]
					elif "reads" in var.keys() and var["reads"]:
						freq = var["reads"]
				level = str(int(S_level(var)))
				result.append({
					"gene_symbol" : gene,
					"var_info" : var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合" if var["bio_category"] in ["Sv", "PSeqRnaSv"] else "扩增" if var["bio_category"] == "Cnv" else var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"],
					"freq" : freq,
					"region" : gene_dict[gene],
					"level" : level,
					"note" : note,
					"note_freq" : note_freq,
					"var_origin" : var["var_origin"],
					"freq_rc" : var["freq_rc"] if "freq_rc" in var.keys() and var["freq_rc"] else "",
					"clinic_num_g" : var["clinic_num_g"] if "clinic_num_g" in var.keys() and var["clinic_num_g"] else "",
					"var_detail" : var
				})
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"region" : gene_dict[gene]
		})
	return sorted(result, key=lambda i:i["gene_symbol"])

# 格式1：适用于CP40，如山东齐鲁/山东省立CP40
# 包含基因、检测内容和检测结果，每个变异一行，包含I/II/III类变异
def format1(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	cp40_dict = CP40Gene(config)
	cp40_summary = getSum_forCP40(var_data_copy, cp40_dict)
	return cp40_summary

# 格式2，拆分10基因和30基因，适用于聊城、厦门市一CP40
# 包含基因、检测内容和检测结果，每个变异一行，包含I/II/III类变异
def format2(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	gene10_dict = CP40_SplitGene(config)["gene10"]
	gene30_dict = CP40_SplitGene(config)["gene30"]
	gene116_other = PAN116_other_gene(config)
	gene10_summary = getSum_forCP40(var_data_copy, gene10_dict)
	gene30_summary = getSum_forCP40(var_data_copy, gene30_dict)
	gene_116_other_summary = getSum_forCP40(var_data_copy, gene116_other)
	return gene10_summary, gene30_summary, gene_116_other_summary

# 函数2：包含基因、检测内容、丰度和变异解读，每个变异一行，体细胞/来源不明变异包含I/II/III类，胚系包含345类
# 基因分为已获批基因（展示时需要标星号）和其他基因
# 主要适用胚系产品（如HRR、HRD C等）
def getSum_forGermlinePanel(var_data, mlpa_data, gene_list_all, gene_list_appr):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_copy = copy.deepcopy(mlpa_data)
	result = []
	var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [3,4,5]) or\
				(var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5])]
	detect_gene = [var["gene_symbol"] for var in var_list]
	# 加入MLPA
	for k, v in mlpa_data_copy.items():
		if v:
			var_list += v
			detect_gene += [var["gene_symbol"] for var in v]
	for var in var_list:
		# 限制变异类型为snvindel，不然会报错
		if var["gene_symbol"] in gene_list_all and var["bio_category"] == "Snvindel":
			# 致病/疑似致病但无用药的归为III类
			var["clinic_num_s"] = S_level(var)
			result.append(var)
		else:
			if "type" in var.keys() and var["type"] and re.search("Loss|Gain", var["type"]) and var["gene_symbol"] in gene_list_all:
				result.append(var)
	for gene in set(gene_list_all) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	# 标记获批基因
	for i in result:
		if i["gene_symbol"] in  gene_list_appr:
			i["appr"] = "T"

	return sorted(result, key = lambda i:i["gene_symbol"])

# 格式3: 适用于HRD C，山东齐鲁
# 包含基因、检测内容、丰度和变异解读，每个变异一行，包含I/II/III类变异
# summary基因中分为已获批基因（展示时需要标星号）和其他基因
def format3(var_data, mlpa_data):
	gene_list_all = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "FANCA", "FANCL", "HDAC2", "PALB2", "PPP2R2A", \
					"PTEN", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "TP53"]
	gene_list_appr = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L"]
	return getSum_forGermlinePanel(var_data, mlpa_data, gene_list_all, gene_list_appr)

# 格式4: 适用于HRR，山东齐鲁 HRR全血、组织和配对，包含体细胞I/II/III类变异和胚系3/4/5类变异
# 包含基因、检测内容、丰度和变异解读，每个变异一行
# summary基因中分为已获批基因（展示时需要标星号）和其他基因
def format4(var_data, mlpa_data):
	gene_list_all = ["BRCA1", "BRCA2", "AR", "ATM", "ATR", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "ESR1", \
					 "FANCA", "FANCL", "HDAC2", "HOXB13", "MRE11A", "NBN", "PALB2", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", \
					 "RAD51D", "RAD54L", "STK11", "TP53", "BRAF", "ERBB2", "KRAS", "NRAS", "PIK3CA"]
	gene_list_appr = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L"]
	return getSum_forGermlinePanel(var_data, mlpa_data, gene_list_all, gene_list_appr)

# 函数3：适用于116系列产品（CRC25组织/组织配对/血液配对、TC21组织/组织配对/血液配对、GA18组织/组织配对/血液配对）和宁波病理
# 展示体细胞/来源不明 I/II/III/肿瘤发生发展相关变异，胚系45类变异
def getSum_forSomaticPanel(var_data, gene_list):
	var_data_copy = copy.deepcopy(var_data)
	result = []
	var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [4,5]) or \
				(var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5])]
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_list:
				# 致病/疑似致病但无用药的归为III类
				var["clinic_num_s"] = S_level(var)
				result.append(var)
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	return sorted(result, key = lambda i:i["gene_symbol"])

#格式5：适用于116系列产品（CRC25组织/组织配对/血液配对、TC21组织/组织配对/血液配对、GA18组织/组织配对/血液配对）和宁波病理
def format5(var_data):
	crc25_list = ["AKT1", "ARAF", "BRAF", "EGFR", "EPCAM", "ERBB2", "FBXW7", "GNAS", "HRAS", "KDR", "KRAS", \
				  "MAP2K1", "MET", "MLH1", "MSH2", "MSH6", "NRAS", "NTRK1", "PIK3CA", "PMS2", "POLE", "PTEN", "SMAD4"]
	tc21_list = ["AKT1", "ALK", "BRAF", "CTNNB1", "EIF1AX", "GNAS", "HRAS", "KRAS", "NRAS", "NTRK1", "NTRK3", \
				 "PAX8", "PDGFRA", "PIK3CA", "PTEN", "RASAL1", "RET", "TERT", "TP53", "TSC2", "TSHR"]
	ga18_list = ["AKT1", "BRAF", "EGFR", "ERBB2", "FGF19", "FGFR1", "FGFR2", "FGFR3", "HRAS", "KIT", "KRAS", "MET", \
				 "MTHFR", "NRAS", "PDGFRA", "PIK3CA", "PTEN", "TP53"]
	lc10_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "RET", "ROS1", "NRAS", "PIK3CA"]
	CRC25_result = getSum_forSomaticPanel(var_data, crc25_list)
	TC21_result = getSum_forSomaticPanel(var_data, tc21_list)
	GA18_result = getSum_forSomaticPanel(var_data, ga18_list)
	LC10_result = getSum_forSomaticPanel(var_data, lc10_list)
	return CRC25_result, TC21_result, GA18_result, LC10_result

#格式6：适用于云南肿瘤LC10
#ROS1/RET/ALK是I/II类融合，MET跳跃突变是MET I类变异，其余匹配hgvs_p
def format6(var_data):
	sv_gene = ["ALK", "ROS1", "RET"]
	spec_gene = ["MET"]
	snvindel_gene = {
		"BRAF" : ["V600E"],
		"EGFR" : ["E746_A750del", "L747_T751del", "E746_S752delinsV", "L747_P753delinsS", "L747_A750delinsP", \
				  "D770_N771insG", "L858R", "L861Q","G719A","T790M"],
		"ERBB2" : ["Y772_A775dup"],
		"KRAS" : ["G12S", "G12D", "G12A", "G12V", "G12C"],
		"NRAS" : ["G12D", "Q61R"],
		"PIK3CA" : ["H1047R"]
	}
	result = {}
	# 处理融合变异, I/II类，（肿瘤发生发展归为III类，此处不展示），注意主基因可能会为“gene1,gene2”
	for gene in sv_gene:
		if gene not in result.keys():
			result.setdefault(gene+"_sv", [])
		result[gene+"_sv"] = [var["var_hgvs"] for var in var_data if var["clinic_num_s"] in [4,5] and var["var_origin"] != "germline" and \
							  "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
							  set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and \
							  var["bio_category"] == "Sv" and gene in re.split(",", var["gene_symbol"])]	
	# 处理MET 14跳跃突变，I类
	result.setdefault("MET_14", [])
	result["MET_14"] = [(var["hgvs_c"],var["hgvs_p"]) for var in var_data if var["clinic_num_s"] == 5 and \
						var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
						set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and \
						var["bio_category"] == "Snvindel" and var["gene_symbol"] == "MET"]
	# 处理其他点突变/插入缺失突变
	for gene, info in snvindel_gene.items():
		for hgvs_p in info:
			if gene+"_"+hgvs_p not in result.keys():
				result.setdefault(gene+"_"+hgvs_p, [])
			result[gene+"_"+hgvs_p] = [var["hgvs_p"] for var in var_data if var["gene_symbol"] == gene and \
									   var["bio_category"] == "Snvindel" and \
									   hgvs_p == var["hgvs_p_ZJZL"].replace("p.", "")]

	return result

#格式7：适用于云南肿瘤HRR
# 胚系4/5类，BRCA MLPA del在模板里判定就好，这边只考虑snvindel
# 体细胞I/II类
def format7(var_data):
	result = {}
	gene_list = ["AR","ATM","BRCA1","BRCA2","CDK12","CHEK1","CHEK2","ERBB2","FANCA","FANCL","KRAS","MRE11A",\
				 "NRAS","PALB2","PIK3CA","PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53"]
	for gene in gene_list:
		if gene not in result.keys():
			result.setdefault(gene, [])
		result[gene] = [var["hgvs_c"] for var in var_data if ((var["clinic_num_s"] in [4,5] and \
						var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
						set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys())) or \
						(var["clinic_num_g"] in [4,5] and var["var_origin"] == "germline")) and \
						gene in re.split(",", var["gene_symbol"]) and var["bio_category"] == "Snvindel"]
	return result

#格式8：适用于云南肿瘤HRD C
# 仅展示I/II类变异
# BRCA1/BRCA2需要区分具体外显子号，其他基因不用
# BRCA Exon2-3：需展示exon2 + intron2 + exon3的变异
# BRCA Exon4：仅展示exon4的变异
def format8(var_data):
	gene_region = ["BRCA1", "BRCA2"]
	other_gene = ["ATM","BARD1","BRIP1","CDH1","CDK12","CHEK1","CHEK2","FANCA","FANCL","HDAC2","PALB2","PPP2R2A",\
				  "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","TP53"]
	var_data_copy = copy.deepcopy(var_data)
	var_list = [var for var in var_data_copy if var["clinic_num_s"] in [4,5] and \
				var["var_origin"] != "germline" and \
				"evi_sum" in var.keys() and \
				var["evi_sum"]["evi_split"] and \
				set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and \
				var["bio_category"] == "Snvindel"]
	result = {}
	# 处理BRCA变异
	for var in var_list:
		if var["gene_symbol"] in gene_region:
			for region in re.split("_", var["gene_region"]):
				region = region.replace("5'UTR", "5UTR").replace("3'UTR", "3UTR")
				if var["gene_symbol"]+"_"+region not in result.keys():
					result.setdefault(var["gene_symbol"]+"_"+region, [])
				result[var["gene_symbol"]+"_"+region].append((var["hgvs_c"], var["hgvs_p"]))
	# 处理其他基因变异
	for gene in other_gene:
		if gene not in result.keys():
			result.setdefault(gene, [])
		result[gene] = [(var["hgvs_c"],var["hgvs_p"]) for var in var_list if var["gene_symbol"] == gene]
	return result

#格式9：适用于福建肿瘤LC10结果解读部分-2022.11.03
#第一段：展示基因介绍
#第二段：配置表证据介绍
#第三段：变异解读
##注意：前三段若变异发生在相同基因上，则第一二段写一遍，第三段合并写。
#第四段：药物总结，所有变异靶向信息都汇总到最后一段
def format9(var_data):
	var_level = s_var_rule(var_data)
	### 处理前三段内容
	# 根据基因将变异进行分类处理
	var_gene_group = {}
	for var in var_level["level_I"] + var_level["level_II"] + var_level["level_onco_nodrug"] + var_level["level_III"]:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in var_gene_group.keys():
				var_gene_group.setdefault(gene, [])
			var_gene_group[gene].append(var)
	### 处理变异解读
	def stran_var_inter(i):
		inter = ""
		if i["bio_category"] == "Snvindel":
			if i["hgvs_p"] != "p.?":
				inter = i["gene_symbol"]+"基因的"+i["hgvs_c"]+" "+i["hgvs_p"]+i["variant_desc_cn"].replace("该变异","")+str(i["variant_interpret_cn"])
			else:
				inter = i["gene_symbol"]+"基因的"+i["hgvs_c"]+i["variant_desc_cn"].replace("该变异","")+str(i["variant_interpret_cn"])
		elif i["bio_category"] == "Cnv":
			inter = "本次实验检测到"+i["gene_symbol"]+"扩增，"+str(i["variant_interpret_cn"])
		elif i["bio_category"] == "Sv":
			if i["variant_interpret_cn"]:
				inter = "本次实验检测到"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"融合，融合类型为"+i["five_prime_gene"]+":"+i["five_prime_cds"]+\
						"-"+i["three_prime_gene"]+":"+i["three_prime_cds"]+"。"+i["variant_interpret_cn"]
			else:
				inter = "本次实验检测到"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"融合，融合类型为"+i["five_prime_gene"]+":"+i["five_prime_cds"]+\
						"-"+i["three_prime_gene"]+":"+i["three_prime_cds"]+"。"
		return inter
	
	### 处理基因介绍
	def stran_gene_inter(i):
		inter = {}
		if not re.search(",", i["gene_symbol"]):
			inter[i["gene_symbol"]] = i["gene_function"]
		else:
			if i["bio_category"] == "Sv" or i["bio_category"] == "PSeqRnaSv":
				inter[i["five_prime_gene"]] = i["five_prime_gene_function"]
				inter[i["three_prime_gene"]] = i["three_prime_gene_function"]
		
		return inter

	# 前三段信息整合
	result = {}
	result["var"] = []
	for gene, var in var_gene_group.items():
		var_info = []
		gene_info = {}
		rpt_info = {}
		for i in var:
			var_info.append(stran_var_inter(i))
			gene_info.update(stran_gene_inter(i))
			if i["var_info_forFJZL"]:
				rpt_info.update(i["var_info_forFJZL"])
		result["var"].append(
			{
				"gene_info" : gene_info.get(gene),
				"rpt_info" : rpt_info.get(gene),
				"var_info" : "".join(var_info)
			}
		)

	## 处理最后一段内容
	# 处理药物总结中的变异描述
	def stran_var_for_drug(i):
		inter = ""
		if i["bio_category"] == "Snvindel":
			type_cn = i["type_cn"] if i["type_cn"] != "--" else i["type"]
			type_cn = type_cn.replace("Intronic", "内含子突变").replace("3'UTR","3'UTR突变").replace("5'UTR", "5'UTR突变")
			if i["hgvs_p"] != "p.?":
				inter = "检出"+i["gene_symbol"]+"基因的"+i["hgvs_c"]+" "+i["hgvs_p"]+type_cn
			else:
				inter = "检出"+i["gene_symbol"]+"基因的"+i["hgvs_c"]+type_cn
		elif i["bio_category"] == "Cnv":
			inter = "检出"+i["gene_symbol"]+"扩增"
		elif i["bio_category"] == "Sv":
			inter = "检出"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"融合"

		return inter

	# 汇总药物总结
	result["drug"] = ""
	regimen_sum = []
	for var in var_level["level_I"] + var_level["level_II"]:
		var_regimen_S = []
		var_regimen_R = []
		var_regimen = []
		for i in ["A","B","C","D"]:
			if var["var_regimen_forFJZL"]["regimen_"+str(i)+"_S"]:
				var_regimen_S.append("、".join(var["var_regimen_forFJZL"]["regimen_"+str(i)+"_S"])+"（证据等级为"+str(i)+"）")
			if var["var_regimen_forFJZL"]["regimen_"+str(i)+"_R"]:
				var_regimen_R.append("、".join(var["var_regimen_forFJZL"]["regimen_"+str(i)+"_R"])+"（证据等级为"+str(i)+"）")

		if var_regimen_S:
			var_regimen.append("可能对"+"，".join(var_regimen_S)+"敏感")
		if var_regimen_R:
			var_regimen.append("可能对"+"，".join(var_regimen_R)+"耐药")
		regimen_sum.append(stran_var_for_drug(var)+"，该突变提示"+"；".join(var_regimen)+"。")
	result["drug"] = "".join(regimen_sum)
	return result

#格式10：适用于云南肿瘤CP40/CRC12/TC21/GA18
# CP40展示I/II/III和肿瘤发生发展相关变异
# CRC12展示I/II类变异，且需要详细变异型，如C12
# TC21/GA18展示I/II类变异（不包含肿瘤发生发展相关）
def format10(var_data):
	var_data_copy = copy.deepcopy(var_data)
	CP40_gene = ["AKT1","ALK","BRAF","CDK4","CTNNB1","DDR2","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","KIT","KRAS","MAP2K1",\
				 "MET","MYC","NRAS","NTRK1","NTRK2","NTRK3","PDGFRA","PIK3CA","PTEN","POLE","RB1","RET","ROS1","TP53"]
	TC21_gene = ["BRAF","KRAS","NRAS","HRAS","GNAS","TERT","RET","ALK","NTRK1","PAX8","PIK3CA","AKT1","EIF1AX","TSC2",\
				 "CTNNB1","PDGFRA","RASAL1","TP53","PTEN","TSHR","NTRK3"]
	GA18_gene = ["AKT1","BRAF","EGFR","ERBB2","FGF19","FGFR1","FGFR2","FGFR3","HRAS","KIT","KRAS","MET","MTHFR",\
				 "NRAS","PDGFRA","PIK3CA","PTEN","TP53"]
	CP40_var_list = [var for var in var_data_copy if var["var_origin"] != "germline" and var["clinic_num_g"] in [3,4,5]]
	Other_var_list = [var for var in var_data_copy if (var["var_origin"] == "germline" and var["clinic_num_g"] in [4,5]) or\
					 (var["var_origin"] != "germline" and var["clinic_num_s"] in [4,5] and \
					  "evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and \
					  set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))]
	# 适用CP40、TC21、GA18
	def summary(var_list, gene_list):
		result = {}
		for gene in gene_list:
			if gene not in result.keys():
				result.setdefault(gene, [])
			result[gene] = [var for var in var_list if gene in re.split(",", var["gene_symbol"])]
		return result

	CP40_result = summary(CP40_var_list, CP40_gene)
	TC21_result = summary(Other_var_list, TC21_gene)
	GA18_result = summary(Other_var_list, GA18_gene)

	# CRC12-完全匹配氨基酸变化
	CRC12_snvindel_complete = {
		"BRAF" : ["V600E"],
		"POLE" : ["P286R", "V411L", "S297F"]
	}
	# CRC12-匹配氨基酸变化开头
	CRC12_snvindel_simple = {
		"KRAS" : ["G12", "G13", "Q61", "A146"],
		"NRAS" : ["G12", "G13", "Q61", "A146"],
		"PIK3CA" : ["E542", "E545", "H1047"]
		}
	# CRC12-扩增，做成列表形式，方便后续新增
	CRC12_cnv = ["ERBB2"]
	# CRC12-融合
	CRC12_sv = ["NTRK1", "NTRK3", "ALK", "RET", "FGFR3"]

	def summary_CRC12(var_list):
		result = {}
		# 处理融合变异, I/II类，注意主基因可能会为“gene1,gene2”
		for gene in CRC12_sv:
			if gene not in result.keys():
				result.setdefault(gene+"_sv", [])
			result[gene+"_sv"] = [var["var_hgvs"] for var in var_list if var["bio_category"] == "Sv" and gene in re.split(",", var["gene_symbol"])]
		# 处理扩增
		for gene in CRC12_cnv:
			if gene not in result.keys():
				result.setdefault(gene+"_cnv", [])
			result[gene+"_cnv"] = ["扩增" for var in var_list if var["bio_category"] == "Cnv" and gene == var["gene_symbol"]]
		# 处理snvindel， 完全匹配hgvs_p
		for gene, hgvs_p_list in CRC12_snvindel_complete.items():
			for hgvs_p in hgvs_p_list:
				if gene+"_"+hgvs_p not in result.keys():
					result.setdefault(gene+"_"+hgvs_p, [])
				result[gene+"_"+hgvs_p] = [var["hgvs_p"] for var in var_list if var["bio_category"] == "Snvindel" and \
										   var["gene_symbol"] == gene and hgvs_p == (var["hgvs_p_ZJZL"].replace("p.", ""))]
		# 处理snvindel，不完全匹配hgvs_p
		for gene, hgvs_p_list in CRC12_snvindel_simple.items():
			for hgvs_p in hgvs_p_list:
				if gene+"_"+hgvs_p not in result.keys():
					result.setdefault(gene+"_"+hgvs_p, [])
				result[gene+"_"+hgvs_p] = [var["hgvs_p"] for var in var_list if var["bio_category"] == "Snvindel" and \
										   var["gene_symbol"] == gene and re.search(hgvs_p, var["hgvs_p"])]
		return result
	CRC12_result = summary_CRC12(Other_var_list)

	return CP40_result, TC21_result, GA18_result, CRC12_result

# 格式11：适用于复旦中山CP40，不在其他字典里添加了，怕出问题。-2022.11.28
# 1. 基因检测列表
# 2. 总结：检测到XX、XX基因突变。
def format11(var_data):
	var_data_copy = copy.deepcopy(var_data)
	cp40_list = ["AKT1","ALK","BRAF","CDK4","CTNNB1","DDR2","EGFR","ERBB2","ESR1","FGFR1",\
				"FGFR2","FGFR3","FGFR4","HRAS","IDH1","IDH2","KEAP1","KIT","KRAS","MAP2K1",\
				"MET","MYC","NFE2L2","NKX2-1","NRAS","NRG1","NTRK1","NTRK2","NTRK3","PDGFRA",\
				"PIK3CA","POLE","PTEN","RB1","RET","ROS1","STK11","TP53"]
	result = []
	detect_gene = []
	var_list = [var for var in var_data_copy if var["var_origin"] != "germline" and var["clinic_num_s"] in [3,4,5]]
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in cp40_list:
				detect_gene.append(gene)
				var["clinic_num_s"] = S_level(var)
				var["gene_symbol"] = gene
				result.append(var)
	for gene in set(cp40_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})	
	return sorted(result, key = lambda i:i["gene_symbol"]), "、".join(sorted(set(detect_gene)))

# 格式12：适用于广东人民gHRR，展示胚系5/4/3类变异
def format12(var_data, mlpa_data):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_copy = copy.deepcopy(mlpa_data)
	gene_list = ["BRCA1","BRCA2","AR","ATM","ATR","BARD1","BRIP1","CDH1","CDK12","CHEK1",\
				 "CHEK2","ESR1","FANCA","FANCL","HDAC2","HOXB13","MRE11","NBN","PALB2","PPP2R2A",\
				 "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53","BRAF","ERBB2","KRAS",\
				 "NRAS","PIK3CA"]
	var_list = [var for var in var_data_copy if var["var_origin"] == "germline" \
												and var["clinic_num_g"] in [3,4,5] \
												and var["gene_symbol"] in gene_list \
												and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for k, v in mlpa_data_copy.items():
		if v:
			var_list += v
			detect_gene+= [var["gene_symbol"] for var in v]

	for gene in set(gene_list) - set(detect_gene):
		var_list.append({
			"gene_symbol" : gene
		})
	
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))

#格式13：适用于云南肿瘤BRCA
# 仅展示I/II/III类和胚系3/4/5类变异
# BRCA1/BRCA2需要区分具体外显子号
# BRCA Exon2-3：需展示exon2 + intron2 + exon3的变异
# BRCA Exon4：仅展示exon4的变异
def format13(var_data, mlpa_data):
	var_data_copy = copy.deepcopy(var_data)
	var_list = [var for var in var_data_copy if ((var["clinic_num_s"] in [3,4,5] and var["var_origin"] != "germline") or \
												(var["clinic_num_g"] in [3,4,5] and var["var_origin"] == "germline")) and \
												var["bio_category"] == "Snvindel" and var["gene_symbol"] in ["BRCA1", "BRCA2"]]
	result = {}
	for var in var_list:
		for region in re.split("_", var["gene_region"]):
			region = region.replace("5'UTR", "5UTR").replace("3'UTR", "3UTR")
			if var["gene_symbol"]+"_"+region not in result.keys():
				result.setdefault(var["gene_symbol"]+"_"+region, [])
			result[var["gene_symbol"]+"_"+region].append((var["hgvs_c"], var["hgvs_p"]))
	
	# MLPA结果
	# value中的格式为：exon1a, exon1b, exon1-exon3, exon5-exon16
	mlpa_data_copy = copy.deepcopy(mlpa_data["B1_Loss"]+mlpa_data["B2_Loss"]+mlpa_data["B1_Gain"]+mlpa_data["B2_Gain"])
	for var in mlpa_data_copy:
		for i in re.split(",", var["value"].replace(" ","")):
			# 处理为[1, 1, 1-3, 5-16]
			i = i.replace("exon","").replace("a","").replace("b", "").replace("up", "")
			exon_num = re.split("-", i)
			# 涉及外显子只有1个时
			if len(exon_num) == 1:
				if var["gene_symbol"]+"_exon"+str(exon_num[0]) not in result.keys():
					result.setdefault(var["gene_symbol"]+"_exon"+str(exon_num[0]), [])
				result[var["gene_symbol"]+"_exon"+str(exon_num[0])].append("mlpa")
			# 涉及多个外显子时
			elif len(exon_num) == 2:
				for num in range(int(exon_num[0]), int(exon_num[1])+1):
					if var["gene_symbol"]+"_exon"+str(num) not in result.keys():
						result.setdefault(var["gene_symbol"]+"_exon"+str(num), [])
					result[var["gene_symbol"]+"_exon"+str(num)].append("mlpa")
			else:
				pass

	return result

#格式14：适用于云南肿瘤BPTM（仅处理POLE和TP53，BRCA可以用格式13）
# 仅展示I/II/III类POLE和TP53
# 需要区分具体外显子号
# BRCA Exon2-3：需展示exon2 + intron2 + exon3的变异
# BRCA Exon4：仅展示exon4的变异
def format14(var_data):
	var_data_copy = copy.deepcopy(var_data)
	var_list = [var for var in var_data_copy if ((var["clinic_num_s"] in [3,4,5] and var["var_origin"] != "germline") or \
												(var["clinic_num_g"] in [3,4,5] and var["var_origin"] == "germline")) and \
												var["bio_category"] == "Snvindel" and var["gene_symbol"] in ["POLE", "TP53"]]
	result = {}
	for var in var_list:
		for region in re.split("_", var["gene_region"]):
			region = region.replace("5'UTR", "5UTR").replace("3'UTR", "3UTR")
			if var["gene_symbol"]+"_"+region not in result.keys():
				result.setdefault(var["gene_symbol"]+"_"+region, [])
			result[var["gene_symbol"]+"_"+region].append((var["hgvs_c"], var["hgvs_p"]))
	
	return result