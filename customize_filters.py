#-*- coding:gbk -*-
import jinja2
import re
from docxtpl import RichText
import jinja2.filters
from libs.rule import nofound_genelist
from libs.rule import decimal_float

'''
自定义jinja2过滤器
'''

### 定制类 ###
# 1. 复旦中山厦门医院CP40，药物“、”.join()展示为一行
def regimen_sum(a):
	result = []
	evi_sum = a["evi_sum"]
	if "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				result.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				result.append("{0}（{1}，{2}级）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
	return "、".join(result) if result else "-"
jinja2.filters.FILTERS["regimen_sum_filter"] = regimen_sum 

# 2. 北大人民检验科HRR-出BRCA报告-检测结论
# gene_symbol gene_region hgvs_c hgvs_p, clinic_num_g.strans
def var_brca_sum(var_list):
	result = []
	clinic_trans = {5 : "致病性变异", 4 : "疑似致病性变异"}
	for var in var_list:
		if var["type"] == "Loss":
			result.append("{0} {1} del, {2}".format(var["gene_symbol"], var["value"], "疑似致病性变异"))

		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}, {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], clinic_trans[var["clinic_num_g"]]))
			else:
				result.append("{0} {1} {2}, {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], clinic_trans[var["clinic_num_g"]]))
	return "；".join(result)
jinja2.filters.FILTERS["var_brca_sum"] = var_brca_sum

# 3. 北大人民检验科HRR-出BRCA报告-变异解读结论
def inter_brca_sum(var_list):
	result_dict = {}
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			var["clinic_num_g"] = 4
		elif var["type"] == "Gain":
			var["clinic_num_g"] = 3

		if var["clinic_num_g"] not in result_dict.keys():
			result_dict.setdefault(var["clinic_num_g"], [])
		result_dict[var["clinic_num_g"]].append(var)
	clinic_trans = {5:"致病性变异", 4:"疑似致病性变异", 3:"意义不明确变异", 2:"疑似良性变异"}
	for i in [5,4,3,2]:
		if i in result_dict.keys():
			result.append(str(len(result_dict[i]))+"个"+clinic_trans[i])
	return "、".join(result)
jinja2.filters.FILTERS["inter_brca_sum"] = inter_brca_sum

# 孙逸仙116-结果小结
def summary_SYX_116(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+" 扩增")
		elif var["bio_category"] == "Sv":
			result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
	return ", ".join(result)
jinja2.filters.FILTERS["summary_SYX_116"] = summary_SYX_116

# 中山人民HRDC检测结果小结-2023.06.02
def summary_ZSRM_hrd(info):
	# 阳性：基于本次送检样本，检出BRCA1基因BRCA1 gene_region hgvs_c hgvs_p I类变异，HRD状态阳性以及TP53基因TP53 gene_region hgvs_c hgvs_p II类变异。
	# 阴性：基于本次送检样本，未检出致病性或可能致病性变异。
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I类", 4 : "II类"}
	
	brca_result = []
	for var in var_dict["ec_type"]["BRCA1_level12"] + var_dict["ec_type"]["BRCA2_level12"]:
		if var["hgvs_p"] != "p.?":
			brca_result.append("{0}基因{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			brca_result.append("{0}基因{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	brca_result_str = "，".join(brca_result)
	
	hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"

	other_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["gene_symbol"] not in ["BRCA1", "BRCA2"]:
			if var["hgvs_p"] != "p.?":
				other_result.append("{0}基因{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
			else:
				other_result.append("{0}基因{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	other_result_str = ",".join(other_result)
	hrd_result_str = "以及".join([hrd_result, other_result_str])

	result = []
	if brca_result_str:
		result.append(brca_result_str)
	if hrd_result_str:
		result.append(hrd_result_str)

	# 组装-阳性/阴性
	if var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		return "基于本次送检样本，检出{0}。".format("，".join(result))
	else:
		return "基于本次送检样本，未检出致病性或疑似致病性变异。"
jinja2.filters.FILTERS["summary_ZSRM_hrd"] = summary_ZSRM_hrd

##############################################################################################################

### 通用类 ###

### 变异列表类 ###
# 1. 体细胞：返回肿瘤发生发展相关变异列表
def return_somatic_onconodrug_list(info):
	var = info[0]
	sample = info[1]
	if re.search("Master", sample["prod_names"]):
		if not sample["control_sample_id"]:
			return var["var_somatic"]["level_onco_nodrug"] + var["var_germline_nodrug"]
		else:
			return var["var_somatic"]["level_onco_nodrug"]
	else:
		return var["var_somatic"]["level_onco_nodrug"]
jinja2.filters.FILTERS["somatic_onconodrug_list"] = return_somatic_onconodrug_list

# 2. 体细胞：返回体细胞/未知来源变异，III类变异列表
def return_somatic_level3_list(info):
	var = info[0]
	sample = info[1]
	if sample["prod_names"] == "BPTM（5基因）":
		return var["ec_type"]["POLE_level3"]+var["ec_type"]["TP53_level3"]+var["ec_type"]["BRCA1_level3"]+var["ec_type"]["BRCA2_level3"]
	elif sample["prod_names"] == "PTM（3基因）":
		return var["ec_type"]["POLE_level3"]+var["ec_type"]["TP53_level3"]
	else:
		return var["var_somatic"]["level_III"]
jinja2.filters.FILTERS["somatic_level3_list"] = return_somatic_level3_list

# 3. 返回林奇变异列表，未检测到变异的基因也要展示
def return_gLS5_list(gLS5):
	result = []
	gene_list = ["EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	for gene in gene_list:
		if gene in gLS5.keys() and gLS5[gene]:
			result.extend(gLS5[gene])
		else:
			result.append({"gene_symbol" : gene})
	return result
jinja2.filters.FILTERS["gLS5_list"] = return_gLS5_list
	
# 4. 返回BRCA基因检测结果，未检测到变异的基因也要展示
def return_brca_list(brca):
	result = []
	brca1_list = brca["snv_s"]["B1_L5"] + brca["snv_s"]["B1_L4"]+brca["mlpa"]["B1_Loss"]+brca["mlpa"]["B1_Gain"]+brca["snv_s"]["B1_L3"]
	brca2_list = brca["snv_s"]["B2_L5"] + brca["snv_s"]["B2_L4"]+brca["mlpa"]["B2_Loss"]+brca["mlpa"]["B2_Gain"]+brca["snv_s"]["B2_L3"]
	if brca1_list:
		result.extend(brca1_list)
	else:
		result.append({"gene_symbol" : "BRCA1"})
	
	if brca2_list:
		result.extend(brca2_list)
	else:
		result.append({"gene_symbol" : "BRCA2"})
	return result
jinja2.filters.FILTERS["brca_list"] = return_brca_list
	

#----------------------------------------------------------------------------------------------------------------

### 变异解读类 ###
# 1. 体细胞：返回体细胞/未知来源变异 I、II类解读变异列表
def return_somatic_level12_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		# var_origin != "germline"
		brca_result = []
		for var in var_brca["snv_s"]["B1_L5"] + var_brca["snv_s"]["B2_L5"] + var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B2_L4"]:
			if var["var_origin"] != "germline":
				brca_result.append(var)
		return brca_result
	elif re.search("BPTM（5基因）", sample["prod_names"]):
		return var["ec_type"]["POLE_level12"] + var["ec_type"]["TP53_level12"] + var["ec_type"]["BRCA1_level12"] + var["ec_type"]["BRCA2_level12"]
	elif re.search("PTM（3基因）", sample["prod_names"]):
		return var["ec_type"]["POLE_level12"] + var["ec_type"]["TP53_level12"]
	# 如果是MP且无配对时，var_origin为germline的且有用药的也要展示
	elif re.search("Master Panel（组织）", sample["prod_names"]) and not sample["control_sample_id"]:
		return var["var_somatic"]["level_I"] + var["var_germline"]["regimen_level_I"] + var["var_somatic"]["level_II"] + var["var_germline"]["regimen_level_II"]
	else:
		return var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"]
jinja2.filters.FILTERS["somatic_level12_inter"] = return_somatic_level12_inter

# 2. 胚系：返回胚系变异4/5类列表，用于解读
def return_germline_level45_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		return var_brca["snv_m"]["B1_G_L5"]+var_brca["snv_m"]["B2_G_L5"]+var_brca["snv_m"]["B1_G_L4"]+var_brca["snv_m"]["B2_G_L4"]+\
			   var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]
	elif re.search("HRR", sample["prod_names"]):
		return var["var_germline"]["level_5"]+var["var_germline"]["level_4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"]
	else:
		return var["var_germline"]["level_5"]+var["var_germline"]["level_4"]
jinja2.filters.FILTERS["germline_level45_inter"] = return_germline_level45_inter

#3. 胚系：返回胚系变异3类列表，用于解读
def return_germline_level3_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		return var_brca["mlpa"]["B1_Gain"]+var_brca["mlpa"]["B2_Gain"]+var_brca["snv_m"]["B1_G_L3"]+var_brca["snv_m"]["B2_G_L3"]
	elif re.search("HRR", sample["prod_names"]):
		return var_brca["mlpa"]["B1_Gain"]+var_brca["mlpa"]["B2_Gain"]+var["var_germline"]["level_3"]
	else:
		return var["var_germline"]["level_3"]
jinja2.filters.FILTERS["germline_level3_inter"] = return_germline_level3_inter

#----------------------------------------------------------------------------------------------------------------

### 变异结果中的内容 ###
# 1. 返回基因，单个或多个
def return_gene_symbol(var):
	result = []
	if "bio_category" in var.keys() and var["bio_category"]:
		if var["bio_category"] in ["Snvindel", "Cnv"]:
			result.append(var["gene_symbol"])
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if "," in var["gene_symbol"]:
				if var["five_prime_gene"] != var["three_prime_gene"]:
					result.append(var["five_prime_gene"])
					result.append(var["three_prime_gene"])
				else:
					result.append(var["five_prime_gene"])
			else:
				result.append(var["gene_symbol"])
		else:
			result.append(var["gene_symbol"])
	else:
		# gene_symbol 可能不存在，这边加个兼容-2023.09.22
		if "gene_symbol" in var.keys():
			result.append(var["gene_symbol"])
		# 兼容完成-2023.09.22
	rt = RichText()
	rt.add("\n".join(result), italic=True)
	return rt
jinja2.filters.FILTERS["gene_symbol"] = return_gene_symbol

# 2. 变异检测结果
def var_info(var):
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
			result.append(var["transcript_primary"])
		elif var["bio_category"] == "Cnv":
			if var['cnv_type'] == 'loss' or var['cnv_type'] == 'Loss':
				result.append("缺失")
			else:
				result.append("扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
				result.append(var["five_prime_transcript"])
			else:
				result.append("{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"]))
				result.append("{0}/{1}".format(var["five_prime_transcript"], var["three_prime_transcript"]))
		elif var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["value"]+" dup")
	rt = RichText()
	rt.add("\n".join(result))
	return rt
	#return result
jinja2.filters.FILTERS["var_info"] = var_info


# 3. 返回变异丰度，适用BRCA等，胚系展示纯合/杂合，体细胞展示频率，MLPA用-
def freq_stran(info):
	'''
	体细胞snvindel：freq_str
	胚系snvindel：厦门项目freq， 上海项目freq_rc，注意MP还要考虑是否有对照样本，无对照样本，来源为germline的直接展示freq_str
	sv：CP40 copies，其他freq_str
	PSeqRnaSv：freq
	'''
	var = info[0]
	sample = info[1]
	if "var_origin" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				if "Master" in sample["prod_names"]:
					if sample["control_sample_id"]:
						return "纯合" if float(var["freq_rc"]) >= 0.85 else "杂合"
					else:
						return var["freq_str"]
					pass 
				elif re.search("116|76|25|21|18", sample["prod_names"]):
					print (var["freq_rc"])
					return "纯合" if var["freq_rc"] and float(var["freq_rc"]) >= 0.85 else "杂合" if var["freq_rc"] and float(var["freq_rc"]) < 0.85 else "未提取到freq_rc！"
				else:
					return "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "Cnv":
			return var["cn_mean"]
		elif var["bio_category"] == "Sv":
			# 新增OncoPro（组织） XW1402
			if re.search("Classic|CRC12|OncoPro（组织）", sample["prod_names"]):
				return str(var["copies"])+" copies"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "PSeqRnaSv":
			return str(var["freq"])+" copies"
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "/"
jinja2.filters.FILTERS["freq_stran"] = freq_stran

# 4. 临床意义-药物
def significance_regimen(var):
	result = []
	#rt = RichText()
	evi_sum = var["evi_sum"]
	if evi_sum and "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				result.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				result.append("{0}（{1}，{2}级）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
		#rt.add("\n".join(result)) 
	if not result:
		result = ["-"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
	#return result
jinja2.filters.FILTERS["significance_regimen"] = significance_regimen

# 5. 返回临床意义，适用BRCA等
def clinic_stran(info):
	var = info[0]
	s_type = info[1]
	g_dict = {
		5 : "致病性",
		4 : "疑似致病性",
		3 : "意义不明确",
		2 : "疑似良性",
		1 : "良性"
	}
	s_dict1 = {
		5 : "I类-强临床意义",
		4 : "II类-潜在临床意义",
		3 : "III类-临床意义不明",
		2 : "IV类-良性/可能良性",
		1 : "IV类-良性/可能良性"
	}
	s_dict2 = {
		5 : "I类",
		4 : "II类",
		3 : "III类",
		2 : "IV类",
		1 : "IV类"
	}
	# 药企输出内容
	s_dict3 = {
		5 : "I类变异",
		4 : "II类变异",
		3 : "III类变异",
		2 : "IV类变异",
		1 : "IV类变异"
	}
	s_dict = s_dict1 if s_type == "s1" else s_dict2 if s_type == "s2" else s_dict3 if s_type == "s3" else {}
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			return g_dict.get(var["clinic_num_g"], "")
		else:
			# 肿瘤发生发展相关归为III类
			if var["clinic_num_s"] in [4, 5] and (not var["evi_sum"]["evi_split"] or \
				(var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))) :
				return s_dict.get(var[3])
			else:
				return s_dict.get(var["clinic_num_s"], "")
	else:
		if var["type"] == "Loss":
			return g_dict.get(4)
		elif var["type"] == "Gain":
			return g_dict.get(3)
jinja2.filters.FILTERS["clinic_stran"] = clinic_stran

# 5.1 返回临床意义，适用MP
def clinic_stran_MP(info):
	var = info[0]
	s_type = info[1]
	sample = info[2]
	g_dict = {
		5 : "致病性",
		4 : "疑似致病性",
		3 : "意义不明确",
		2 : "疑似良性",
		1 : "良性"
	}
	s_dict1 = {
		5 : "I类-强临床意义",
		4 : "II类-潜在临床意义",
		3 : "III类-临床意义不明",
		2 : "IV类-良性/可能良性",
		1 : "IV类-良性/可能良性"
	}
	s_dict2 = {
		5 : "I类",
		4 : "II类",
		3 : "III类",
		2 : "IV类",
		1 : "IV类"
	}
	s_dict = s_dict1 if s_type == "s1" else s_dict2 if s_type == "s2" else {}
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			if sample["control_sample_id"]:
				return g_dict.get(var["clinic_num_g"], "")
			else:
				return s_dict.get(var["clinic_num_s"], "")
		else:
			# 肿瘤发生发展相关归为III类
			if var["clinic_num_s"] in [4, 5] and (not var["evi_sum"]["evi_split"] or \
				(var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))) :
				return s_dict.get(var[3])
			else:
				return s_dict.get(var["clinic_num_s"], "")
	else:
		if var["type"] == "Loss":
			return g_dict.get(4)
		elif var["type"] == "Gain":
			return g_dict.get(3)
jinja2.filters.FILTERS["clinic_stran_MP"] = clinic_stran_MP


# 6. 变异来源转化
def var_origin_stran(var):
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			return "胚系"
		elif var["var_origin"] == "somatic":
			return "体细胞"
		else:
			return "待定"
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "胚系"
jinja2.filters.FILTERS["var_origin_stran"] = var_origin_stran

# 7. 结果小结表格表头-丰度-适用于BRCA和林奇
def title_freq_brca(prod_names):
	if prod_names in ["BRCA1/BRCA2（全血）", "BRCA1/2（扩增子）", "HRR（全血）", "林奇综合征"]:
		return  "基因型"
	elif prod_names in ["BRCA1/BRCA2（组织）"]:
		return "丰度"
	elif prod_names in ["BRCA1/BRCA2（组织 全血）"]:
		return "基因型/丰度"
jinja2.filters.FILTERS["title_freq_brca"] = title_freq_brca

# 8. III 类变异表头-丰度-体细胞
def title_freq_III(sample):
	if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/2（扩增子）", "HRR（全血）", "林奇综合征"]:
		return  "基因型"
	elif sample["prod_names"] in ["BRCA1/BRCA2（组织）","BPTM（5基因）", "PTM（3基因）", "HRR（组织）", "HRD Complete（组织）", "10基因（血液）", \
								  "Pan116（血液）", "LC76（血液）", "CRC25（血液）", "TC21（血液）", "GA18（血液）"]:
		return "丰度"
	elif sample["prod_names"] in ["BRCA1/BRCA2（组织 全血）", "HRR（组织 全血）"]:
		return "丰度"
	elif sample["prod_names"] in ["10基因（组织）", "Classic Panel", "CRC12-MSI",\
								  "Pan116（组织）", "LC76（组织）", "CRC25（组织）", "TC21（组织）", "GA18（组织）",\
								  "Master Panel（组织）", "Master Panel（血液）"]:
		return "丰度/拷贝数"
jinja2.filters.FILTERS["title_freq_III"] = title_freq_III

# 9. 肿瘤发生发展相关变异-丰度
# 目前涉及的项目有HRR组织/HRR配对/HRD/LC10组织/LC10血液/PAN116组织/PAN116血液/LC76组织/LC76血液/CRC25组织/CRC25血液/TC21组织/TC21血液/GA18组织/GA18血液
# Master组织/Master血液/CP40/CRC12
# 丰度/拷贝数
def title_freq_onconodrug(sample):
	if sample["prod_names"] in ["HRR（组织）", "HRR（组织 全血）", "HRD Complete（组织）", "10基因（血液）", \
							   "Pan116（血液）", "LC76（血液）", "CRC25（血液）", "TC21（血液）", "GA18（血液）"]:
		return "丰度"
	elif sample["prod_names"] in ["Pan116（组织）", "LC76（组织）", "CRC25（组织）", "TC21（组织）", "GA18（组织）",\
								 "Master Panel（组织）", "Master Panel（血液）", "Classic Panel", "CRC12-MSI", "10基因（组织）"]:
		return "丰度/拷贝数"
jinja2.filters.FILTERS["title_freq_onconodrug"] = title_freq_onconodrug

# 10. 靶向治疗表-丰度
# 适用于MP/116/LC10/CP40/CRC12
def title_freq_targetRegimen(sample):
	freq = ["丰度"]
	if re.search("组织", sample["prod_names"]):
		freq.append("拷贝数")
	if sample["prod_names"] in ["Classic Panel", "CRC12-MSI"]:
		freq.append("拷贝数")
	if sample["prod_names"] == "Master Panel（血液）":
		freq.append("拷贝数")
	if sample["control_sample_id"]:
		freq.append("基因型")
	return "/".join(freq)
jinja2.filters.FILTERS["title_freq_targetRegimen"] = title_freq_targetRegimen

# 11 PARP结果汇总表-丰度表头
def title_parp(sample):
	if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "HRR（全血）"]:
		return "基因型"
	elif sample["prod_names"] in ["HRR（组织）", "HRD Complete（组织）"]:
		return "丰度"
	elif sample["prod_names"] in ["HRR（组织 全血）"]:
		return "基因型/丰度"
jinja2.filters.FILTERS["title_parp"] = title_parp
#----------------------------------------------------------------------------------------------------------------

### 结果小结类 ###
# 1. PD-L1
def pdl1_summary(pdl1):
	return "阴性。" if pdl1["result"] == "阴性" else "阳性，{0}为{1}。".format(pdl1["type"], pdl1["value"])
jinja2.filters.FILTERS["pdl1_summary"] = pdl1_summary

# 2. MSI
def msi_summary(msi):
	return "微卫星稳定（MSS）。" if msi["var_id"] == "MSS" else "微卫星不稳定（MSI-H）。" if msi["var_id"] == "MSI-H" else "未返回MSI结果！"
jinja2.filters.FILTERS["msi_summary"] = msi_summary

# 3. TMB
def tmb_summary(tmb):
	TMB_result = "低" if tmb["var_id"] == "TMB-L" else "高"
	return "{0} Muts/Mb, 肿瘤突变负荷较{1}（{2}）。".format(tmb["TMB_value"], TMB_result, "TMB-H" if tmb["var_id"] == "TMB-H" else "TMB-L")
jinja2.filters.FILTERS["tmb_summary"] = tmb_summary

# 4. GEP
def gep_summary(gep):
	return "GEP分值为{0}分".format(gep["gep_score"])
jinja2.filters.FILTERS["gep_summary"] = gep_summary

# 5. TME
def tme_summary(tme):
	tme_dict = {
		"IE/F" : "免疫富集/纤维化亚型(IE/F)",
		"IE" : "免疫富集/非纤维化亚型(IE)",
		"F" : "纤维化亚型(F)",
		"D" : "免疫荒漠型(D)"
	}
	return "TME分型为{0}".format(tme_dict.get(tme["tme_type"], tme["tme_type"]))
jinja2.filters.FILTERS["tme_summary"] = tme_summary

# 6. 免疫检查点抑制剂疗效相关基因
def io_summary(io):
	result = []
	if io["io_p_summary"]:
		result.append(io["io_p_summary"]+"（疗效正相关）")
	if io["io_n_summary"]:
		result.append(io["io_n_summary"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	rt = RichText()
	rt.add("；\n".join(result)+"。")
	return rt
jinja2.filters.FILTERS["io_summary"] = io_summary

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

# 7. HRD检测结果
def hrd_summary(info):
	gss = info[0]
	sample = info[1]
	result = []
	hrd_result = "HRD阳性" if gss["BRCA1"] or gss["BRCA2"] or float(gss["gss"]["gsscore"]) >= 50 else "HRD阴性"
	# 加结果仅供参考的一些情况-2023.08.29
	# 未检测到BRCA致病/疑似致病变异+以下任意一个条件
	# 1）baf_noise > 0.055
	# 2) depth_noise > 0.35
	# 3) 无肿瘤细胞含量
	# 4）有肿瘤细胞含量但是格式不正确的
	# 5）有肿瘤细胞含量，格式正确，但是数值<30
	note = "（结果仅供参考）" if not gss["BRCA1"] and not gss["BRCA2"] and \
							(float(gss["gss"]["baf_noise"]) > 0.055 or \
							 float(gss["gss"]["depth_noise"]) > 0.35 or \
							 not sample["tumor_content"] or \
							 (sample["tumor_content"] and not is_number(sample["tumor_content_num"])) or \
							 (sample["tumor_content"] and is_number(sample["tumor_content_num"]) and float(sample["tumor_content_num"]) < 30)) \
							else ""
	result.append(hrd_result + note)
	if gss["summary"]:
		result.append("HRR通路相关基因突变：{0}".format(gss["summary"]))
	rt = RichText()
	rt.add("；\n".join(result)+"。")
	return rt
jinja2.filters.FILTERS["hrd_summary"] = hrd_summary

# 8. GA检测结果
def ga_summary(info):
	ga = info[0]
	msi = info[1]
	result = []
	if ga["ebv_type"]["ebv_type"] == "P" and msi["var_id"] == "MSI-H":
		if ga["ebv_sum"]:
			result.append("EB病毒感染型（EBV），{0}".fromat(ga["ebv_sum"]))
		else:
			result.append("EB病毒感染型（EBV）")
		result.append("微卫星不稳定型（MSI）")
	elif ga["ebv_type"]["ebv_type"] == "P":
		result.append("EB病毒感染型（EBV）")
	elif msi["var_id"] == "MSI-H":
		result.append("微卫星不稳定型（MSI）")
	else:
		if ga["gs_sum"]:
			result.append("基因组稳定型（GS）")
		if ga["cin_sum"]:
			result.append("染色体不稳定型（CIN）")
	if not result:
		result = ["-"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["ga_summary"] = ga_summary

# 9. 子宫内膜癌分子分型结果
def ec_summary(ec_type):
	ec_dict = {
		"POLE-ultramutated type EC" : "POLE突变型（POLE mutation，POLE mut）",
		"MSI-H type EC" : "错配修复功能缺陷（Mismatch repair deficiency，MMRd）",
		"CNH type EC" : "TP53基因突变（p53 abnormality，p53 abn）",
		"CNL type EC" : "非特异性分子谱（Non-specific molecular profile，NSMP）"
	}
	return ec_dict.get(ec_type, ec_type)
jinja2.filters.FILTERS["ec_summary"] = ec_summary
	
# 10. 胚系变异（需要展示具体变异）
def var_g_summary(var_list):
	g_var_list = var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"]
	result = []
	if g_var_list:
		for var in g_var_list:
			hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p"]
			clinic_g_cn = "致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"
			result.append("检出{0} {1}，为{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	else:
		result = ["在检测范围内，未检出致病性/疑似致病性变异"]
	return "；".join(result)+"。"
jinja2.filters.FILTERS["var_g_summary"] = var_g_summary

# 10. 体细胞/来源不明变异
def var_s_summary(info):
	var_list = info[0]
	sample = info[1]
	result = []
	def sum_var(vlist):
		v_result = []
		for var in vlist:
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					v_result.append(var["gene_symbol"]+" "+var["hgvs_p"])
				else:
					v_result.append(var["gene_symbol"]+" "+var["hgvs_c"])
			elif var["bio_category"] == "Cnv":
				v_result.append(var["gene_symbol"]+" 扩增")
			elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					v_result.append("MET exon14 跳跃")
				else:
					if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in v_result:
						v_result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
		return ", ".join(v_result)
	# 有对照样本时
	c_var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"])
	c_var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	c_var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"])
	c_var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	# 无对照样本时
	var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"] +\
				 var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"])
	var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
							var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_germline_nodrug"])
	var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
								  var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	if sample["control_sample_id"]:
		str_1 = "检出{0}个体细胞变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(c_var_all_num), str(c_var_onco_drug_num), str(c_var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(c_var_onco_drug_str)
		if c_var_all_num != 0:
			if c_var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出体细胞变异。"
	else:
		str_1 = "检出{0}个变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(var_all_num), str(var_onco_drug_num), str(var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(var_onco_drug_str)
		if var_all_num != 0:
			if var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出变异。"
jinja2.filters.FILTERS["var_s_summary"] = var_s_summary

#----------------------------------------------------------------------------------------------------------------

### HRD ###
# 1. BRCA检测结果
def hrd_brca_result(var_list):
	result = []
	if var_list:
		for var in var_list:
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
	else:
		result = ["未检出致病性或疑似致病性变异"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["hrd_brca_result"] = hrd_brca_result

### 治疗策略-解读 ###
def get_evi_sum(evi_sum):
	result = []
	if evi_sum and "evi_split" in evi_sum.keys():
		if "Predictive_merge" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Predictive_merge"]:
			result.extend([{"regimen" : a["regimen_name"], "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Predictive_merge"]])
		if "Prognostic" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Prognostic"]:
			result.extend([{"regimen" : "预后相关", "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Prognostic"]])
		if "Diagnostic" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Diagnostic"]:
			result.extend([{"regimen" : "辅助诊断相关", "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Diagnostic"]])
	if not result:
		result = [{"regimen" : "", "inter": "目前关于该变异的临床治疗实践尚不明确。"}]
	return result
jinja2.filters.FILTERS["evi_sum"] = get_evi_sum

### 其他类 ###
# 1. 治疗方案介绍-生物标志物
def approval_regimen_biomarker(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], i["hgvs_p"]))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF 野生型")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD阴性")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD阳性")
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	#print (result_redup)
	#print (judge_mergeMET)
	if judge_mergeMET:
		if "MET exon14 跳跃" in result_redup:
			result_redup.remove("MET exon14 跳跃")
	if not result_redup:
		result_redup = ["-"]
	rt = RichText()
	rt.add("\n".join(result_redup))
	return rt
	#return result_redup
jinja2.filters.FILTERS["approval_regimen_biomarker"] = approval_regimen_biomarker

# 2. PARP结果汇总表，根据产品（HRR和HRD）选择对应展示列表
def choose_parp_list(info):
	hrr_list = info[0]["cdx"]["format5_forHRR"]
	hrd_list = info[0]["cdx"]["format4_forHRDC"]
	sample = info[1]
	if sample["prod_names"] in ["HRR（全血）", "HRR（组织）", "HRR（组织 全血）"]:
		return hrr_list
	elif sample["prod_names"] in ["HRD Complete（组织）"]:
		return hrd_list
jinja2.filters.FILTERS["choose_parp_list"] = choose_parp_list

# 体细胞等级转化
def somatic_class_stran(clinic_num_s):
	if clinic_num_s == 5:
		return "I类"
	elif clinic_num_s == 4:
		return "II类"
jinja2.filters.FILTERS["somatic_class_stran"] = somatic_class_stran


# 5. 检测详细结果-检测结果
def detect_result(var):
	result = []
	if var["bio_category"] == "Snvindel":
		if var["hgvs_p"] != "p.?":
			result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
		else:
			result.append(var["gene_region"]+" "+var["hgvs_c"])
		result.append(var["transcript_primary"])
	elif var["bio_category"] == "Cnv":
		result.append("扩增")
	elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
		if var["five_prime_gene"] == var["three_prime_gene"] and var["three_prime_gene"] == "MET":
			result.append("MET exon14 跳跃")
			result.append(var["three_prime_transcript"])
		else:
			result.append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
			result.append(var["five_prime_transcript"]+"/"+var["three_prime_transcript"])
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["detect_result"] = detect_result

# 6. 检测详细结果-频率-体细胞
def detect_freq_s(var):
	freq = ""
	if var["bio_category"] == "Snvindel":
		freq = var["freq_str"]
	elif var["bio_category"] == "Cnv":
		freq = var["cn_mean"]
	elif var["bio_category"] == "Sv":
		freq = var["freq_str"]
	elif var["bio_category"] == "PSeqRnaSv":
		freq = str(var["freq"])+" copies"
	return freq
jinja2.filters.FILTERS["freq_s"] = detect_freq_s

# 7. 检测详细结果-频率-胚系
# 适用上海项目，有对照样本的情况
def detect_freq_g_SH(var):
	return "未提取到数据！" if not var["freq_sc"] else "纯合" if float(var["freq_sc"]) >= 0.85 else "杂合"
jinja2.filters.FILTERS["freq_g_SH"] = detect_freq_g_SH

# io展示整理为方法
def io_stran(io_list):
	if not io_list:
		io_list = ["-"]
	rt = RichText()
	rt.add("\n".join(io_list))
	return rt
jinja2.filters.FILTERS["io_stran"] = io_stran

# NCCN指南推荐基因检测结果，检测结果展示-MP
def cdx_type1(info):
	cdx_list = info[0]
	sample = info[1]
	rt = RichText()
	result = []
	if cdx_list:
		for var in cdx_list:
			if var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					result.append("MET exon14 跳跃，{0}".format(var["freq"]))
				else:
					if "var_info_M" in var.keys() and var["var_info_M"]:
						result.append(var["var_info_M"]+"，"+var["freq"])
					else:
						result.append(var["var_info"]+"，"+var["freq"])
			else:
				if var["gene_symbol"] in ["BRCA1", "BRCA2"]:
					if sample["control_sample_id"]:
						if var["var_origin"] == "germline":
							result.append("{0}，{1}".format(var["var_info"], var["freq_rc_str"]))
							result.append("胚系变异")
						else:
							result.append("{0}，{1}".format(var["var_info"], var["freq"]))
							result.append("体细胞变异")
					else:
						result.append("{0}，{1}".format(var["var_info"], var["freq"]))
				else:
					result.append("{0}，{1}".format(var["var_info"], var["freq"]))
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type1"] = cdx_type1

# NCCN指南推荐基因检测结果，检测结果展示-116变异检测结果
def cdx_type2_var_info(var):
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				return var["hgvs_p"]
			else:
				return var["hgvs_c"]
		elif var["bio_category"] == "Cnv":
			return "扩增"
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				return "MET exon14 跳跃"
			else:
				return "{0}-{1} 融合".format(var["five_prime_gene"], var["three_prime_gene"])
	else:
		return "未检测到"
jinja2.filters.FILTERS["cdx_type2_var_info"] = cdx_type2_var_info

# NCCN指南推荐基因检测结果，检测结果展示-116丰度
def cdx_type2_freq(info):
	var = info[0]
	sample = info[1]
	rt = RichText()
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				result.append("未提取到freq_rc！" if not var["freq_rc"] else "纯合" if float(var["freq_rc"]) >= 0.85 else "杂合")
				#rt.add("未提取到freq_rc！" if not var["freq_rc"] else "纯合" if float(var["freq_rc"]) >= 0.85 else "杂合", \
				#	   size = 17, color = "#000000", font = "Source Han Sans Normal")
			else:
				result.append(var["freq_str"])
				#rt.add(var["freq_str"], size = 17, color = "#000000", font = "Source Han Sans Normal")
			if var["gene_symbol"] in ["BRCA1", "BRCA2"] and sample["control_sample_id"]:
				result.append("胚系变异" if var["var_origin"] == "germline" else "体细胞变异")
				#rt.add("\n胚系变异" if var["var_origin"] == "germline" else "\n体细胞变异", size = 17, color = "#000000", font = "Source Han Sans Normal")
		elif var["bio_category"] == "Cnv":
			result.append(var["cn_mean"])
			#rt.add(var["cn_mean"], size = 17, color = "#000000", font = "Source Han Sans Normal")
		elif var["bio_category"] == "Sv":
			result.append(var["freq_str"])
			#rt.add(var["freq_str"], size = 17, color = "#000000", font = "Source Han Sans Normal")
	else:
		result.append("-")
		#rt.add("-", size = 17, color = "#000000", font = "Source Han Sans Normal")
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type2_freq"] = cdx_type2_freq

# NCCN指南推荐基因检测结果，检测结果展示-116临床意义
def cdx_type2_class(var):
	stran = {
		"germline" : {
			5 : "致病性",
			4 : "疑似致病性"
		},
		"somatic" : {
			5 : "I类",
			4 : "II类",
			3 : "III类"
		}
	}
	if "bio_category" in var.keys():
		if var["var_origin"] == "germline":
			return stran["germline"].get(var["clinic_num_g"])
		else:
			return stran["somatic"].get(var["clinic_num_s"])
	else:
		return "-"
jinja2.filters.FILTERS["cdx_type2_class"] = cdx_type2_class

# NCCN指南推荐基因检测结果-CP40
def cdx_type3_var_info(var):
	result = []
	if "var_info" in var.keys() and var["var_info"]:
		if "MET-MET" in var["var_info"]:
			result.append("MET exon14 skipping")
		elif var["merge_sv_list"]:
			for a in var["merge_sv_list"]:
				result.append(a+"融合")
		else:
			result.append(var["var_info"])
		if var["note"]:
			result.append("(MET exon14 skipping)")
	else:
		result = ["未检测到"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type3_var_info"] = cdx_type3_var_info

# NCCN指南推荐基因检测结果-丰度
def cdx_type3_freq(var):
	result = []
	if "freq" in var.keys() and var["freq"]:
		if "融合" in var["var_info"]:
			result.append(str(var["freq"])+" copies")
		else:
			result.append(var["freq"])
		if var["note_freq"]:
			result.append("("+var["note_freq"]+")")
	else:
		result.append("-")
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type3_freq"] = cdx_type3_freq

# NCCN指南推荐基因检测结果-临床意义
def cdx_type3_class(var):
	if "level" in var.keys() and var["level"]:
		if var["level"] == "5":
			return "I类"
		elif var["level"] == "4":
			return "II类"
		else:
			return "III类"
	else:
		return "-"
jinja2.filters.FILTERS["cdx_type3_class"] = cdx_type3_class

# 选择CRC25、TC21、GA18检测结果汇总
def choose_116_list(info):
	CRC25_list = info[0]["CRC25"]
	GA18_list = info[0]["GA18"]
	TC21_list = info[0]["TC21"]
	prod_names = info[1]
	if prod_names in ["CRC25（组织）", "CRC25（血液）"]:
		return CRC25_list
	elif prod_names in ["GA18（组织）", "GA18（血液）"]:
		return GA18_list
	elif prod_names in ["TC21（组织）", "TC21（血液）"]:
		return TC21_list
jinja2.filters.FILTERS["choose_116_list"] = choose_116_list

### 产品声明 ###
# 1. 序号
def product_state_index(info):
	var_brca = info[0]
	sample = info[1]
	if sample["prod_names"] in ["BRCA1/BRCA2（全血）",  "BRCA1/2（扩增子）"]:
		if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"] or\
																	   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
			if re.search("健康", sample["tumor_names_cn"]):
				return "5"
			else:
				return "7"
		else:
			if re.search("健康", sample["tumor_names_cn"]):
				return "4"
			else:
				return "6"
	elif sample["prod_names"] in ["BRCA1/BRCA2（组织）"]:
		if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"]:
			return "7"
		else:
			return "6"
	elif sample["prod_names"] in ["BRCA1/BRCA2（组织 全血）"]:
		if var_brca["snv_m"]["B1_G_L5"] or var_brca["snv_m"]["B2_G_L5"] or var_brca["snv_m"]["B1_G_L4"] or var_brca["snv_m"]["B2_G_L4"] or\
																	   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
			return "7"
		else:
			return "6"
	elif sample["prod_names"] in ["林奇综合征"]:
		return "4"
	# 加一个HRR全血，健康人4，患者6
	elif sample["prod_names"] in ["HRR（全血）"]:
		if re.search("健康", sample["tumor_names_cn"]):
			return "4"
		else:
			return "6"
	else:
		return "6"
jinja2.filters.FILTERS["product_state_index"] = product_state_index

### sample ###
# 1. 样本类型
# 116血液项目，样本类型会将样本和对照的写一起
def sample_type(sample):
	result = ""
	if sample["prod_names"] in ["Pan116（血液）", "TC21（血液）", "GA18（血液）", "LC76（血液）", "CRC25（血液）", "Master Panel（血液）"] and sample["control_sample_id"]:
		if sample["sample_type"] and sample["control_sample_type"]:
			if sample["sample_type"] != sample["control_sample_type"]:
				result = sample["sample_type"]+","+sample["control_sample_type"]
			else:
				result = sample["sample_type"]
		elif sample["sample_type"] and not sample["control_sample_type"]:
			result = sample["sample_type"]
		elif not sample["sample_type"] and sample["control_sample_type"]:
			result = sample["control_sample_type"]
	else:
		result = sample["sample_type"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["sample_type"] = sample_type

# 2. 样本数量
# 116血液项目，样本数量会将样本和对照的写一起
def sample_amount(sample):
	result = ""
	if sample["prod_names"] in ["Pan116（血液）", "TC21（血液）", "GA18（血液）", "LC76（血液）", "CRC25（血液）", "Master Panel（血液）"] and sample["control_sample_id"]:
		if sample["sample_type"] and sample["control_sample_type"]:
			if sample["sample_type"] != sample["control_sample_type"]:
				result = sample["sample_amount"]+","+sample["control_sample_amount"]
			else:
				result = sample["sample_amount"]
		elif sample["sample_type"] and not sample["control_sample_type"]:
			result = sample["sample_amount"]
		elif not sample["sample_type"] and sample["control_sample_type"]:
			result = sample["control_sample_amount"]
	else:
		result = sample["sample_amount"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["sample_amount"] = sample_amount

# 3. 样本采集日期
# 血液项目采集日期从blood_collection_date中提取，其他从gather_data中提取
def gather_data(sample):
	result = ""
	if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "HRR（全血）", "10基因（血液）", "61遗传基因", "Pan116（血液）", "林奇综合征",\
							    "TC21（血液）", "GA18（血液）", "LC76（血液）", "CRC25（血液）", "BPTM（全血）", "Master Panel（血液）", "BRCA1/2（扩增子）", "MRD（LC10）"]:
		result = sample["blood_collection_date"]
	else:
		result = sample["gather_data"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["gather_data"] = gather_data


### 治疗方案介绍-过滤 ###
# MP组织，仅卵巢癌、乳腺癌、前列腺癌时展示HRD解读，故其他癌种在治疗方案部分，需要过滤掉HRD+和HRD-的内容
def approval_regimen_filter_hrd(info):
	regimen_list = info[0]
	sample = info[1]
	hrd_p = {"biomarker_type" : "HRD+"}
	hrd_n = {"biomarker_type" : "HRD-"}
	result = []
	if "Master" in sample["prod_names"] and not set(["乳腺癌", "卵巢癌", "前列腺癌"])&set(sample["tumor_list"]) and "实体瘤" not in sample["tumor_names_cn"]:
		for regimen in regimen_list:
			for a in regimen["var"]:
				if a == hrd_p or a == hrd_n:
					regimen["var"].remove(a)
			if regimen["var"]:
				result.append(regimen)
	else:
		result = regimen_list
	return result
jinja2.filters.FILTERS["approval_regimen_filter_hrd"] = approval_regimen_filter_hrd


### QC ###
# 1. HRR全血
def qc_gHRR(qc, lib_quality_control):
	ngs_qc = qc["dna_data_qc"] if "dna_data_qc" in qc.keys() else {}
	lib_qc = lib_quality_control["lib_dna_qc"] if "lib_dna_qc" in lib_quality_control.keys() else {}
	if float(ngs_qc["cleandata_q30_num"]) >= 0.75 and float(ngs_qc["depth_ssbc_num"]) >= 100:
		if lib_qc and "dna_qty" in lib_qc.keys() and lib_qc["dna_qty"] and "library_qty" in lib_qc.keys() and lib_qc["library_qty"]:
			if float(lib_qc["dna_qty"]) >= 20 and float(lib_qc["library_qty"]) >= 200:
				return "合格"
			else:
				return "风险"
		else:
			return "合格（质控项有缺失，请补齐数据后自行评估）"
	else:
		return "风险"
jinja2.filters.FILTERS["qc_gHRR"] = qc_gHRR

# 2. HRR组织

qc_stand = {
	"gHRR" : {},
	"tHRR" : {
		"qc" : {
			"dna_data_qc" : {
				"cleandata_q30_num" : 0.75,
				"depth_ssbc_num" : 300
			}
		},
		"lib_quality_control" : {
			"lib_dna_qc" : {
				"dna_qty" : 30,
				"library_qty" : 200
			}
		},
		"tumor_content" : 20
	}
}

def qc_tHRR(qc, lib_quality_control, sample, prod):
	stand = qc_stand.get(prod)
	# qc结果
	qc_result = []
	if "dna_data_qc" in qc.keys() and qc["dna_data_qc"]:
		for item in stand["qc"]["dna_data_qc"].keys():
			if item in qc["dna_data_qc"].keys() and float(qc["dna_data_qc"][item]) >= stand["dna_data_qc"][item]:
				qc_result.append("合格")
			else:
				qc_result.append("风险")
	# 湿实验结果
	if "lib_dna_qc" in lib_quality_control.keys() and lib_quality_control["lib_dna_qc"]:
		for item in stand["lib_quality_control"]["dna_data_qc"].keys():
			pass

# v4临检通用模板使用
# MP小结展示具体变异-体细胞变异
def var_sum_s(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+" 扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			# 融合可能会有重复（rna exon相同，断点不同的情况）
			if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in result:
				result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
	return ", ".join(result)
jinja2.filters.FILTERS["var_sum_s_filter"] = var_sum_s

# MP小结展示具体变异-胚系变异
def var_sum_g(var_list):
	result = []
	for var in var_list:
		hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p"]
		clinic_g_cn = "致病变异" if var["clinic_num_g"] == 5 else "疑似致病变异"
		result.append("检出{0} {1}，为{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	return "；".join(result)
jinja2.filters.FILTERS["var_sum_g_filter"] = var_sum_g


######################################## 定制模板过滤器#################################################
# 重庆西南CP40，检测结果中要列出【研究证据】参考文献
# 规则如下：
# 1. A或C3等级，如果匹配到NCCN的，展示“NCCN临床实践指南”
# 2. 有PMID的展示PMID号，没有的展示NCT号，都没有的展示斜杠，都有的展示PMID号
# 3. NCCN和（PMID or NCT）可能会同时出现，NCCN要展示，（PMID or NCT）放其中1个
def get_refer_CQXN(evi_sum):
	# 获取PMID
	def getPMID_from_inter(inter):
		pmid_list = []
		mat = re.compile(r"PMID.\s?\d+")
		for i in mat.findall(str(inter)):
			if re.search(":|: |：|： ", i):
				pmid = (re.split(":|: |：|： ", i))[1].replace(" ", "")
			else:
				pmid = (re.split("PMID", i))[1]
			pmid_list.append(pmid)
		return pmid_list
	# 获取NCT
	def getNCT_from_inter(inter):
		nct_list = []
		mat = re.compile(r"NCT\d+")
		for i in mat.findall(str(inter)):
			nct_list.append(i)
		return nct_list
	# 汇总研究证据，包含治疗、辅助诊断和预后
	evi_list = []
	for i in ["Predictive", "Prognostic", "Diagnostic"]:
		if i in evi_sum["evi_split"].keys() and evi_sum["evi_split"][i]:
			evi_list.extend(evi_sum["evi_split"][i])
	# 提取参考文献
	refer = []
	for evi in evi_list:
		# 判断是否有nccn
		if "A" in evi["evi_conclusion"] or "C3" in evi["evi_conclusion"]:
			if re.search("NCCN", evi["evi_interpretation"]):
				refer.append("NCCN临床实践指南")
		# 判断pmid和NCT
		pmid_list = getPMID_from_inter(evi["evi_interpretation"])
		nct_list = getNCT_from_inter(evi["evi_interpretation"])
		if pmid_list:
			for i in pmid_list:
				refer.append("PMID: "+str(i))
		else:
			if nct_list:
				refer.extend(nct_list)
	# 去重
	refer_redup = []
	for i in refer:
		if i not in refer_redup:
			refer_redup.append(i)
	if not refer_redup:
		refer_redup = ["/"]
	rt = RichText()
	rt.add("\n".join(refer_redup))
	return rt
jinja2.filters.FILTERS["refer_CQXN"] = get_refer_CQXN

def ZDY_HRR_germline_summary(var_list):
	var_dict = {}
	for var in var_list:
		if var["type"] == "Loss":
			var_info = "{0} {1} del".format(var["gene_symbol"], var["value"])
			var["clinic_num_g"] = 4
		else:
			if var["hgvs_p"] != "p.?":
				var_info = "{0} {1}:{2}:{3}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"], var["hgvs_p"])
			else:
				var_info = "{0} {1}:{2}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"])
		if var["clinic_num_g"] not in var_dict.keys():
			var_dict.setdefault(var["clinic_num_g"], [])
		var_dict[var["clinic_num_g"]].append(var_info)
	result = []
	if 5 in var_dict.keys():
		result.append("检出{0}个致病性变异，{1}".format(str(len(var_dict[5])), ", ".join(var_dict[5])))
	if 4 in var_dict.keys():
		result.append("检出{0}个疑似致病性变异，{1}".format(str(len(var_dict[4])), ", ".join(var_dict[4])))
	if not var_dict:
		result.append("未检出致病或疑似致病性胚系变异")
	return "；".join(result)
jinja2.filters.FILTERS["ZDY_HRR_germline_summary"] = ZDY_HRR_germline_summary

# 中山人民HRDC检测结果小结-更新-2023.07.12
def summary_ZSRM_hrd_v2(info):
	# 阳性：基于本次送检样本，检出HRD状态【阳性/阴性】、BRCA gene_region hgvs_c hgvs_p I类变异、TP53 gene_region hgvs_c hgvs_p。
	# 阴性：基于本次送检样本，检出HRD【阳性/阴性】，未检出I类、II类变异。
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I类", 4 : "II类"}

	var_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["hgvs_p"] != "p.?":
			var_result.append("{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			var_result.append("{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	
	hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"

	if var_result:
		return "基于本次送检样本，检出"+hrd_result+"、"+"、".join(var_result)+"。"
	else:
		return "基于本次送检样本，检出"+hrd_result+"，未检出I类、II类变异。"
jinja2.filters.FILTERS["summary_ZSRM_hrd_v2"] = summary_ZSRM_hrd_v2

# 重庆西南-三类变异删除解读为3的非编码区变异（剪接变异除外）
def filter_III_var_CQXN(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] in ["Sv", "Cnv"]:
			result.append(var)
		elif var["bio_category"] == "Snvindel":
			if var["clinic_num_s"] == 3 and var["type"] in ["Intronic", "3'UTR", "5'UTR", "FlankingRegion3", "FlankingRegion5"]:
				pass
			else:
				result.append(var)
	return result
jinja2.filters.FILTERS["filter_III_var_CQXN"] = filter_III_var_CQXN

# 深圳二院HRR，可能获益临床试验仅展示BRCA基因-2023.07.17
def filter_clinicalTrial_SZEY(clinical_list):
	return [a for a in clinical_list if a["gene_symbol"] in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["filter_clinicalTrial_SZEY"] = filter_clinicalTrial_SZEY

# 深圳二院HRR，药物介绍仅展示BRCA基因相关的-2023.07.17
def filter_drug_SZEY(drug_list):
	for drug in drug_list:
		drug["var_filter"] = [a for a in drug["var"] if ("gene_symbol" in a.keys() and a["gene_symbol"] and a["gene_symbol"] in ["BRCA1", "BRCA2"]) or\
							  "biomarker_type" in a.keys() and "BRCA1" in a["biomarker_type"] or \
							  "biomarker_type" in a.keys() and "BRCA2" in a["biomarker_type"]]
	return [a for a in drug_list if a["var_filter"]]
jinja2.filters.FILTERS["filter_drug_SZEY"] = filter_drug_SZEY

# 河北省人民CP40，I类变异仅展示AB药物，II类变异仅展示CD药物
# info格式为[msi, knb, var_level_I, var_level_II, regimen_approval_list]
def filter_drug_HNRM(info):
	# 处理变异列表，返回指定列表、指定等级药物的结果，格式为[("阿法替尼", var1), ("阿法替尼", var2), ("厄洛替尼", var1)]
	def get_list(var_list, level_list):
		result = []
		for var in var_list:
			bio = []
			if "bio_category" in var.keys():
				if var["bio_category"] == "Snvindel":
					bio.append({
						"bio_category" : var["bio_category"],
						"gene_symbol" : var["gene_symbol"],
						"hgvs_c" : var["hgvs_c"],
						"hgvs_p" : var["hgvs_p"]
					})
				elif var["bio_category"] == "Cnv":
					bio.append({"bio_category" : var["bio_category"], "gene_symbol" : var["gene_symbol"]})
				elif var["bio_category"] == "Sv":
					if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
						bio.append({"bio_category" : var["bio_category"], "hgvs" : "MET exon14 skipping"})
					else:
						for i in var["merge_sv_list"]:
							bio.append({"bio_category" : var["bio_category"], "hgvs" : i+"融合"})
			elif "var_id" in var.keys():
				if var["var_id"] == "KRAS/NRAS/BRAF WT":
					bio.append({"var_id" : "KRAS/NRAS/BRAF 野生型"})
				elif var["var_id"] == "MSI-H":
					bio.append({"var_id" : "MSI-H"})
			if "evi_sum" in var.keys() and var["evi_sum"] and "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				for evi in var["evi_sum"]["evi_split"]["Predictive"]:
					if evi["evi_conclusion_simple"] in level_list:
						result.append((evi["regimen_name"], bio))
		return result
	
	# 1. 汇总MSI-H AB、KNB AB、I类AB和II类CD药物
	regimen_list_filter = []
	regimen_list_filter.extend(get_list([info[0]], ["A", "B"]))
	regimen_list_filter.extend(get_list([info[1]], ["A", "B"]))
	regimen_list_filter.extend(get_list(info[2], ["A", "B"]))
	regimen_list_filter.extend(get_list(info[3], ["C", "D"]))

	# 2. 格式转化{"阿法替尼" : [var1, var2], "厄洛替尼" : [var1]}
	regimen_dict_filter = {}
	for i in regimen_list_filter:
		regimen_dict_filter.setdefault(i[0], [])
		for j in i[1]:
			regimen_dict_filter[i[0]].append(j)

	# 对治疗方案介绍列表进行过滤，并新增过滤后的变异列表var_hnrm
	regimen_approval_list = info[4]
	regimen_approval_filter = []
	for regimen in regimen_approval_list:
		regimen_name = regimen["regimen_cn"] if regimen["regimen_cn"] else regimen["regimen_en"]
		if regimen_name in regimen_dict_filter.keys():
			regimen["var_hnrm"] = regimen_dict_filter[regimen_name]
			regimen_approval_filter.append(regimen)

	return regimen_approval_filter
jinja2.filters.FILTERS["filter_drug_HNRM"] = filter_drug_HNRM

# 杭州区域肺癌需要判断EGFR基因检出情况
def judge_EGFR(info):
	var_list = info[0]
	sample = info[1]
	gene_list = [var["gene_symbol"] for var in var_list]
	if "EGFR" not in gene_list and "肺癌" in sample["tumor_list"] and "实体瘤" not in sample["tumor_names_cn"]:
		return "EGFRnone"
jinja2.filters.FILTERS["judge_EGFR"] = judge_EGFR

# 杭州区域肺癌需要判断EGFR T790M突变情况
def sort_var_forhz_egfr(info):
	var_list = info[0]
	sample = info[1]
	for var in var_list:
		var["egfr_sort"] = 0
	t790m = {"gene_symbol" : "EGFR", "t790m_var_info" : "未检出T790M变异", "egfr_sort" : 1}
	egfr_hgvs_p_list = [var["hgvs_p"].replace("p.","").replace("(","").replace(")","") for var in var_list if var["gene_symbol"] == "EGFR" and var["bio_category"] == "Snvindel"]
	
	for var in var_list:
		if var["bio_category"] == "Cnv" and var["gene_symbol"] == "EGFR":
			egfr_hgvs_p_list.append("扩增")

	if "T790M" not in egfr_hgvs_p_list and egfr_hgvs_p_list and "肺癌" in sample["tumor_list"] and "实体瘤" not in sample["tumor_names_cn"]:
		# 计算EGFR最后一个变异所在index
		egfr_last_index = 0
		for var in var_list:
			if var["gene_symbol"] == "EGFR":
				egfr_last_index = var_list.index(var)
		# 插入T790M检测结果
		var_list.insert(egfr_last_index+1, t790m)

	return var_list
jinja2.filters.FILTERS["sort_var_forhz_egfr"] = sort_var_forhz_egfr

# 复旦中山CP40-小结部分展示化疗位点结果，格式上基因名需要合并单元格-2023.07.31
def chemo_stran_fdzs(chemo_list):
	result = {}
	for i in chemo_list:
		result.setdefault(i["gene_symbol"], [])
		result[i["gene_symbol"]].append({
			"dbsnp" : i["dbsnp"],
			"genotype" : i["genotype"]
		})
	chemo_result = []
	for k, v in result.items():
		v_sort = sorted(v, key=lambda i:i["dbsnp"])
		chemo_result.append({
			"gene_symbol" : k,
			"info" : v_sort
		})
	return sorted(chemo_result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["chemo_stran_fdzs"] = chemo_stran_fdzs

# 福建附一新增参考文献-2023.07.31
def refer_nccn_fjfy(sample):
	refer = {
		"非小细胞肺癌" : [
			"中国非小细胞肺癌患者EGFR+T790M基因突变检测专家共识（2018年版）",
			"二代测序技术在NSCLC中的临床应用中国专家共识（2020版）",
			"非小细胞肺癌分子病理检测临床实践指南（2021版）",
			"中国非小细胞肺癌RET基因融合临床检测专家共识（2021年版）",
			"非小细胞肺癌MET临床检测中国专家共识（2022年版）",
			"非小细胞肺癌细针穿刺细胞学标本基因检测专家共识（2022年版）",
			"非小细胞肺癌恶性浆膜腔积液分子病理检测中国专家共识（2022年版）",
			"非小细胞肺癌融合基因检测临床实践中国专家共识（2023年版）",
			"中国晚期非小细胞肺癌BRAF突变诊疗专家共识（2023年版）"
		],
		"结直肠癌" : [
			"结直肠癌及其他相关实体瘤微卫星不稳定性检测中国专家共识（2019年版）",
			"结直肠癌分子标志物临床检测中国专家共识（2021年版）",
			"结直肠癌分子检测高通量测序中国专家共识（2021年版）"
		],
		"卵巢癌" : [
			"卵巢上皮性癌BRCA基因检测的中国专家讨论（2017年版）",
			"基于下一代测序技术的BRCA1_2基因检测指南(2019年版)",
			"上皮性卵巢癌PARP抑制剂相关生物标志物检测的中国专家共识（2020年版）",
			"BRCA1_2数据解读中国专家共识(2021年版)"
		],
		"前列腺癌" : [
			"中国前列腺癌患者基因检测专家共识（2020年版）",
			"前列腺癌同源重组修复基因检测及变异解读专家共识（2022年版）"
		],
		"乳腺癌" : [
			"中国乳腺癌患者BRCA基因检测与临床应用专家共识（2018年版）",
			"基于下一代测序技术的BRCA1_2基因检测指南(2019年版)",
			"复发/转移性乳腺癌标志物临床应用专家共识（2019年版）",
			"晚期乳腺癌基因检测热点问题中国专家共识（2021年版）",
			"基于靶标指导乳腺癌精准治疗标志物临床应用专家共识（2022年版）"
		],
		"胃癌" : [
			"胃癌HER2检测指南（2016年版）",
			"HER2阳性晚期胃癌分子靶向治疗的中国专家共识（2016年版）",
			"胃癌高通量测序临床应用中国专家共识（2022年版）"
		]
	}
	result = []
	if "实体瘤" not in sample["tumor_names_cn"]:
		set_list = set(sample["tumor_list"]) & set([i for i in refer.keys()])
		for tumor in set_list:
			result.extend(refer[tumor])
	return result
jinja2.filters.FILTERS["refer_nccn_fjfy"] = refer_nccn_fjfy

# 删除变异中的BCL2L11变异-2023.08.02
def remove_bcl2l11(var_list):
	for var in var_list:
		if var["gene_symbol"] == "BCL2L11" and "hgvs_c" in var.keys() and var["hgvs_c"] and var["hgvs_c"] == "c.394+1479_394+4381del":
			var_list.remove(var)
	return var_list
jinja2.filters.FILTERS["remove_bcl2l11"] = remove_bcl2l11

# 西安交大一LC10/PAN116/CP40，删除KRAS、NRAS、BRAF、KNB中含有“CSCO”的证据
# 无需考虑删除证据后对变异等级的影响-2023.08.15
def filter_csco(evi_list):
	return [i for i in evi_list if "CSCO" not in i["evi_interpretation"]]
jinja2.filters.FILTERS["filter_csco"] = filter_csco

# 安徽省立BRCA-结果小结部分，需要基因斜体，并且变异之间用、间隔-2023.08.28
# 在每个变异（除了最后一个）后面加个顿号，就能在模板中使用for循环展示变异了
def ahsl_brca_sum(info):
	var_list = info[0]
	judge_intumor = info[1]
	var_g = {5 : "致病性变异", 4 : "疑似致病性变异"}
	for var in var_list:
		if var["type"] == "Loss":
				var["clinic_num_g_stran"] = "疑似致病性变异"
		else:
			if var["var_origin"] == "germline":
				var["clinic_num_g_stran"] = var_g.get(var["clinic_num_g"])
			else:
				if var["clinic_num_g"] in [4, 5]:
					var["clinic_num_s_stran"] = "I类变异" if judge_intumor == "intumor" else "II类变异"
				else:
					var["clinic_num_s_stran"] = "III类变异"
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			if var["type"] == "Loss" or var["var_origin"] == "germline":
				var["clinic_num_g_stran"] = var["clinic_num_g_stran"]+"、"
			else:
				var["clinic_num_s_stran"] = var["clinic_num_s_stran"]+"、"
	return var_list
jinja2.filters.FILTERS["ahsl_brca_sum"] = ahsl_brca_sum

# 南方医院CP40， NCCN指南推荐标志物检测结果汇总表，需要拆分为一/二/三级变异展示-2023.09.01
# 该过滤器可以返回指定等级结果
def split_cdx_nfyy(info):
	var_list = [var for var in info[0] if "var_info" in var.keys()]
	level = info[1]
	return [var for var in var_list if var["level"] == str(int(level))]
jinja2.filters.FILTERS["split_cdx_nfyy"] = split_cdx_nfyy

# 上海仁济CP40，化疗需要展示变异/野生型-2023.09.18
# wt 为野生型，只有一个为ref为杂合变异型，两个都不为ref为纯合变异型
def judge_chemo_genotype(info):
	dbsnp = info[0]
	alt = info[1]
	chemo_dict = {
		"rs3918290" : {"ref": "C", "wt": "C/C"},
		"rs55886062" : {"ref": "A", "wt": "A/A"},
		"rs67376798" : {"ref": "T", "wt": "T/T"},
		"rs75017182" : {"ref": "G", "wt": "G/G"},
		"rs10929302" : {"ref": "G", "wt": "G/G"},
		"rs4148323" : {"ref": "G", "wt": "G/G"},
		"rs8175347" : {"ref": "(TA)6", "wt": "(TA)6/(TA)6"}
		}
	split_alt = re.split("/", alt)
	if dbsnp in chemo_dict.keys():
		if split_alt[0] == chemo_dict.get(dbsnp)["ref"] and split_alt[1] == chemo_dict.get(dbsnp)["ref"]:
			return "野生型"
		else:
			if len(set(split_alt)) == 1:
				return "纯合变异型"
			else:
				return "杂合变异型"
		#elif split_alt[0] != chemo_dict.get(dbsnp)["ref"] and split_alt[1] != chemo_dict.get(dbsnp)["ref"]:
		#	return "纯合变异型"
		#else:
		#	return "杂合变异型"
	else:
		return "位点不在配置中"
jinja2.filters.FILTERS["judge_chemo_genotype"] = judge_chemo_genotype

# 孙逸仙116，子宫内膜癌时，结果小结和变异解读要把POLE、TP53变异升为I类，且变异解读第一段要加相关证据。
def syx_ec_var_class(info):
	var_list = info[0]
	class_num = info[1]
	level_I = [var for var in var_list if var["clinic_num_s"] == 5 or (var["clinic_num_s"] == 4 and var["gene_symbol"] in ["POLE", "TP53"])]
	level_II = [var for var in var_list if var["clinic_num_s"] == 4 and var["gene_symbol"] not in ["POLE", "TP53"]]
	if class_num == 1:
		return level_I
	elif class_num == 2:
		return level_II
	else:
		return []
jinja2.filters.FILTERS["syx_ec_var_class"] = syx_ec_var_class

# 中山人民HRDC检测结果小结-更新-2023.10.30
def summary_ZSRM_hrd_v3(info):
	# 阳性：基于本次送检样本，检出HRD状态【阳性/阴性】、BRCA gene_region hgvs_c hgvs_p I类变异、TP53 gene_region hgvs_c hgvs_p。
	# 阴性：基于本次送检样本，检出HRD【阳性/阴性】，未检出I类、II类变异。
	# 2023.10.30更新内容：兼容生信v1.3.0版本，hrd结果在gss字段中的情况。
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I类", 4 : "II类"}
	gss = info[2]

	var_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["hgvs_p"] != "p.?":
			var_result.append("{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			var_result.append("{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	
	#hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"
	if hrd and "var_id" in hrd.keys():
		hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"
	elif gss and "var_id" in gss.keys():
		hrd_result = "HRD状态阳性" if gss["var_id"] == "HRD+" else "HRD状态阴性"
	else:
		hrd_result = "未获取到HRD结果！"

	if var_result:
		return "基于本次送检样本，检出"+hrd_result+"、"+"、".join(var_result)+"。"
	else:
		return "基于本次送检样本，检出"+hrd_result+"，未检出I类、II类变异。"
jinja2.filters.FILTERS["summary_ZSRM_hrd_v3"] = summary_ZSRM_hrd_v3

# 获取变异列表中的Snvindel-2023.11.16
def filter_snvindel(var_list):
	return [var for var in var_list if var["bio_category"] == "Snvindel"]
jinja2.filters.FILTERS["filter_snvindel"] = filter_snvindel

# 获取变异列表中的RNA SV-2023.11.16
def filter_sv(var_list):
	return [var for var in var_list if var["bio_category"] in ["Sv", "PSeqRnaSv"]]
jinja2.filters.FILTERS["filter_sv"] = filter_sv

# 获取变异列表中的EGFR -2023.11.22
def filter_egfr(var_list):
	return [var for var in var_list if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR"]
jinja2.filters.FILTERS["filter_egfr"] = filter_egfr

#  XW5301-汇总基因检测结果，包含未检测到变异的基因-2023.11.27
def getSum_forXW5301(var_list, gene_rule, gene_transcript):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				var["transcript_primary_xw5031"] = gene_transcript.get(gene)
				result.append(var)
	for gene in set(gene_rule) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"transcript_primary_xw5031" : gene_transcript.get(gene)
		})
	return sorted(result, key = lambda i:gene_rule.index(i["gene_symbol"]))

# XW5301-snvindel排序-2023.11.27
def sort_for5301_snvindel(var_list):
	gene_rule = ["MET", "EGFR", "ALK", "KRAS", "ROS1", "RET", "ERBB2", "BRAF", "NRAS", "PIK3CA"]
	gene_transcript = {
		"MET" : "NM_000245",
		"EGFR" : "NM_005228",
		"ALK" : "NM_004304",
		"KRAS" : "NM_033360",
		"ROS1" : "NM_002944",
		"RET" : "NM_020975",
		"ERBB2" : "NM_004448",
		"BRAF" : "NM_004333",
		"NRAS" : "NM_002524",
		"PIK3CA" : "NM_006218"
	}
	return getSum_forXW5301(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for5301_snvindel"] = sort_for5301_snvindel

# XW5301-sv排序-2023.11.27
def sort_for5301_sv(var_list):
	gene_rule = ["ALK", "ROS1", "RET"]
	gene_transcript = {
		"ALK" : "NM_004304",
		"ROS1" : "NM_002944",
		"RET" : "NM_020975"
	}
	return getSum_forXW5301(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for5301_sv"] = sort_for5301_sv

# AD4701-tBRCA 检测结果小结-2023.11.30
def ad4701_tbrca_sum(var_list):
	clinic_5_count = len([var for var in var_list if var["clinic_num_g"] == 5])
	clinic_4_count = len([var for var in var_list if var["clinic_num_g"] == 4])
	if clinic_5_count and clinic_4_count:
		return "本次实验检出{0}个致病性变异和{1}个疑似致病性变异。".format(str(clinic_5_count), str(clinic_4_count))
	elif clinic_5_count and not clinic_4_count:
		return "本次实验检出{0}个致病性变异。".format(str(clinic_5_count))
	elif not clinic_5_count and clinic_4_count:
		return "本次实验检出{0}个疑似致病性变异。".format(str(clinic_4_count))
	else:
		return "本次实验未检出致病性或疑似致病性变异。"
jinja2.filters.FILTERS["ad4701_tbrca_sum"] = ad4701_tbrca_sum

# XW2402 CP 获取G12C结果 2023年12月27日
def filter_G12C(var_list):
	return [var for var in var_list if var['gene_symbol'] == 'KRAS' and var['hgvs_p'] == 'p.(G12C)']
jinja2.filters.FILTERS['filter_G12C'] = filter_G12C

# 获取变异列表中的CNV-2023年12月27日
def filter_cnv(var_list):
	return [var for var in var_list if var["bio_category"] == "Cnv"]
jinja2.filters.FILTERS["filter_cnv"] = filter_cnv

# XW5101 汇总结果
def getSum_forXW5101(var_list, gene_rule):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				result.append(var)
	for gene in set(gene_rule) - set(detect_gene):
		result.append({'gene_symbol' : gene})
	return sorted(result, key = lambda i:gene_rule.index(i["gene_symbol"]))

def sort_5101_snv(var_list):
    gene_rule = ['FLT3', 'KMT2A', 'NPM1', 'NUP98', 'TP53']
    return getSum_forXW5101(var_list, gene_rule)
jinja2.filters.FILTERS['sort_5101_snv'] = sort_5101_snv

def sort_5101_sv(var_list):
    gene_rule = ['KMT2A', 'NPM1', 'NUP98']
    return getSum_forXW5101(var_list, gene_rule)
jinja2.filters.FILTERS['sort_5101_sv'] = sort_5101_sv

def getSum_forXW4205(var_list, gene_rule, gene_transcript):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				var["transcript_primary_xw4205"] = gene_transcript.get(gene)
				result.append(var)
	for gene in sorted(set(gene_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"transcript_primary_xw4205" : gene_transcript.get(gene)
		})
	return result

def sort_for4205_snvindel(var_list):
	#gene_rule = ["MET", "EGFR", "ALK", "KRAS", "ROS1", "RET", "ERBB2", "BRAF", "NRAS", "PIK3CA"]
	gene_rule = ['ALK', 'BRAF', 'EGFR', 'ERBB2', 'KRAS', 'MET', 'NRAS', 'PIK3CA', 'RET', 'ROS1']
	gene_transcript = {
		"ALK" : "NM_004304",
		"BRAF" : "NM_004333",
		"EGFR" : "NM_005228",
		"ERBB2" : "NM_004448",
		"KRAS" : "NM_033360",
		"MET" : "NM_000245",
		"NRAS" : "NM_002524",
		"PIK3CA" : "NM_006218",
		"RET" : "NM_020975",
		"ROS1" : "NM_002944"
	}
	return getSum_forXW4205(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for4205_snvindel"] = sort_for4205_snvindel

def filter_snv_cnv(var_list):
	return [var for var in var_list if var['bio_category'] in ['Snvindel', 'Cnv']]
jinja2.filters.FILTERS['filter_snv_cnv'] = filter_snv_cnv

def filter_hd_xw1402(info):
	qc = info[0]
	hd = info[1]
	result = []
	detect_gene = []
	for var in hd:
		if var['var_auto_result'] == 'T':
			detect_gene.append(var['gene_symbol'])
	if not detect_gene:
		result.append('HD阴性')
	else:
		if 'MTAP' in detect_gene and 'CDKN2A' in detect_gene:
			result.append('HD阳性，检出MTAP纯合缺失、CDKN2A纯合缺失')
		elif 'MTAP' in detect_gene:
			result.append('HD阳性，检出MTAP纯合缺失，未检出CDKN2A纯合缺失')
		elif 'CDKN2A' in detect_gene:
			result.append('HD阳性，检出CDKN2A纯合缺失，未检出MTAP纯合缺失')

	if qc['dna_data_qc']['snp_cover_ratio_num'] < 0.9:
		result = ['N/A']
	else:
		if not hd:
			result = ['HD阴性']
	
	sum = 'N/A' if result == ['N/A'] else 'HD阴性' if result == ['HD阴性'] else ''.join(result)

	return sum
jinja2.filters.FILTERS['filter_hd_xw1402'] =  filter_hd_xw1402

def filter_xw1402_dna(var_list):
	dna_rule = ['ALK', 'BRAF', 'BRCA1', 'BRCA2', 'CLDN18', 'CDKN2A', 'EGFR', 'ERBB2', 'FGFR2','IDH1', 'KEAP1', 'KRAS', 'MTAP', 'MET', 'MLH1', 'MSH2', 'MSH6', 'NTRK1', 'NTRK2', 'NTRK3', 'PMS2', 'RET', 'ROS1', 'STK11']
	result = []
	detect_gene = []
	for var in var_list:
		if var['gene_symbol'] in dna_rule:
			detect_gene.append(var['gene_symbol'])
			result.append(var)
	for gene in sorted(set(dna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	
	return sorted(result, key = lambda i:dna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1402_dna'] = filter_xw1402_dna
	

def filter_xw1402_rna(var_list):
	rna_rule = ['ALK', 'BRAF', 'CLDN18', 'EGFR', 'FGFR2', 'MET', 'NRG1', 'NTRK1', 'NTRK2', 'NTRK3', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(',', var['gene_symbol'])):
			detect_gene.append(gene)
			if gene in rna_rule:
				var['gene_symbol'] = gene
				if gene == 'MET':
					var['gene_symbol'] = 'MET'
				result.append(var)
	for gene in sorted(set(rna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	return sorted(result, key = lambda i:rna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1402_rna'] = filter_xw1402_rna

def filter_xw0258_hd(var_list):
	result = []
	if var_list:
		for var in var_list:
			if var['var_auto_result'] == "T" and var['type'] == "homodel":
				result.append(var)
	return result
jinja2.filters.FILTERS['xw0258_hd'] = filter_xw0258_hd

def filter_xw1404_dna(var_list):
	dna_rule = ['ALK', 'BRAF', 'CDKN2A', 'EGFR', 'ERBB2', 'KRAS', 'MTAP', 'MET', 'NTRK1', 'NTRK2', 'NTRK3', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		if var['gene_symbol'] in dna_rule:
			detect_gene.append(var['gene_symbol'])
			result.append(var)
	for gene in sorted(set(dna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	
	return sorted(result, key = lambda i:dna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1404_dna'] = filter_xw1404_dna
	

def filter_xw1404_rna(var_list):
	rna_rule = ['ALK', 'BRAF', 'EGFR', 'MET', 'NTRK1', 'NTRK2', 'NTRK3', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(',', var['gene_symbol'])):
			detect_gene.append(gene)
			if gene in rna_rule:
				if gene == 'MET':
					var['gene_symbol'] = 'MET'
				result.append(var)
	for gene in sorted(set(rna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	return sorted(result, key = lambda i:rna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1404_rna'] = filter_xw1404_rna

def filter_xw1405_dna(var_list):
	dna_rule = ['ALK', 'BRAF', 'BRCA1', 'BRCA2', 'CLDN18', 'CDKN2A', 'EGFR', 'ERBB2', 'FGFR2', 'HRAS', 'IDH1', 'IDH2', 'KEAP1', 'KRAS', 'MTAP', 'MET', 'MLH1', 'MSH2', 'MSH6', 'NRAS', 'NTRK1', 'NTRK2', 'NTRK3', 'PALB2', 'PMS2', 'RET', 'ROS1', 'STK11']
	result = []
	detect_gene = []
	for var in var_list:
		if var['gene_symbol'] in dna_rule:
			detect_gene.append(var['gene_symbol'])
			result.append(var)
	for gene in sorted(set(dna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	
	return sorted(result, key = lambda i:dna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1405_dna'] = filter_xw1405_dna

def filter_xw1405_rna(var_list):
	rna_rule = ['ALK', 'BRAF', 'CLDN18', 'EGFR', 'FGFR2', 'MET', 'NTRK1', 'NTRK2', 'NTRK3', 'NRG1', 'RET', 'ROS1']
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(',', var['gene_symbol'])):
			detect_gene.append(gene)
			if gene in rna_rule:
				if gene == 'MET':
					var['gene_symbol'] = 'MET'
				result.append(var)
	for gene in sorted(set(rna_rule) - set(detect_gene)):
		result.append({
			"gene_symbol" : gene,
			"bio_category" : ""
		})
	return sorted(result, key = lambda i:rna_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS['xw1405_rna'] = filter_xw1405_rna

def xw1402_hd(hd):
	target = ['MTAP', 'CDKN2A']
	if hd:
		detect_gene = [var['gene_symbol'] for var in hd] if hd else []
		notdetect = list(set(target) - set(detect_gene))
		detect_gene = list(set(detect_gene))
		if len(detect_gene) == 2:
			return hd
		elif len(detect_gene) == 1:
			hd.append({
				'gene_symbol' : ''.join(notdetect),
				'region' : '/',
				'var_auto_result' : 'F'
			})
	hd = sorted(hd, key = lambda i:target.index(i["gene_symbol"]))
	return hd
jinja2.filters.FILTERS['xw1402_hd'] = xw1402_hd

def ad3101_inter(info):
	var = info[0]
	sample = info[1]
	if sample['product_name'] == 'AD3101' and sample['prod_names'] == 'Master Panel（组织）':
		return var["var_somatic"]["level_I"] + var["var_germline"]["level_5"] + var["var_somatic"]["level_II"] + var["var_somatic"]["level_onco_nodrug"] + var["var_germline"]["level_4"]
jinja2.filters.FILTERS['ad3101_inter'] = ad3101_inter

def approval_regimen_biomarker_v2(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], i["hgvs_p"]))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			# 2024.08.30-CNV 区分Loss，其他的写扩增
			if i["cnv_type"] == "Loss" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "Gain" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段重复".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HeteDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 杂合大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HomoDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 纯合大片段缺失".format(i["gene_symbol"]))
			else:
				result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
					# 加一个兼容-2024.01.25
					# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
					if not i["five_prime_region"]:
						i["five_prime_region"] = re.split("_", re.split("-", i["hgvs"])[0])[-1]
						i["three_prime_region"] = re.split("_", i["hgvs"])[-1]
					# 兼容完成-2024.01.25
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF 野生型")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD阴性")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD阳性")
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	# 2024.02.27更新-MET DNA/RNA共检时都要展示
	#if judge_mergeMET:
	#	if "MET exon14 跳跃" in result_redup:
	#		result_redup.remove("MET exon14 跳跃")
	# 2024.02.27更新完成
	if not result_redup:
		result_redup = ["-"]
	rt = RichText()
	rt.add("\n".join(result_redup))
	return rt
jinja2.filters.FILTERS["approval_regimen_biomarker_v2"] = approval_regimen_biomarker_v2

def get_drug_mechanism_cn(drug_details):
	result = []
	for i in drug_details:
		if i["drug_mechanism_cn"] and i["drug_mechanism_cn"] not in result and i["drug_name"] != "化疗;Chemotherapy":
			result.append(i["drug_mechanism_cn"])
	return result
jinja2.filters.FILTERS["get_drug_mechanism_cn"] = get_drug_mechanism_cn

def mrd_summary(var_list):
	result = []
	var_list = var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"]
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"] + " " + var["freq_str"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"] + " " + var["freq_str"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+" 扩增" + " " + var["cn_mean"])
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃"+" "+var["freq_str"])
			else:
				if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in result:
					result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合"+" "+var["freq_str"])
	if result:
		return ", ".join(result)
	else:
		return "本次未检测到相关驱动突变"
jinja2.filters.FILTERS["mrd_summary"] = mrd_summary
	