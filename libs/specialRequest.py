#-*- coding:gbk -*-
import re
from libs.AAstrans import strans
from libs.rule import s_var_rule
#from libs.getConfig import get_FDZS_database
#from libs.getConfig import get_FJZL_database
from libs.getConfig import getconfigxlsx, get_shrj_geneinfo
from libs.rule import S_level
import copy
from collections import defaultdict
'''
	处理模板中各种特殊要求 
'''

# 该函数可更新
def varInter_FJZL(hgvs_p, var_type, gene_region, hgvs_c):
	'''
	福建省肿瘤医院（只有单独组织和单独全血）
	1. 变异描述：该患者外周血基因组DNA样品检测出BRCA2基因 < 14号外显子发生c.7409_7410insT，该突变导致基因编码蛋白第2471位氨基酸由苏氨酸突变为组氨酸并于2474位发生提前终止，
				可能形成功能损伤或失活的蛋白。>
	< >内的描述由该函数实现
	'''
	#转换gene_region，如exon转换为12号外显子
	if re.search("exon", gene_region):
		gene_region_cn = gene_region.replace("exon", "")+"号外显子"
	elif re.search("intron", gene_region):
		gene_region_cn = gene_region.replace("intron", "")+"号内含子"
	else:
		gene_region_cn = gene_region

	inter1 = gene_region_cn+"发生"+hgvs_c
	inter2 = strans(hgvs_p, var_type)
	inter = inter1+"，"+inter2 if inter2 else inter1
	return inter

# 该函数可在报告模板中实现
def var_summary_FJZL(var_brca):
	'''
	福建省肿瘤医院（只有单独组织和单独全血）
	1. 该患者外周血基因组DNA样品检测出< BRCA1基因致病性变异 >
	'''	
	# 统计致病性、非致病性、BRCA1、BRCA
	gene_rule = ["BRCA1", "BRCA2"]
	gene_list = [i["gene_symbol"] for i in var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"]+var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"]]
	level_rule = ["疑似致病性变异", "致病性变异", "意义不明确"]
	level_dict = {4 : "疑似致病性变异", 5 : "致病性变异", 3 : "意义不明确"}
	level_list = [level_dict.get(i["clinic_num_g"]) for i in var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"]+var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]]
	if var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
		level_list.append(level_dict.get(4))
	gene_str = "和".join(sorted(list(set(gene_list)), key=gene_rule.index))
	level_str = "和".join(sorted(list(set(level_list)), key=level_rule.index))

	return gene_str+"基因"+level_str

# 该函数可在报告模板中实现
def var_summary_GDRM(var_brca):
	'''
	广东人民（只有单独组织和单独全血）
	全血：
		阳性：基于本次送检样本，本次检验检测到X个致病性变异和X个疑似致病性变异
		阴性：基于本次送检样本，未检测到致病性或疑似致病性变异
	组织：
		阳性：基于本次送检样本，本次检验检测到X个具有临床意义的变异
		阴性：基于本次送检样本，未检测到具有临床意义的变异
	'''
	level45_num_list_G = []
	level5_num = len(var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"])
	level4_num = len(var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"])
	if level5_num > 0:
		level45_num_list_G.append(str(level5_num)+"个致病性变异")
	if level4_num > 0:
		level45_num_list_G.append(str(level4_num)+"个疑似致病性变异")
	
	summary_G = "基于本次送检样本，本次检验检测到"+"和".join(level45_num_list_G) if level45_num_list_G else "基于本次送检样本，未检测到致病性或疑似致病性变异"
	summary_S = "基于本次送检样本，本次检验检测到"+str(level5_num+level4_num)+"个具有临床意义的变异" if level5_num+level4_num > 0 else "基于本次送检样本，未检测到具有临床意义的变异"
	
	return summary_G, summary_S
	
# 可在报告模板中实现
def var_GZZL(var_brca):
	'''
	贵州肿瘤（只有单独组织和单独全血）
	全血：
		本次实验检测到一些（意义不明确变异、疑似良性变异、良性变异）
	组织：
		本次实验检测到一些（临床意义不明、疑似良性变异、良性变异）
	'''
	level1_3_G = []
	level1_3_S = []
	if var_brca["snv_s"]["B1_L3"] or var_brca["snv_s"]["B2_L3"]:
		level1_3_G.append("意义不明确变异")
		level1_3_S.append("临床意义不明")
	if var_brca["snv_s"]["B1_L2"] or var_brca["snv_s"]["B2_L2"]:
		level1_3_G.append("疑似良性变异")
		level1_3_S.append("疑似良性变异")
	if var_brca["snv_s"]["B1_L1"] or var_brca["snv_s"]["B2_L1"]:
		level1_3_G.append("良性变异")
		level1_3_S.append("良性变异")
	level13_str_G = "、".join(level1_3_G)	
	level13_str_S = "、".join(level1_3_S)

	return level13_str_G, level13_str_S

def var_summary_SHRJ(var_brca):
	'''
	上海仁济（单独组织、单独全血）
	检测结果需要拼接
	'''
	varInter_list = []
	for var in var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"]+var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]:
		if var["hgvs_p"] != "p.?":
			inter = "检测到"+var["gene_symbol"]+"基因的"+var["hgvs_c"]+"（"+var["hgvs_p"]+"）突变，"+var["variant_desc_cn"]
		else:
			inter = "检测到"+var["gene_symbol"]+"基因的"+var["hgvs_c"]+"突变，"+var["variant_desc_cn"]
		if re.search("Splicing|Intronic", var["type"]):
			inter = inter[0:-1]+"，可能造成异常剪接。"
		if re.search("FrameShift|Nonsense", var["type"]):
			inter = inter[0:-1]+"，可能形成功能损伤或失活的蛋白。"
		varInter_list.append(inter)

	for var in var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"]:
		varInter_list.append("检测到"+var["gene_symbol"]+"基因的大片段缺失，该类变异为疑似致病性变异，可能影响"+var["gene_symbol"]+"的蛋白功能。")
	
	return "".join(varInter_list)

def varInfo_XAJDY(gene_region, type_cn):
	'''
	西安交大一
	变异需增加：第XX号外显子/内含子XX突变
	'''
	region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}

	region_list_en = re.split("_", gene_region)
	region_list_cn = []
	for i in region_list_en:
		if re.search("exon", i):
			region_list_cn.append("第"+i.replace("exon", "")+"号外显子")
		elif re.search("intron", i):
			region_list_cn.append("第"+i.replace("intron", "")+"号内含子")
		else:
			region_cn = region_dict[i] if i in region_dict.keys() else i
			region_list_cn.append(region_cn)
	# 先兼容为空的情况
	type_cn = type_cn if type_cn else ""
	return "到".join(region_list_cn)+type_cn

# 看起来可在模板中实现
def var_SD(var_brca):
	'''
	山东省立/山东齐鲁（单独组织、单独全血、配对）
	全血：
		本次实验检测到一些（疑似良性变异、良性变异）
	组织：
		本次实验检测到一些（良性/疑似良性变异）
	配对：
		本次实验检测到一些（疑似良性变异、良性变异、良性/疑似良性变异）
	'''

	level1_2_G = []
	level1_2_S = []
	g_stran = {2 : "疑似良性变异", 1 : "良性变异"}
	g_rule = ["疑似良性变异", "良性变异"]
	s_stran = {2 : "良性/疑似良性变异", 1 : "良性/疑似良性变异"}
	s_rule = ["疑似良性变异", "良性/疑似良性变异", "良性变异"]

	# 胚系 适用于单独全血或配对中的全血报告
	for var in var_brca["snv_m"]["B1_G_L2"]+var_brca["snv_m"]["B2_G_L2"]+var_brca["snv_m"]["B1_G_L1"]+var_brca["snv_m"]["B2_G_L1"]:
		level1_2_G.append(g_stran.get(var["clinic_num_g"]))
	level1_2_G_str = "、".join(sorted(list(set(level1_2_G)), key=g_rule.index))

	# 体细胞 适用于单独组织或配对中的组织报告
	for var in var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"]:
		if var["var_origin"] == "germline":
			level1_2_S.append(g_stran.get(var["clinic_num_g"]))
		else:
			level1_2_S.append(s_stran.get(var["clinic_num_s"]))
	level1_2_S_str = "、".join(sorted(list(set(level1_2_S)), key=s_rule.index))

	return level1_2_G_str, level1_2_S_str

def var_summary_WHXH(var_brca):
	'''
	武汉协和（单独组织、单独全血、配对），均用胚系的等级
	1. 总结：
		阴性：基于本次送检样本，未检出有害突变或疑似有害突变。
		阳性：本次实验检测到X个XX突变。
	2. level_3总结：基于数据库信息均为 < 意义不明确突变、疑似无害突变、无害突变、多态性改变 >
	'''
	# 4/5类总结
	level4_5 = []
	level_4_count = len(var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B2_L4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"])
	level_5_count = len(var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L5"])
	if level_5_count:
		level4_5.append(str(level_5_count)+"个有害突变")
	if level_4_count:
		level4_5.append(str(level_4_count)+"个疑似有害突变")
	level4_5_summary = "和".join(level4_5)

	#1/2/3类总结-加一个仅胚系变异的1/2/3总结，适用山东模板
	level_stran = {3 : "意义不明确突变", 2 : "疑似无害突变", 1 : "无害突变"}
	sort_rule = ["意义不明确突变", "疑似无害突变", "无害突变", "多态性改变"]
	def level_sum(var1_3_list):
		level1_3 = []
		for var in var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"]:
			if var["clinic_num_g"] in [2,3]:
				level1_3.append(level_stran.get(var["clinic_num_g"]))
			else:
				if "tag" in var.keys() and re.search("Polymorphism", var["tag"]):
					level1_3.append("多态性改变")
				else:
					level1_3.append(level_stran.get(var["clinic_num_g"]))
		return "、".join(sorted(list(set(level1_3)), key=sort_rule.index))
	all_origin_var = var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"]
	germline_var = var_brca["snv_m"]["B1_G_L3"]+var_brca["snv_m"]["B2_G_L3"]+var_brca["snv_m"]["B1_G_L2"]+var_brca["snv_m"]["B2_G_L2"]+var_brca["snv_m"]["B1_G_L1"]+var_brca["snv_m"]["B2_G_L1"]
	level1_3_summary = level_sum(all_origin_var)
	level1_3_summary_G = level_sum(germline_var)

	return level4_5_summary, level1_3_summary, level1_3_summary_G

# 看起来可在模板中实现	
def var_level123_SDZL(var_brca):
	'''
	山东肿瘤（单独组织、单独全血、配对）
	level1_3总结
	单独全血：均为 < 意义不明确变异、疑似良性变异、良性变异 >
	单独组织：意义不明确变异、良性/可能良性变异
	配对：意义不明确变异、疑似良性变异、良性变异、良性/可能良性变异

	'''
	g_stran = {3 : "意义不明确变异", 2 : "疑似良性变异", 1 : "良性变异"}
	s_stran = {3 : "意义不明确变异", 2 : "良性/可能良性变异", 1 : "良性/可能良性变异"}
	sort_rule = ["意义不明确变异", "疑似良性变异", "良性/可能良性变异", "良性变异"]
	# 可用于BRCA全血
	level1_3_G = [g_stran.get(var["clinic_num_g"]) for var in var_brca["snv_m"]["B1_G_L3"]+var_brca["snv_m"]["B2_G_L3"]+var_brca["snv_m"]["B1_G_L2"]+var_brca["snv_m"]["B2_G_L2"]+var_brca["snv_m"]["B1_G_L1"]+var_brca["snv_m"]["B2_G_L1"]]
	level1_3_G_str = "、".join(sorted(list(set(level1_3_G)), key=sort_rule.index))
	# 可用于BRCA组织或配对
	level1_3_S = [s_stran.get(var["clinic_num_s"]) for var in var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"] if var["var_origin"] != "germline"]
	level1_3_S += [g_stran.get(var["clinic_num_g"]) for var in var_brca["snv_s"]["B1_L3"]+var_brca["snv_s"]["B2_L3"]+var_brca["snv_s"]["B1_L2"]+var_brca["snv_s"]["B2_L2"]+var_brca["snv_s"]["B1_L1"]+var_brca["snv_s"]["B2_L1"] if var["var_origin"] == "germline"]
	level1_3_S_str = "、".join(sorted(list(set(level1_3_S)), key=sort_rule.index))

	return level1_3_G_str, level1_3_S_str

def brca_judge_inApprovalTumor(var_brca):
	'''
	判断BRCA项目样本是否在获批癌种里（主要用于体细胞变异，判断最高等级是I类还是II类）
	判断规则：
	若治疗方案中最高等级为A、B，判定为在获批癌种里，C、D判定为不在里面
	'''
	var = var_brca["snv_s"]["B1_L4"]+var_brca["snv_s"]["B1_L5"]+var_brca["snv_s"]["B2_L4"]+var_brca["snv_s"]["B2_L5"]
	level_list = []
	for i in var:
		# 这边若治疗方案为空，会被判定成潜在临床意义！！
		if "evi_sum" in i.keys() and i["evi_sum"] and "evi_split" in i["evi_sum"].keys() and i["evi_sum"]["evi_split"] and "Predictive" in i["evi_sum"]["evi_split"].keys():
			level_list += [j["evi_conclusion_simple"] for j in i["evi_sum"]["evi_split"]["Predictive"]]
	if level_list:
		if set(["A", "B"]) & set(level_list):
			return "intumor"
		else:
			return "outtumor"

def judge_GA_tumor_KNB(var_data, tumor_list):
	'''
	复旦中山CP40：癌种为消化道相关肿瘤时，需额外判断KNB
	'''
	# 更新KNB判定规则：KRAS/NRAS不为致癌/疑似致癌，且不在exon2/3/4, 未检测到BRAF V600E-2022.11.21
	GA_tumor_list = ["食管癌","胃癌","胃食管交界癌","贲门癌","小肠癌","胆道癌","胆管癌","胆囊癌","肝癌","胃肠道间质瘤","胰腺癌","阑尾腺癌",\
					 "壶腹部腺癌","十二指肠腺癌","食管鳞状细胞癌","食管腺癌","胃腺癌","贲门癌"]
	level1_2_var = [var for var in var_data if var["clinic_num_s"] in [5, 4] and var["bio_category"] == "Snvindel"]
	for var in level1_2_var:
		if any([var["gene_symbol"] in ["KRAS", "NRAS"] and var["gene_region"] in ["exon2", "exon3", "exon4"],
				var["gene_symbol"] == "BRAF" and var["hgvs_p"] == "p.(V600E)"]):
			break
	else:
		if set(tumor_list) & set(GA_tumor_list):
			return "yes"
		else:
			return 0

def sv_shfk(var_data):
	'''
	上海肺科CP40：检测到融合时，其他说明中要写"本次检测到gene1-gene2融合。具体的融合型为gene1:exon-gene2:exon"
	2022.07.27 新增，MET 14 skipping也要展示。DNA和RNA均检测到，则展示DNA结果解读+"本次实验在RNA水平也检测到MET exon14 skipping。"；只有RNA检测到时则放RNA解读，仅DNA检测到暂时无法判定，不放解读。
	'''
	sv_var = [var for var in var_data if var["bio_category"] == "Sv" and not (var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET")]
	RNA_met = [var for var in var_data if var["bio_category"] == "Sv" and var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET"]
	DNA_RNA_met = [var for var in var_data if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "MET" and "judge_mergeMET" in var.keys() and var["judge_mergeMET"]]
	sv_str_list = []
	for var in sv_var:
		if "var_desc_merge" in var.keys():
			sv_str_list.append("检测到"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合，具体的融合型为"+var["var_desc_merge"])
	rna_str_list = []
	for var in RNA_met:
		rna_str_list.append(var["variant_desc_cn"]+var["variant_interpret_cn"])
	dna_rna_str_list = []
	for var in DNA_RNA_met:
		if var["hgvs_p"] == "p.?":
			dna_rna_str_list.append("本次实验检测到"+var["gene_symbol"]+" "+var["hgvs_c"]+"变异，"+var["variant_desc_cn"]+var["variant_interpret_cn"]+"本次实验在RNA水平也检测到MET exon14 skipping。")
		else:
			dna_rna_str_list.append("本次实验检测到"+var["gene_symbol"]+" "+var["hgvs_p"]+"变异，"+var["variant_desc_cn"]+var["variant_interpret_cn"]+"本次实验在RNA水平也检测到MET exon14 skipping。")
	
	sv_str = "本次"+"；".join(sv_str_list)+"。" if sv_str_list else ""
	rna_str = ";".join(rna_str_list)
	dna_rna_str = ";".join(dna_rna_str_list)
	return sv_str + rna_str + dna_rna_str

def CP40_split_gene(var_data):
	'''
	聊城、厦门市一CP40，报告中需要展示成10gene+30gene（只包含I/II类变异）
	special["var_cp40_split"] = {
		gene10_level_I_II : [],
		gene30_level_I_II : []
	}
	'''
	gene10_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	var_levelI_II = var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]
	var_levelIII = var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"]
	# 新增III类-适用复旦中山厦门医院
	gene10_level_I_II = [var for var in var_levelI_II if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_I_II = [var for var in var_levelI_II if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene10_level_III = [var for var in var_levelIII if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_III = [var for var in var_levelIII if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	# III类拆分为肿瘤发生发展相关和III类变异-适用内江二院CP40-2023.06.09
	gene10_level_nj_onco_nodrug = [var for var in var_data["var_somatic"]["level_onco_nodrug"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_nj_onco_nodrug = [var for var in var_data["var_somatic"]["level_onco_nodrug"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene10_level_nj_III = [var for var in var_data["var_somatic"]["level_III"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_nj_III = [var for var in var_data["var_somatic"]["level_III"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	return gene10_level_I_II, gene30_level_I_II, gene10_level_III, gene30_level_III, gene10_level_nj_onco_nodrug, gene30_level_nj_onco_nodrug, gene10_level_nj_III, gene30_level_nj_III

def CP40_FJFY_summary(level_I, level_II, level_III, level_onco_nodrug):
	'''
	福建医科附一CP40，报告中需要展示I/II类变异详情
	I类：XX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为X，对XX、XX可能敏感/可能耐药（A级）
	II类：XX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为X

	snvindel：XX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为X
	cnv：XX基因检测到扩增，拷贝数为X
	sv：XX基因检测到XX-XX融合，拷贝数为X copies
	'''
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}
		if var["bio_category"] == "Cnv":
			var_info = var["gene_symbol"]+"基因检测到扩增，拷贝数为"+str(var["cn_mean"])
		elif var["bio_category"] == "Sv":
			var_info = var["gene_symbol"]+"基因检测到"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合，拷贝数为"+str(var["copies"])+" copies"
		elif var["bio_category"] == "Snvindel":
			# 加一个MET exon 14跳跃突变-2023.07.13
			if var["gene_symbol"] == "MET" and "exon14 skipping" in var["variant_interpret_cn"]:
				if var["hgvs_p_FJFY"] != "p.?":
					var_info  = "MET基因14号外显子跳跃突变"+var["hgvs_c"]+"("+var["hgvs_p_FJFY"]+")，突变丰度为"+str(var["freq_str"])
				else:
					var_info = "MET基因14号外显子跳跃突变"+var["hgvs_c"]+"，突变丰度为"+str(var["freq_str"])
			# MET exon14跳跃添加结束-2023.07.13
			else:
				region_list_en = re.split("_", var["gene_region"])
				region_list_cn = []
				for i in region_list_en:
					if re.search("exon", i):
						region_list_cn.append(i.replace("exon", "")+"号外显子")
					elif re.search("intron", i):
						region_list_cn.append(i.replace("intron", "")+"号内含子")
					else:
						region_cn = region_dict[i] if i in region_dict.keys() else i
						region_list_cn.append(region_cn)
				if var["hgvs_p_FJFY"] != "p.?":
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_p_FJFY"]+"，突变丰度为"+str(var["freq_str"])
				else:
					if var["type_cn"] != "--":
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到"+ var["type_cn"] + var["hgvs_c"]+"突变，突变丰度为"+str(var["freq_str"])
					else:
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到"+ var["hgvs_c"]+"突变，突变丰度为"+str(var["freq_str"])
		return var_info

	def regimen_info(var):
		regimen_info = ""
		sense_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if \
							re.search("Sensitive", i["clinical_significance"]) and \
							i["evi_conclusion_simple"]=="A"] if \
								"Predictive" in var["evi_sum"]["evi_split"].keys() else \
								[]
		resis_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if \
							re.search("Resistant", i["clinical_significance"]) and \
							i["evi_conclusion_simple"]=="A"] if \
								"Predictive" in var["evi_sum"]["evi_split"].keys() else \
								[]

		sense_regimen_str = "对" + "、".join(sense_regimen) + "可能敏感（A级）" if sense_regimen else ""
		resis_regimen_str = "对" + "、".join(resis_regimen) + "可能耐药（A级）" if resis_regimen else ""

		return sense_regimen_str, resis_regimen_str
	# I类变异
	var_levelI_sum = []

	# 新增开始-2022.11.08
	lc10_gene = ["BRAF", "EGFR", "ERBB2", "KRAS", "MET", "ALK", "ROS1", "RET", "NRAS", "PIK3CA"]
	var_levelI_sum_lc10_list = []
	var_levelI_sum_other_list = []
	# 新增结束-2022.11.08

	for var in level_I:
		tmp_list = []
		var_info = var_info_stran(var)
		sense_regimen_str, resis_regimen_str = regimen_info(var)
		tmp_list.append(var_info)
		if sense_regimen_str:
			tmp_list.append(sense_regimen_str)
		if resis_regimen_str:
			tmp_list.append(resis_regimen_str)
		var_levelI_sum.append("，".join(tmp_list))

		# 新增开始-2022.11.08
		for gene in re.split(",",var["gene_symbol"]):
			if gene in lc10_gene:
				if "，".join(tmp_list) and "，".join(tmp_list) not in var_levelI_sum_lc10_list:
					var_levelI_sum_lc10_list.append("，".join(tmp_list))
			else:
				if "，".join(tmp_list) and "，".join(tmp_list) not in var_levelI_sum_other_list:
					var_levelI_sum_other_list.append("，".join(tmp_list))
		# 新增结束-2022.11.08

	level_I_sum = "；".join(var_levelI_sum) if var_levelI_sum else ""
	var_levelI_sum_lc10 = "；".join(var_levelI_sum_lc10_list) if var_levelI_sum_lc10_list else ""
	var_levelI_sum_other = "；".join(var_levelI_sum_other_list) if var_levelI_sum_other_list else ""

	# II类变异
	level_II_sum = "；".join([var_info_stran(var) for var in level_II]) if level_II else ""
	# 新增开始-2022.11.08
	level_II_sum_lc10_list = []
	level_II_sum_other_list = []
	for var in level_II:
		for gene in re.split(",", var["gene_symbol"]):
			if gene in lc10_gene:
				if var_info_stran(var) and var_info_stran(var) not in level_II_sum_lc10_list:
					level_II_sum_lc10_list.append(var_info_stran(var))
			else:
				if var_info_stran(var) and var_info_stran(var) not in level_II_sum_other_list:
					level_II_sum_other_list.append(var_info_stran(var))

	level_II_sum_lc10 = "；".join(level_II_sum_lc10_list) if level_II_sum_lc10_list else ""
	level_II_sum_other = "；".join(level_II_sum_other_list) if level_II_sum_other_list else level_II_sum_other_list
	# 其他基因结果还需要包含三类和肿瘤发生发展相关
	level_III_sum_other = []
	for var in level_III+level_onco_nodrug:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in lc10_gene:
				level_III_sum_other.append(var)
	# 新增结束-2022.11.08

	result = {
		"level_I_sum" : level_I_sum,
		"level_II_sum" : level_II_sum,
		"level_I_lc10" : var_levelI_sum_lc10,
		"level_II_lc10" : level_II_sum_lc10,
		"level_I_other" : var_levelI_sum_other,
		"level_II_other" : level_II_sum_other,
		"level_III_other" : level_III_sum_other
	}
	
	return result
		
def Master_sum(var_data):
	# 统计体细胞变异个数（DNA/RNA分开统计，取消！！）
	result = s_var_rule(var_data)

	#level_I = [var for var in var_data if var["clinic_num_s"] == 5 and var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"]]
	#level_II = [var for var in var_data if var["clinic_num_s"] == 4 and var["var_origin"] != "germline" and "evi_sum" in var.keys() and var["evi_sum"]["evi_split"]]
	#level_onco_nodrug = [var for var in var_data if var["clinic_num_s"] in [5, 4] and var["var_origin"] != "germline" and "evi_sum" in var.keys() and not var["evi_sum"]["evi_split"]]
	#level_III = [var for var in var_data if var["clinic_num_s"] == 3 and var["var_origin"] != "germline"]

	def var_count(var_list):
		num = 0
		for var in var_list:
			num += 1
			# DNA/RNA融合外显子一样，计数为1-2023.02.08
	#		if "rna_detect" in var.keys():
	#			num += 1 
		return num
	return var_count(result["level_I"]), var_count(result["level_II"]), var_count(result["level_onco_nodrug"]), var_count(result["level_III"])

def varInfo_SYX(gene_region, type_cn, type, hgvs_c, hgvs_p):
	'''
	孙逸仙
	变异需增加：XX号外显子/内含子XX hgvs_p/hgvs_c XX突变
	'''
	region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}

	type_stran = {
		"Intronic" : "内含子突变",
		"3'UTR" : "3'UTR突变",
		"5'UTR" : "5'UTR突变",
		"3'FLANKING" : "非编码区突变",
		"5'FLANKING" : "非编码区突变"
	}
	region_list_en = re.split("_", gene_region)
	region_list_cn = []
	for i in region_list_en:
		if re.search("exon", i):
			region_list_cn.append(i.replace("exon", "")+"号外显子")
		elif re.search("intron", i):
			region_list_cn.append(i.replace("intron", "")+"号内含子")
		else:
			region_cn = region_dict[i] if i in region_dict.keys() else i
			region_list_cn.append(region_cn)
	# 先兼容为空的情况
	type_cn = type_cn if type_cn and type_cn != "--" else type_stran.get(type, "突变") 
	var_info = hgvs_p if hgvs_p != "p.?" else hgvs_c
	return "到".join(region_list_cn) + var_info + type_cn

def varInfo_FDZS(gene_symbol, gene_region, variant_desc_cn, config):
	'''
	复旦中山
	1）基因描述使用配置表
	2）遗传疾病使用配置表
	3）遗传疾病描述信息使用配置表
	4）风险使用配置表
	5）gene_region转化为中文，如"exon1"转化为1号外显子
	6）变异描述，截取部分信息，仅展示"导致基因编码蛋白第1056位氨基酸由谷氨酰胺突变为终止密码子"
	'''
	# 1. 获取配置表信息
	#FDZS_gene_dict = get_FDZS_database(config)
	FDZS_gene_dict, FJZL_dict, io_dict = getconfigxlsx(config)
	result = FDZS_gene_dict[gene_symbol] if gene_symbol in FDZS_gene_dict.keys() and FDZS_gene_dict[gene_symbol] else {}
	# 2. 获取gene_region
	result["gene_region_cn"] = ""
	if gene_region:
		region_dict = {
			"exon" : "外显子",
			"intron" : "内含子",
			"3'UTR" : "3'UTR",
			"5'UTR" : "5'UTR",
			"3'FLANKING" : "非编码区",
			"5'FLANKING" : "非编码区"
			}

		region_list_en = re.split("_", gene_region)
		region_list_cn = []
		for i in region_list_en:
			if re.search("exon", i):
				region_list_cn.append(i.replace("exon", "")+"号外显子")
			elif re.search("intron", i):
				region_list_cn.append(i.replace("intron", "")+"号内含子")
			else:
				region_cn = region_dict[i] if i in region_dict.keys() else i
				region_list_cn.append(region_cn)
		result["gene_region_cn"] = "到".join(region_list_cn)
	# 3. 处理变异描述
	var_desc = re.split("，",variant_desc_cn) if variant_desc_cn else []
	result["var_desc"] = var_desc[1].replace("。","") if var_desc and len(var_desc) >= 2 else ""

	return result

def var_summary_CQFY(var_brca):
	'''
	重庆附一进院BRCA（仅单独全血）
	***不考虑MLPA***
	1. 若两个基因均有检出且不同变异等级，分开写
	       该样本检测到BRCA1基因为2类--可能良性（likely benign）的变异；BRCA2基因为1类--良性（benign）的变异。
	2. 若两个基因均有检出且为同一变异等级，写一起
	       该样本检测到BRCA1、BRCA2基因为2类--可能良性（likely benign）的变异。
	3. 若同一个基因存在多个3类以上变异位点的情况，均列出，且4,5类标红展示
	       该样本检测到BRCA1基因有5类--致病（pathogenic）和3类--意义不明（uncertain significance）的变异；BRCA2基因为1类--良性（benign）的变异。
	'''	
	# 等级转化
	level_dict = {
		"5" : "5类--致病（pathogenic）",
		"4" : "4类--可能致病（likely pathogenic）",
		"3" : "3类--意义未明(uncertain significance)",
		"2" : "2类--可能良性（likely benign）",
		"1" : "1类--良性（benign）"
		}

	result = {}

	# 获取BRCA1变异的所有等级（不包含MLPA）
	brca1_level_list = [str(int(i["clinic_num_g"])) for i in var_brca["snv_s"]["B1_L5"] + \
															 var_brca["snv_s"]["B1_L4"] + \
															 var_brca["snv_s"]["B1_L3"] + \
															 var_brca["snv_s"]["B1_L2"] + \
															 var_brca["snv_s"]["B1_L1"]]
	# 获取BRCA2变异的所有等级（不包含MLPA）
	brca2_level_list = [str(int(i["clinic_num_g"])) for i in var_brca["snv_s"]["B2_L5"] + \
															 var_brca["snv_s"]["B2_L4"] + \
															 var_brca["snv_s"]["B2_L3"] + \
															 var_brca["snv_s"]["B2_L2"] + \
															 var_brca["snv_s"]["B2_L1"]]
	# L21仅提取12中的最高等级，字符串
	result["B1_L21"] = level_dict["2"] if "2" in brca1_level_list else level_dict["1"]
	result["B2_L21"] = level_dict["2"] if "2" in brca2_level_list else level_dict["1"]
	# L3 仅提取3，字符串
	result["B1_L3"] = level_dict["3"] if "3" in brca1_level_list else ""
	result["B2_L3"] = level_dict["3"] if "3" in brca2_level_list else ""
	# L54 提取检出的45，字符串（去重后拼接）
	result["B1_L54"] = "、".join(sorted(set([level_dict[i] for i in brca1_level_list if i in ["4","5"]]), reverse=True))
	result["B2_L54"] = "、".join(sorted(set([level_dict[i] for i in brca2_level_list if i in ["4","5"]]), reverse=True))
	# L543提取345，列表（去重），用于判断描述中用"有"还是"为"，仅1个用"为"，多个用"有"
	result["B1_L543"] = list(set([level_dict[i] for i in brca1_level_list if i in ["3","4","5"]]))
	result["B2_L543"] = list(set([level_dict[i] for i in brca2_level_list if i in ["3","4","5"]]))

	return result

def varInfo_FJZL(var, tumor_list,config):
	'''
	福建肿瘤，变异解读第二段，需要展示变异发生频率
	第四段，需要合并治疗方案
	'''
	# 获取变异发生频率信息
	#fjzl_database = get_FJZL_database(config)
	fdzs_dict, fjzl_database, io_dict = getconfigxlsx(config)
	var_freq_info = {}
	for gene in re.split(",", var["gene_symbol"]):
		for tumor in tumor_list:
			# 若匹配到基因 癌种，则使用该信息
			if (gene, tumor) in fjzl_database.keys():
				if gene not in var_freq_info.keys():
					var_freq_info.setdefault(gene, "")
				var_freq_info[gene] = fjzl_database[(gene, tumor)]
				
		# 若未匹配到基因 癌种，则匹配基因
		if gene not in var_freq_info.keys() and (gene, "") in fjzl_database.keys():
			var_freq_info.setdefault(gene, "")
			var_freq_info[gene] = fjzl_database[(gene, "")]
	
	# 治疗方案小结处理
	# 将药物按敏感/耐药，等级进行拆分，便于后续汇总
	regimen_info = {}
	def getinfo_S(regimen_list, conclusion_list):
		return [a["regimen_name"] for a in regimen_list if a["evi_conclusion_simple"] in conclusion_list]

	for i in ["A", "B", "C", "D"]:
		regimen_info["regimen_"+str(i)+"_S"] = getinfo_S(var["evi_sum"]["regimen_S"], [i])
		regimen_info["regimen_"+str(i)+"_R"] = getinfo_S(var["evi_sum"]["regimen_R"], [i])
	
	return var_freq_info, regimen_info

def BRCA_FJFY_summary(var_data):
	'''
	福建医科附一BRCA，报告中需要展示I/II类变异详情-2022.11.11
	*BRCA全血
		5类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，基因型为纯合/杂合（丰度为XX），对XX、XX可能敏感/可能耐药（A级）
		4类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，基因型为纯合/杂合（丰度为XX），对XX、XX可能敏感/可能耐药（A级）
		3类：待定
	*BRCA组织
		I类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，突变丰度为XX，对XX、XX可能敏感/可能耐药（A级）
		II类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，突变丰度为XX
		III类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，突变丰度为XX

	snvindel：XX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为X
	mlpa：XX基因value检测到大片段缺失/大片段重复
	'''
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}
		if var["type"] == "Loss":
			var_info = var["gene_symbol"]+"基因"+var["value"]+"检测到大片段缺失"
		elif var["type"] == "Gain":
			var_info = var["gene_symbol"]+"基因"+var["value"]+"检测到大片段重复"
		elif var["bio_category"] == "Snvindel":
			region_list_en = re.split("_", var["gene_region"])
			region_list_cn = []
			for i in region_list_en:
				if re.search("exon", i):
					region_list_cn.append(i.replace("exon", "")+"号外显子")
				elif re.search("intron", i):
					region_list_cn.append(i.replace("intron", "")+"号内含子")
				else:
					region_cn = region_dict[i] if i in region_dict.keys() else i
					region_list_cn.append(region_cn)
			if var["hgvs_p"] != "p.?":
				if var["var_origin"] == "germline":
					gene_type = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_p"]+"，基因型为"+gene_type+"（丰度为"+var["freq_str"]+"）"
				else:
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_p"]+"，突变丰度为"+str(var["freq_str"])
			else:
				if var["var_origin"] == "germline":
					gene_type = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_c"]+"，基因型为"+gene_type+"（丰度为"+var["freq_str"]+"）"
				else:
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到"+ var["type_cn"] + var["hgvs_c"]+"，突变丰度为"+str(var["freq_str"])
		return var_info

	def regimen_info(var):
		sense_regimen = []
		resis_regimen = []
		if var["evi_sum"]:
			regimen_info = ""
			#print (var["evi_sum"])
			sense_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] \
							if re.search("Sensitive", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] \
							if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []
			resis_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] \
							if re.search("Resistant", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] \
							if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []

		sense_regimen_str = "对" + "、".join(sense_regimen) + "可能敏感（A级）" if sense_regimen else ""
		resis_regimen_str = "对" + "、".join(resis_regimen) + "可能耐药（A级）" if resis_regimen else ""

		return sense_regimen_str, resis_regimen_str

	# 单个位点总结
	def var_inter(var):
		var_info = var_info_stran(var)
		sense_regimen_str, resis_regimen_str = regimen_info(var)
		tmp_list = []
		tmp_list.append(var_info)
		if sense_regimen_str:
			tmp_list.append(sense_regimen_str)
		if resis_regimen_str:
			tmp_list.append(resis_regimen_str)
		return "，".join(tmp_list)

	result = {}
	# BRCA全血
	result["gBRCA_5"] = "；".join([var_inter(var) for var in var_data["snv_s"]["B1_L5"] + var_data["snv_s"]["B2_L5"]])
	result["gBRCA_4"] = "；".join([var_inter(var) for var in var_data["snv_s"]["B1_L4"] + var_data["snv_s"]["B2_L4"] + var_data["mlpa"]["B1_Loss"]+var_data["mlpa"]["B2_Loss"]])
	result["gBRCA_3"] = "；".join([var_inter(var) for var in var_data["snv_s"]["B1_L3"] + var_data["snv_s"]["B2_L3"] + var_data["mlpa"]["B1_Gain"]+var_data["mlpa"]["B2_Gain"]])
	# BRCA组织
	result["tBRCA_I"] = "；".join([var_inter(var) for var in var_data["snv_s"]["B1_L5"] + var_data["snv_s"]["B2_L5"]+ var_data["snv_s"]["B1_L4"] + var_data["snv_s"]["B2_L4"] if var["clinic_num_s"] == 5])
	result["tBRCA_II"] = "；".join([var_inter(var) for var in var_data["snv_s"]["B1_L5"] + var_data["snv_s"]["B2_L5"]+ var_data["snv_s"]["B1_L4"] + var_data["snv_s"]["B2_L4"] if var["clinic_num_s"] == 4])
	result["tBRCA_III"] = "；".join([var_inter(var) for var in var_data["snv_s"]["B1_L3"] + var_data["snv_s"]["B2_L3"]])
		
	return result

def HRR_FJFY_summary(var_data):
	'''
	福建医科附一HRR，报告中需要展示I/II类变异详情-2022.11.17
	*HRR全血
		BRCA1/BRCA2 5类：XXX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，基因型为纯合/杂合，对XX、XX可能敏感/可能耐药（A级）
		BRCA1/BRCA2 4类：XXX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，基因型为纯合/杂合，对XX、XX可能敏感/可能耐药（A级）
		BRCA1/BRCA2 3类：XXX基因X号外显子/内含子检测出XX突变hgvs_p/hgvs_c，基因型为纯合/杂合
		其他基因 5/4/3类：仅展示是否有检测到，无需详细信息
	*HRR组织
		BRCA1/BRCA2 I类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为XX，对XX、XX可能敏感/可能耐药（A级）
		BRCA1/BRCA2 II类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，突变丰度为XX
		BRCA1/BRCA2 III类+肿瘤发生发展相关：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，突变丰度为XX
		其他基因 I/II/III类+肿瘤发生发展相关：仅展示是否有检测到，无需详细信息
	*HRR配对
	*（只需开发这个，单独全血或组织的可用配对的结果）
		BRCA1/BRCA2 5类：XXX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，基因型为纯合/杂合，对XX、XX可能敏感/可能耐药（A级）
		BRCA1/BRCA2 4类：XXX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，基因型为纯合/杂合，对XX、XX可能敏感/可能耐药（A级）
		BRCA1/BRCA2 I类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为XX，对XX、XX可能敏感/可能耐药（A级）
		BRCA1/BRCA2 II类：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，突变丰度为XX
		BRCA1/BRCA2 3类：XXX基因X号外显子/内含子检测出XX突变hgvs_p/hgvs_c，基因型为纯合/杂合
		BRCA1/BRCA2 III类+肿瘤发生发展相关：BRCAX基因X号外显子/内含子检测到XX突变hgvs_p(P.XX)/hgvs_c，突变丰度为XX
		其他基因 5/4/3类和I/II/III+肿瘤发生发展相关：仅展示是否有检测到，无需详细信息

	snvindel：XX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为X
	mlpa：XX基因value检测到大片段缺失/大片段重复
	'''
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}
		if var["type"] == "Loss":
			var_info = var["gene_symbol"]+"基因"+var["value"]+"检测到大片段缺失"
		elif var["type"] == "Gain":
			var_info = var["gene_symbol"]+"基因"+var["value"]+"检测到大片段重复"
		elif var["bio_category"] == "Snvindel":
			region_list_en = re.split("_", var["gene_region"])
			region_list_cn = []
			for i in region_list_en:
				if re.search("exon", i):
					region_list_cn.append(i.replace("exon", "")+"号外显子")
				elif re.search("intron", i):
					region_list_cn.append(i.replace("intron", "")+"号内含子")
				else:
					region_cn = region_dict[i] if i in region_dict.keys() else i
					region_list_cn.append(region_cn)
			if var["hgvs_p"] != "p.?":
				if var["var_origin"] == "germline":
					gene_type = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_p"]+"，基因型为"+gene_type
				else:
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_p"]+"，突变丰度为"+str(var["freq_str"])
			else:
				if var["var_origin"] == "germline":
					gene_type = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_c"]+"，基因型为"+gene_type
				else:
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到"+ var["type_cn"] + var["hgvs_c"]+"，突变丰度为"+str(var["freq_str"])
		return var_info

	def regimen_info(var):
		sense_regimen = []
		resis_regimen = []
		if var["evi_sum"]:
			regimen_info = ""
			#print (var["evi_sum"])
			sense_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if re.search("Sensitive", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []
			resis_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if re.search("Resistant", i["clinical_significance"]) and i["evi_conclusion_simple"]=="A"] if var["evi_sum"]["evi_split"] and "Predictive" in var["evi_sum"]["evi_split"].keys() else []

		sense_regimen_str = "对" + "、".join(sense_regimen) + "可能敏感（A级）" if sense_regimen else ""
		resis_regimen_str = "对" + "、".join(resis_regimen) + "可能耐药（A级）" if resis_regimen else ""

		return sense_regimen_str, resis_regimen_str

	# 单个位点总结
	def var_inter(var):
		var_info = var_info_stran(var)
		sense_regimen_str, resis_regimen_str = regimen_info(var)
		tmp_list = []
		tmp_list.append(var_info)
		if sense_regimen_str:
			tmp_list.append(sense_regimen_str)
		if resis_regimen_str:
			tmp_list.append(resis_regimen_str)
		return "，".join(tmp_list)

	result = {}
	# 胚系结果
	result["G_B_5"] = "；".join([var_inter(var) for var in var_data["var_germline"]["level_5"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel"])

	result["G_B_4"] = "；".join([var_inter(var) for var in var_data["var_germline"]["level_4"] + var_data["mlpa"]["B1_Loss"]+var_data["mlpa"]["B2_Loss"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and (var["bio_category"] == "Snvindel" or ("type" in var.keys() and var["type"] == "Loss"))])

	result["G_B_3"] = "；".join([var_inter(var) for var in var_data["var_germline"]["level_3"] + var_data["mlpa"]["B1_Gain"]+var_data["mlpa"]["B2_Gain"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and (var["bio_category"] == "Snvindel" or ("type" in var.keys() and var["type"] == "Gain"))])

	result["G_other"] = "；".join([var_inter(var) for var in var_data["var_germline"]["level_5"]+var_data["var_germline"]["level_4"]+var_data["var_germline"]["level_3"] \
															if var["gene_symbol"] not in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel"])
	# 体细胞结果
	result["S_B_I"] = "；".join([var_inter(var) for var in var_data["var_somatic"]["level_I"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])

	result["S_B_II"] = "；".join([var_inter(var) for var in var_data["var_somatic"]["level_II"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])

	result["S_B_III"] = "；".join([var_inter(var) for var in var_data["var_somatic"]["level_III"] + var_data["var_somatic"]["level_onco_nodrug"] \
															if var["gene_symbol"] in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])

	result["S_other"] = "；".join([var_inter(var) for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]+\
															 var_data["var_somatic"]["level_III"] + var_data["var_somatic"]["level_onco_nodrug"] \
															if var["gene_symbol"] not in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel" ])														

		
	return result

def varS3_splitfor301(var_data):
	'''
	将3类变异拆分为一行两个变异展示
	'''
	# 若为DNA/RNA共检融合，且RNA中有多个时，需要分开展示-2023.04.13
	var_list = []
	for var in var_data:
		var_list.append(var)
		if "dup_rna_detect" in var.keys() and len(var["dup_rna_detect"]) > 1 :
			for rna_var in var["dup_rna_detect"]:
				var_list.append(rna_var)
		#else:
		#	var_list.append(var)

	result = []
	if var_list:
		for i in range(0, len(var_list)-2, 2):
			tmp_dict = {}
			for j in range(1, 3):
				tmp_dict["var"+str(j)] = var_list[i+j-1]
			result.append(tmp_dict)

		rest_num= len(var_list) % 2
		rest_tmp_dict = {}
		for j in range(1, 3):
			rest_tmp_dict["var"+str(j)] = ""

		num = 1
		last_row_num = len(var_list)-rest_num if rest_num != 0 else len(var_list)-rest_num-2
		for j in range(last_row_num, len(var_list)):
			rest_tmp_dict["var"+str(num)] = var_list[j]
			num += 1
		result.append(rest_tmp_dict)
	return result

def sort_var_for301(var_data):
	'''
	301定制MP，检测结果小结和解读部分变异排序
	排序： 按基因、变异等级（I、II）、类型、频率（频率在各部分处理的时候已经排序了）
	10基因排前面
	'''
	LC10_gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	for var in var_data:
		gene = re.split(",", var["gene_symbol"])
		var["sort_301"] = 0 if set(gene) & set(LC10_gene_list) else 1
	return sorted(var_data, key=lambda i : (i["sort_301"], i["gene_symbol"]))

def sort_var_for301_2(var_data):
	'''
	301定制MP，排序规则更新
	1. 按基因（按I/II/III类变异排，变异等级越高的基因排在前面，等级相同的按首字母排）
	2. 点突变、融合、CNV
	'''
	# 1. 初步排序，根据I/II/III类；点突变、融合、CNV
	var_data_copy = copy.deepcopy(var_data)
	type_rule = {"Snvindel" : 0, "Sv" : 1, "PSeqRnaSv" : 2, "Cnv" : 3}
	for var in var_data_copy:
		# s_level, I类对应5，II类对应4，肿瘤发生发展相关对应3，用于报告填充
		var["s_level"] = S_level(var)
		# 证据最高等级转化为数字，便于排序，I类标记为0，II类标记为1， 肿瘤发生发展相关标记为2
		var["top_level_forsort"] = 0 if var["top_level"] in ["A", "B"] else 1 if var["top_level"] in ["C", "D"] else 2
		# sv的gene_symbol，若两个基因都在检测范围里，则返回five_prime_gene,three_prime_gene (目前返回的是three_prime_gene,five_prime_gene)
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			var["gene_symbol"] = var["five_prime_gene"]+","+var["three_prime_gene"]
	var_data_copy_sort = sorted(var_data_copy, key=lambda i : (i["top_level_forsort"], type_rule.get(i["bio_category"])))

	# 2. 按基因进行拆分，便于后面组装
	var_split = {}
	for var in var_data_copy_sort:
		if var["gene_symbol"] not in var_split.keys():
			var_split.setdefault(var["gene_symbol"], [])
		var_split[var["gene_symbol"]].append(var)

	# 3. 获取基因对应变异的最高等级，按顺序获取基因列表
	gene_level_sort = sorted([{"gene_symbol" : var["gene_symbol"], "top_level_forsort" : var["top_level_forsort"]} for var in var_data_copy],\
							 key=lambda i:(i["top_level_forsort"], i["gene_symbol"]))
	gene_list = []
	for i in gene_level_sort:
		if i["gene_symbol"] not in gene_list:
			gene_list.append(i["gene_symbol"])
	
	# 4. 根据已排序的基因列表，将结果进行组装
	result = []
	for gene in gene_list:
		result.extend(var_split[gene])

	return result

# 这个考虑可在报告模板中实现
def get_summary_GDRM_gHRR(var_data, somatic_var_data):
	'''
	用于广东人民gHRR检测结果总结
	基于本次送检样本，检测到XXX基因gene_region hgvs_c hgvs_p致病/疑似致病性变异。
	2023.04.26新增体细胞/来源不明变异：基于本次送检样本，检测到XXX基因gene_region hgvs_c hgvs_p I类/II类变异
	'''
	clinic_num_g_stran = {5 : "致病性变异", 4 : "疑似致病性变异"}
	var_list = []
	for var in var_data:
		if var["type"] == "Loss":
			var_list.append("{0}基因{1} del {2}".format(var["gene_symbol"], var["value"], "疑似致病性变异")) 
		else:
			if var["hgvs_p"] != "p.?":
				var_list.append("{0}基因{1} {2} {3}{4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"],clinic_num_g_stran.get(var["clinic_num_g"])))
			else:
				var_list.append("{0}基因{1} {2}{3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"],clinic_num_g_stran.get(var["clinic_num_g"])))
		
	# 2023.04.26新增体细胞/来源不命名变异总结
	clinic_num_s_stran = {5 : "I类变异", 4 : "II类变异"}
	somatic_var_list = []
	for var in somatic_var_data:
		if var["hgvs_p"] != "p.?":
			somatic_var_list.append("{0}基因{1} {2} {3} {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"],clinic_num_s_stran.get(var["clinic_num_s"])))
		else:
			somatic_var_list.append("{0}基因{1} {2} {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"],clinic_num_s_stran.get(var["clinic_num_s"])))

	return ";".join(var_list), ";".join(somatic_var_list), ";".join(somatic_var_list+var_list)

def PAN116_split_gene(var_data):
	'''
	厦门市一116，报告中需要展示成10gene+其他基因（只包含I/II类变异）以及胚系4/5类变异
	special["var_pan116_split"] = {
		gene10_level_I_II : [],
		other_level_I_II : [],
		gene10_germline_45 : [],
		other_germline_45 : []
	}
	'''
	gene10_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
#	gene10_level_I_II = [var for var in var_data if set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==5] + [var for var in var_data if set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==4]
#	gene30_level_I_II = [var for var in var_data if not set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==5] + [var for var in var_data if not set(re.split(",",var["gene_symbol"])) & set(gene10_list) and var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"]==4]
	# 新增III类-适用复旦中山厦门医院
	gene10_level_I_II = [var for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_I_II = [var for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene10_level_III = [var for var in var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"] if set(re.split(",",var["gene_symbol"])) & set(gene10_list)]
	gene30_level_III = [var for var in var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"] if not set(re.split(",",var["gene_symbol"])) & set(gene10_list)]

	gene10_germline_45 = [var for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"] if set(re.split(",", var["gene_symbol"])) & set(gene10_list)]
	other_germline_45 = [var for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"] if not set(re.split(",", var["gene_symbol"])) & set(gene10_list)]


	return gene10_level_I_II, gene30_level_I_II, gene10_level_III, gene30_level_III, gene10_germline_45, other_germline_45

# 统计体细胞/来源未明，非多态性变异数-适用郑大一HRR-2023.07.03
def varnum_ZDY_HRR(var_list):
	num = 0
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			var["tag"] = var["tag"] if var["tag"] else ""
			if not re.search("Polymorphism", var["tag"]):
				num += 1
	return num

# 上海仁济CP40，重要靶向基因结果汇总，要按固定的基因列表和基因顺序展示-2023.09.18
# summary: 仅返回结果汇总，包含I/II/III类变异
# var_inter：返回I/II类 + 其他基因I/II类结果
# gene_inter：返回I/II类 + 其他基因I/II类，需要去重
def shrj_cp40_var(var_list, config):
	var_dict = {
		"KRAS" : "点突变/插入/缺失",
		"NRAS" : "点突变/插入/缺失",
		"BRAF" : "点突变/插入/缺失",
		"POLE" : "点突变/插入/缺失",
		"ERBB2" : "点突变/插入/缺失/拷贝数变异",
		"NTRK1" : "点突变/插入/缺失/融合",
		"NTRK2" : "点突变/插入/缺失/融合",
		"NTRK3" : "点突变/插入/缺失/融合",
	}
	gene_list = ["KRAS", "NRAS", "BRAF", "POLE", "ERBB2", "NTRK1", "NTRK2", "NTRK3"]
	gene_info = get_shrj_geneinfo(config)
	detect_gene = []

	summary_result = []
	var_inter = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			var["level"] = str(int(S_level(var)))
			#var["gene_inter_shrj"] = gene_info.get(var["gene_symbol"])
			var["gene_num"] = gene_list.index(gene) if gene in gene_list else 999
			# summary添加
			if gene in var_dict.keys():
				detect_gene.append(gene)
				summary_result.append(var)
			# I/II类变异汇总
			if var["level"] in ["5", "4"]:
				if var not in var_inter:
					var_inter.append(var)

	# 未检出变异的重要基因
	for gene in set(gene_list) - set(detect_gene):
		summary_result.append({
			"gene_symbol" : gene,
			"var_type" : var_dict[gene],
			"gene_num" : gene_list.index(gene)
		})

	# 排序
	summary_result_sort = sorted(summary_result, key=lambda i:i["gene_num"])
	var_inter_sort = sorted(var_inter, key=lambda i:i["gene_num"])
	gene_inter_sort = []
	# 新增一个带gene的基因介绍-2023.09.28
	gene_inter_sort2 = []
	# 新增完成-2023.09.28
	for var in var_inter_sort:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene_info.get(gene) not in gene_inter_sort:
				gene_inter_sort.append(gene_info.get(gene))
		#if var["gene_inter_shrj"] not in gene_inter_sort:
		#	gene_inter_sort.append(var["gene_inter_shrj"])
			# 新增一个带gene的基因介绍-2023.09.28
			if {"gene" : gene, "gene_info" : gene_info.get(gene)} not in gene_inter_sort2:
				gene_inter_sort2.append({
					"gene" : gene,
					"gene_info" : gene_info.get(gene)
				})
			# 新增完成-2023.09.28

	return summary_result_sort, var_inter_sort, gene_inter_sort, gene_inter_sort2

# XW0417
def sort_for_xw0417(var_list):
	result = {
		'BICC1' : '',
		'CASP7' : '',
		'TACC3v1' : '',
		'TACC3v3' : '',
		'BAI' : '',
		'R248' : '',
		'G370' : '',
		'S249' : '',
		'Y373' : ''
	}
	for var in var_list:
		if var['bio_category'] in ['Sv', 'PSeqRnaSv']:
			if var['gene_symbol'] == 'FGFR2' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'BICC1' and var['three_prime_cds'] == 'exon3':
				result["BICC1"] = 'BICC1'
			elif var['gene_symbol'] == 'FGFR2' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'CASP7' and var['three_prime_cds'] == 'exon2':
				result['CASP7'] = 'CASP7'
			elif var['gene_symbol'] == 'FGFR3' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'TACC3' and var['three_prime_cds'] == 'exon11':
				result['TACC3v1'] = 'v1'
			elif var['gene_symbol'] == 'FGFR3' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'TACC3' and var['three_prime_cds'] == 'exon10':
				result['TACC3v3'] = 'v3'
			elif var['gene_symbol'] == 'FGFR3' and var['five_prime_cds'] == 'exon17' and var['three_prime_gene'] == 'BAIAP2L1' and var['three_prime_cds'] == 'exon2':
				result['BAI'] = 'BAI'
		
#		if var['gene_symbol'] == 'FGFR3' and var['bio_category'] == 'Snvindel':
#			if var['hgvs_c'] == 'c.742C>T':
#				result['R248'] = 'T'
#			elif var['hgvs_c'] == 'c.1108G>T':
#				result['G370'] = 'T'
#			elif var['hgvs_c'] == 'c.746C>G':
#				result['S249'] = 'T'
#			elif var['hgvs_c'] == 'c.1118A>G':
#				result['Y373'] = 'T'
		# V4更新
		if var['gene_symbol'] == 'FGFR3' and var['bio_category'] == 'Snvindel':
			if var['hgvs_p'] == 'p.(R248C)':
				result["R248"] = var['hgvs_c']
			elif var['hgvs_p'] == 'p.(G370C)':
				result["G370"] = var['hgvs_c']
			elif var['hgvs_p'] == 'p.(S249C)':
				result["S249"] = var['hgvs_c']
			elif var['hgvs_p'] == 'p.(Y373C)':
				result["Y373"] = var['hgvs_c']
	
	tmp_dict = {key:value for key, value in result.items() if value}
	result['final'] = 'T' if tmp_dict else 'F'

	return result

def judgeLBD(var_list):
	conclusion = []
	var_list = [var for var in var_list if var['gene_symbol'] == 'AR']
	if not var_list:
		return ''
	cnv = [var for var in var_list if var['bio_category'] == 'Cnv']
	#snv = [var for var in var_list if var['bio_category'] == 'Snvindel' and var['type'] == 'nonSynonymous_Substitution']
	snv = [var for var in var_list if var['bio_category'] == 'Snvindel']
	if cnv:
		for var in cnv:
			if var['cnv_type'] == 'Loss':
				conclusion.append('AR Loss')
			else:
				conclusion.append('AR Amp')
	if snv:
		lbd_present = False
		non_lbd_present = False
		for var in snv:
			# 646-919 之间的错义突变为LBD ，其他变异为非LBD突变，p.(V666M)
			num = re.findall(r'\d+', var['hgvs_p']) if var['hgvs_p'] != 'p.?' else []
			if not num:
				non_lbd_present = True
			elif 646 <= int(num[0]) <= 919 and var['type'] == 'nonSynonymous_Substitution':
				lbd_present = True
			else:
				non_lbd_present = True
		if lbd_present:
			conclusion.append('AR（LBD突变）')
		elif non_lbd_present:
			conclusion.append('AR（非LBD突变）')
	seen = set()
	unique_conclusion = []
	for item in conclusion:
		if item not in seen:
			seen.add(item)
			unique_conclusion.append(item)
	
	return '、'.join(unique_conclusion) if unique_conclusion else ''

def _process_hrr_genes(var_list: list, hrr_genelist: list) -> dict:
    """
    Helper function to process HRR genes from variant list.
    
    Args:
        var_list: List of variant dictionaries
        hrr_genelist: List of HRR genes to check against
        
    Returns:
        Dictionary containing processed gene lists
    """
    hrr_genes = set()
    other_genes = set()
    
    for var in var_list:
        genes = set(re.split(",", var["gene_symbol"]))
        hrr_genes.update(genes & set(hrr_genelist))
        other_genes.update(genes - set(hrr_genelist))
    ar_lbd = judgeLBD(var_list) if 'AR' in hrr_genes else ''

    result = {
        'AR': ar_lbd,
        'HRR': '、'.join(sorted(hrr_genes - {'AR'})) if hrr_genes - {'AR'} else '',
        'other': '、'.join(sorted(other_genes)) if other_genes else ''
    }
    
    return result

def getSum_for_tXW6002(var_list: list) -> dict:
    """
    Process tissue XW6002 HRR genes.
    
    Args:
        var_list: List of variant dictionaries
        
    Returns:
        Dictionary containing AR, HRR and other gene information
    """
    t_XW6002_HRR_genelist = [
        "ABRAXAS1", "AKT1", "AKT2", "AR", "ATM", "ATR", "ATRX", "AURKA", 
        "BARD1", "BRAF", "BRCA1", "BRCA2", "BRIP1", "CD274", "CHEK1", 
        "CHEK2", "EGFR", "ERBB2", "ERCC1", "FANCA", "FANCC", "FANCD2", 
        "FANCL", "FGFR1", "GEN1", "GSTP1", "IFNGR1", "KRAS", "MRE11", 
        "MSH6", "MYB", "NBN", "PALB2", "POLD1", "POLE", "PRKACA", "RAD50", 
        "RAD51B", "RAD51C", "RAD51D", "RAD54L", "RBM10", "TP53", "TSC1", 
        "TSC2", "VEGFA"
    ]
    return _process_hrr_genes(var_list, t_XW6002_HRR_genelist)

def getSum_for_gXW6002(var_list: list) -> dict:
    """
    Process germline XW6002 HRR genes.
    
    Args:
        var_list: List of variant dictionaries
        
    Returns:
        Dictionary containing AR, HRR and other gene information
    """
    g_XW6002_HRR_genelist = [
        "ABRAXAS1", "AKT1", "AKT2", "AR", "ATM", "ATR", "AURKA", "BARD1", 
        "BRAF", "BRCA1", "BRCA2", "BRIP1", "CD274", "CHEK1", "CHEK2", 
        "EGFR", "ERBB2", "ERCC1", "ERCC3", "ERCC4", "FANCA", "FANCD2", 
        "FANCL", "FANCM", "FGFR1", "GEN1", "KRAS", "MAPK1", "MRE11", 
        "MSH6", "NBN", "NPM1", "PALB2", "POLD1", "POLE", "RAD50", "RAD51", 
        "RAD51B", "RAD51C", "RAD51D", "RAD52", "RAD54L", "SETD2", "SLX4", 
        "TP53", "TSC1", "TSC2", "XRCC1", "XRCC2"
    ]
    return _process_hrr_genes(var_list, g_XW6002_HRR_genelist)

def getSum_for_XW5902(var_list, hd):
	result = ''
	if not var_list and not hd:
		result = '未检出'
	elif var_list and hd:
		result = '检出SNVIndel/HD'
	elif hd:
		result = '检出HD'
	elif var_list:
		result = '检出SNVIndel'
	
	return result

def var_ad3101_summary(var):
	somatic = var['var_somatic']
	germline = var['var_germline']
	variants = [
		*somatic['level_I'],
		*germline['level_5'],
		*somatic['level_II'],
		*germline['level_4'],
		*somatic['level_onco_nodrug'],
		*somatic['level_III'],
		*germline['level_3']
	]

	counter = defaultdict(int)
	for item in variants:
		g = item.get('clinic_ad3101_g', 0)
		s = item.get('clinic_ad3101_s', 0)

		if g >= 3 and not s:
			key = f'g{g}'
		elif s >= 3:
			key = f's{s}'
		else:
			continue
		counter[key] += 1

	return {
		'all' : len(variants),
		'num_5' : counter.get('g5', 0) + counter.get('s5', 0),
		'num_4' : counter.get('g4', 0) + counter.get('s4', 0),
		'num_3' : counter.get('g3', 0) + counter.get('s3', 0)
	}

def XW6701_summary(var_list, main_gene, hrr_genelist):
	'''
	用于XW6701检测结果总结
	'''
	result = {'BRCA': [], 'CP': [],'hrr': [], 'other': []}
	allowed_bio = ['Snvindel', 'Cnv', 'Sv', 'PSeqRnaSv', 'PHd']
	seen = {'BRCA': set(), 'hrr': set(), 'other': set(), 'CP': set()}  # 跟踪已添加条目，避免重复
	clin_map = {'Pathogenic': '致病性变异', 'Likely pathogenic': '疑似致病性变异'}
	func_map = {'Oncogenic': '致癌性变异', 'Likely oncogenic': '疑似致癌性变异'}
	sort_order = {
        '致病性变异': 0,
        '致癌性变异': 1,
        '疑似致病性变异': 2,
        '疑似致癌性变异': 3
    }

	def get_var_info(item):
		bio_category = item['bio_category']
		if bio_category == 'Snvindel':
			return f"{item['hgvs_c']} {item['hgvs_p']}"
		elif bio_category == 'Cnv':
			return item['cnv_type']
		elif bio_category in ['Sv', 'PSeqRnaSv']:
			return f"{item['five_prime_gene']}-{item['three_prime_gene']}"
		elif bio_category == 'PHd':
			return "HD"
		else:
			return ''

	for item in var_list:
		gene_symbol = item['gene_symbol']
		bio_category = item['bio_category']
		clin_sig = item['clinical_significance']
		func_class = item['function_classification']

		if bio_category not in allowed_bio:
			continue

		processed = False
		for category, genes in [('BRCA', main_gene), ('hrr', hrr_genelist), ('CP', ['CDK12', 'PALB2'])]:
			# 7601 hrr新增了PIK3CA，解读只返回了体系
			if gene_symbol in genes and (clin_sig in clin_map or func_class in func_map):
			    # if gene_symbol in genes and clin_sig in clin_map:
				inter = clin_map[clin_sig] if clin_sig != '-' else func_map[func_class]
				var_info = get_var_info(item)
				entry_key = (gene_symbol, inter, var_info)

				if entry_key not in seen[category]:
					seen[category].add(entry_key)
					result[category].append({
						'gene_symbol': gene_symbol,
						'inter': inter,
						'var_info': var_info
					})
				processed = True
				break
		if processed:
			continue
		
		inter = None
		if clin_sig != '-' and func_class == '-' and clin_sig in clin_map:
			inter = clin_map[clin_sig]
		elif func_class in func_map:
			inter = func_map[func_class]
		if inter:
			var_info = get_var_info(item)
			entry_key = (gene_symbol, inter, var_info)
			if entry_key not in seen['other']:
				seen['other'].add(entry_key)
				result['other'].append({
					'gene_symbol': gene_symbol,
					'inter': inter,
					'var_info': var_info
				})
	
	for category in result:
		result[category] = sorted(result[category], key=lambda x: (sort_order[x['inter']], x['gene_symbol'].upper()))
		
	return result

def getSum_for_tXW6701(var_list, hd, tumor_list):
	var = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug'] + hd if \
		"前列腺癌" in tumor_list else var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug']
	main_gene = ['BRCA1', 'BRCA2']
	hrr_genelist = [
    "ARID1A", "ARID1B", "ARID2", "ATM", "ATR", "BAP1", "BARD1","BRIP1", "CCND1", "CCNE1", "CDK4",
	"CHEK1", "CHEK2", "ERCC1", "ERCC2", "FANCA", "FANCC", "FANCD2","FANCL", "GEN1", "HDAC2", "MDM2", "MDM4", "MLH1", "MLH3",
    "MSH2", "MSH3", "MSH6", "MTOR", "MUTYH", "NBN","PBRM1", "PMS2", "POLD1", "POLE", "RAD50", "RAD51B", "RAD51C",
    "RAD51D", "RAD54L", "SMARCA2", "SMARCA4", "SMARCB1", "TERT"]
	result = XW6701_summary(var, main_gene, hrr_genelist)

	return result

def getSum_for_gXW6701(var_list, hd, tumor_list):
	var = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug'] + hd if \
		"前列腺癌" in tumor_list else var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug']
	main_gene = ['BRCA1', 'BRCA2']
	hrr_genelist = [
    "ARID1A", "ATM", "ATR", "BAP1", "BARD1", "BRIP1", "CCND1", "CCND2", "CCNE1", "CDK4", "CDKN2A",
    "CHEK1", "CHEK2", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "FANCA", "FANCD2", "FANCL", "FANCM", "GEN1", "HDAC2", "MDM2", "MDM4",
    "MLH1", "MSH2", "MSH6", "MTOR", "MUTYH", "NBN", "NPM1", "PMS2", "POLD1", "POLE", "RAD50", "RAD51", "RAD51B",
    "RAD51C", "RAD51D", "RAD52", "RAD54L", "SLX4", "SMARCA4", "TERT", "TP53", "XRCC1", "XRCC2"]
	result = XW6701_summary(var, main_gene, hrr_genelist)

	return result

def getSum_for_XW6701(var_list, hd, tumor_list):
	var = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug'] + hd if \
		"前列腺癌" in tumor_list else var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug']
	main_gene = ['BRCA1', 'BRCA2']
	hrr_genelist = ['ATM', 'BARD1', 'BRIP1', 'CHEK1', 'CHEK2', 'FANCA', 'FANCL', 'HDAC2',
				  'PPP2R2A', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD54L']
	result = XW6701_summary(var, main_gene, hrr_genelist)

	return result

def getSum_for_XW6003(var_list, hd):
	result = {'MTAP' : set(), 'CDKN2A' : set(), 'CDKN2B' : set(), 'other' : set()}
	target_genes = {'MTAP', 'CDKN2A', 'CDKN2B'}
	var_list = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug'] + var_list['var_somatic']['level_III']

	for var in var_list:
		gene = var['gene_symbol']
		if gene in target_genes:
			result[gene].add(var['bio_category'])
		else:
			result['other'].add(gene)
	if hd:
		for item in hd:
			gene = item['gene_symbol']
			if gene in target_genes:
				result[gene].add('HD（纯合缺失）')
	
	return {
		'MTAP': '/'.join(sorted(result['MTAP'])),
		'CDKN2A': '/'.join(sorted(result['CDKN2A'])),
		'CDKN2B': '/'.join(sorted(result['CDKN2B'])),
		'other': '、'.join(sorted(result['other']))
	}

# 临沂肿瘤医院116
def PAN116_LYZL_summary(level_I, level_II, level_onco_nodrug, level_5, level_4):
	'''
	临沂肿瘤医院-SF2
	首页展示A/B类证据相关变异
	EGFR基因19号外显子p.(L858R)突变 MET基因扩增 EML4-ALK基因融合
	检出XXX基因XX外显子p.(XXXXX)突变,为致病性变异;
	检出XXX基因XX外显子p.(XXXXX)突变,为疑似致病性变异。
	非配对样本，把配对样本检测的胚系基因筛选出来
	'''
	result = {
		"level_I_sum" : "",
		"nc_germline_sum" : "",
		"germline_sum" : ""
	}
	# 检测胚系变异的40基因
	gene_list = [
    "VHL", "TSC2", "TSC1", "TP53", "TERT", "STK11", "SMARCA4", "SMAD4", "RET",
    "RB1", "PTEN", "PTCH1", "POLE", "POLD1", "PMS2", "PDGFRA", "PALB2", "NF2",
    "NF1", "MSH6", "MSH2", "MLH1", "MET", "KIT", "HRAS", "FLCN", "FANCA", "EPCAM",
    "EGFR", "CDKN2A", "CDK4", "CDK12", "BRCA2", "BRCA1", "BRAF", "BAP1", "ATR",
    "ATM", "APC", "ALK"]

	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}
		if var["bio_category"] == "Cnv":
			var_info = var["gene_symbol"]+"基因拷贝数扩增"
		elif var["bio_category"] == "Sv":
			var_info = var["five_prime_gene"]+"-"+var["three_prime_gene"]+"基因融合"
		elif var["bio_category"] == "Snvindel":
			if var["gene_symbol"] == "MET" and "exon14 skipping" in var["variant_interpret_cn"]:
				if var["hgvs_p"] != "p.?":
					var_info  = "MET基因14号外显子跳跃突变"+var["hgvs_p"]
				else:
					var_info = "MET基因14号外显子跳跃突变"+var["hgvs_c"]
			else:
				region_list_en = re.split("_", var["gene_region"])
				region_list_cn = []
				for i in region_list_en:
					if re.search("exon", i):
						region_list_cn.append(i.replace("exon", "")+"号外显子")
					elif re.search("intron", i):
						region_list_cn.append(i.replace("intron", "")+"号内含子")
					else:
						region_cn = region_dict[i] if i in region_dict.keys() else i
						region_list_cn.append(region_cn)
				if var["hgvs_p"] != "p.?":
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn) + " " + var["hgvs_p"]+"突变"
				else:
					if var["type_cn"] != "--":
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn) + " " + var["hgvs_c"]+"突变"
					else:
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn) + " " + var["hgvs_c"]+"突变"
		return var_info

	# I类 具有临床意义的基因变异
	level_I_sum = [var_info_stran(var) for var in level_I] if level_I else []
	result["level_I_sum"] = "；\n".join(level_I_sum) if level_I_sum else ""
	# 非配对样本胚系变异
	var_list = level_I + level_II + level_onco_nodrug
	var_list = [var for var in var_list if var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"] if var_list else []
	var_list = sorted(var_list, key=lambda x: (x["clinic_num_g"], x["gene_symbol"]), reverse= True)
	if not var_list:
		result["nc_germline_sum"] = []
	else:
		nc_germline_sum = ["检出" + var_info_stran(var) + "，为致病性变异，待验证" if var["clinic_num_g"] == 5 else "检出" + var_info_stran(var) + "，为疑似致病性变异，待验证" for var in var_list]
		result["nc_germline_sum"] = "；\n".join(nc_germline_sum) if nc_germline_sum else ""
	# 配对样本胚系变异
	germline_sum = ["检出" + var_info_stran(var) + "，为致病性变异" if var["clinic_num_g"] == 5 else "检出" + var_info_stran(var) + "，为疑似致病性变异" for var in level_5 + level_4] if level_5 + level_4 else []
	result["germline_sum"] = "；\n".join(germline_sum) if germline_sum else ""

	return result

def getSum_for_gXW7601(var_list):
	var = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug']
	main_gene = ['BRCA1', 'BRCA2']
	hrr_genelist = ['AKT1', 'AR', 'ATM', 'ATR', 'BARD1', 'BRIP1', 'CDK12', 'CHEK1', 'CHEK2', 'FANCA', 'FANCL', 'HDAC2', 
				 'MLH1', 'MRE11', 'NBN', 'PALB2', 'PIK3CA', 'PTEN', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD54L']
	result = XW6701_summary(var, main_gene, hrr_genelist)
	return result

def getSum_for_tXW7601(var_list, hd, tumor_list):
	var = var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug'] + hd if \
		"前列腺癌" in tumor_list else var_list['var_somatic']['level_I'] + var_list['var_somatic']['level_II'] + var_list['var_somatic']['level_onco_nodrug']
	main_gene = ['BRCA1', 'BRCA2']
	hrr_genelist = ['ATM', 'BARD1', 'BRIP1', 'CDK12', 'CHEK1', 'CHEK2', 'FANCA', 'FANCL', 'HDAC2', 'PALB2', 'PPP2R2A', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD54L']
	result = XW6701_summary(var, main_gene, hrr_genelist)
	return result
