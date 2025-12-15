#-*- coding:gbk -*-
import re
from libs.processTheSnvindel import process_snvindel
from libs.processTheCNV import process_cnv
from libs.processTheSV import process_sv
from libs.processTheKNB import process_knb
from libs.processTheEC_type import process_ec_type
from libs.getInterRef import getRef_from_inter
from libs.NCCN_Gene_detect import getNCCN_detect
from libs.mergeMET14 import mergeVar
from libs.detectResultSummary import getSummary_detect
import copy
from libs.specialRequest import judge_GA_tumor_KNB
from libs.specialRequest import sv_shfk
#from libs.io_detect import io_detect
#from libs.io_detect import io_detect_for_116
#from libs.io_detect import io_detect_for_ZJZL
from libs.io_detect import get_io_detect
from libs.specialRequest import CP40_split_gene, PAN116_split_gene
from libs.specialRequest import CP40_FJFY_summary
from libs.processTheMLPA import process_mlpa
from libs.processTheGA_type import process_ga_type
from libs.specialRequest import Master_sum
from libs.getHRD_Master import getgss
from libs.rule import s_var_rule, s_var_rule_gene, s_var_rule_genelist, g_var_rule, g_var_regimen_rule, var_regimen_rule, g_var_rule_genelist
from libs.rule import process_result_regimen_ZJZL, var_bptm_rule, var_lyn5_rule, nofoundPath_genelist, nofound_genelist
from libs.specialRequest import HRR_FJFY_summary, varS3_splitfor301, sort_var_for301, sort_var_for301_2, get_summary_GDRM_gHRR
from libs.specialRequest import varInter_FJZL, var_summary_GDRM, var_GZZL, var_summary_SHRJ, var_summary_FJZL, var_SD
from libs.specialRequest import varInfo_XAJDY, var_summary_CQFY, BRCA_FJFY_summary, var_summary_WHXH, var_level123_SDZL, brca_judge_inApprovalTumor, varnum_ZDY_HRR
from libs.specialRequest import shrj_cp40_var
from libs.specialRequest import sort_for_xw0417
from libs.specialRequest import getSum_for_tXW6002
from libs.specialRequest import getSum_for_gXW6002
from libs.specialRequest import getSum_for_XW5902
from libs.specialRequest import var_ad3101_summary, getSum_for_tXW6701, getSum_for_gXW6701, getSum_for_XW6003, getSum_for_XW6701, PAN116_LYZL_summary, getSum_for_gXW7601, getSum_for_tXW7601

def getVar(jsonDict, config, report_name):
	data = {}
	snvindel = copy.deepcopy(process_snvindel(jsonDict, config))
	cnv = copy.deepcopy(process_cnv(jsonDict, config))
	sv_combination, rna_sv, rna_sv_only = copy.deepcopy(process_sv(jsonDict, config))
	mlpa, mlpa_image, mlpa_image_del = process_mlpa(jsonDict, config)

	# 更新为非肺癌，3类变异中也不展示BCL2L11-2022.09.06
	# 进院模板要展示BCL2L11，放3类，这边再加个临检的限制-2022.09.14
	# 这边再加个临检通用模板的限制
	# 2023.06.01-最新的临检116/76模板，附录最后不放BCL2L11了，该部分注释掉。如果报告中不想展示BCL2L11变异，可在数据审核步骤过滤掉。
	#if re.search("Pan116|LC76", jsonDict["sample_info"]["prod_names"]) and jsonDict["sample_info"]["report_module_type"] == "rummage" and re.search("rummage", report_name):
	#	for var in snvindel:
	#		if var["gene_symbol"] == "BCL2L11" and var["hgvs_c"] and var["hgvs_c"] == "c.394+1479_394+4381del":
	#			snvindel.remove(var)

	

	### 变异排序
	# 新增浙江肿瘤Master的变异排序，送检单位：浙江省肿瘤医院，产品：Master，业务类型：进院
	# 排序规则：按胚系/体细胞、I/II/肿瘤发生发展相关/III类、相同等级按治疗方案最高等级、snvindel>sv>cnv、频率降
	if re.search("浙江省肿瘤医院", jsonDict["sample_info"]["company"]) and re.search("Master|3231基因", jsonDict["sample_info"]["prod_names"]) and jsonDict["sample_info"]["report_module_type"] == "hospital":
		var_origin_rule = {"somatic" : 0, "germline" : 1}
		top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
		var_type_rule = {"Snvindel" : 0, "Sv" : 1, "PSeqRnaSv" : 2, "Cnv" : 3}
		var_data_without_rnasv = sorted(snvindel + cnv + sv_combination, key=lambda i : (var_origin_rule.get(i["var_origin"], 0), top_level_rule.get(i["top_level"]), var_type_rule.get(i["bio_category"])))
		var_data_rna_sv = sorted(rna_sv, key=lambda i : (var_origin_rule.get(i["var_origin"], 0), top_level_rule.get(i["top_level"]), var_type_rule.get(i["bio_category"])))
		var_data =  sorted(snvindel + cnv + sv_combination + rna_sv_only, key=lambda i : (var_origin_rule.get(i["var_origin"], 0), top_level_rule.get(i["top_level"]), var_type_rule.get(i["bio_category"])))
	
	# 其他模板都按下面的排序规则
	else:
		# 排序 按I、II、基因、类型、频率（频率在各部分处理的时候已经排序了）
		rule = {"Snvindel" : 0, "Cnv" : 1, "Sv" : 2, "PSeqRnaSv" : 3}
		# 更新排序规则，按I/II（同等级时按治疗方案最高等级进行排序）、基因、类型、频率（频率在各部分处理的时候已经排序了）
		# 同等级时按治疗方案最高等级进行排序 A > B > C > D > N
		top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
		#var_data_without_rnasv = sorted(snvindel + cnv + sv_combination, key=lambda i : (i["clinic_num_s"], top_level_rule.get(i["top_level"]), i["gene_symbol"], rule.get(i["bio_category"])))
		#var_data_rna_sv = sorted(rna_sv, key=lambda i : (i["clinic_num_s"], top_level_rule.get(i["top_level"]), i["gene_symbol"], rule.get(i["bio_category"])))
		#var_data = sorted(snvindel + cnv + sv_combination + rna_sv_only, key=lambda i : (i["clinic_num_s"], top_level_rule.get(i["top_level"]), i["gene_symbol"], rule.get(i["bio_category"])))
		# 更新下排序，把clinic_num_s排到top_level后面，并且clinic_num_s做个倒序转化-2023.02.01
		clinic_num_s_rule = {5 : 0, 4 : 1, 3 : 2, 2 : 3, 1 : 4}
		var_data_without_rnasv = sorted(snvindel + cnv + sv_combination, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"], rule.get(i["bio_category"])))
		var_data_rna_sv = sorted(rna_sv, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"], rule.get(i["bio_category"])))
		var_data = sorted(snvindel + cnv + sv_combination + rna_sv_only, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"], rule.get(i["bio_category"])))
	### 

	# 合并MET 14跳跃突变
	if jsonDict['sample_info']['product_name'] != 'XW1402':
		var_data_without_rnasv, judge_mergeMET_without_rnasv = mergeVar(var_data_without_rnasv)
		var_data, judge_mergeMET = mergeVar(var_data)

	# 体细胞结果整理
	data["var_somatic_without_rnasv"] = s_var_rule(var_data_without_rnasv)
	data["var_somatic_rna_sv"] = s_var_rule(var_data_rna_sv)
	data["var_somatic"] = s_var_rule(var_data)

	# 胚系结果整理
	data["var_germline"] = {**g_var_rule(var_data), **g_var_regimen_rule(var_data)}
	
	# 胚系+体细胞有用药建议/预后/辅助诊断变异整理
	data["var_for_regimen"] = var_regimen_rule(var_data)
	data["var_for_regimen_without_rnasv"] = var_regimen_rule(var_data_without_rnasv)

	# 浙江肿瘤胚系+体细胞有用药建议/预后/辅助诊断变异整理（报告中仅展示用药，无用药且有预后/辅助诊断的写“暂无靶向用药提示”，模板代码中较难实现，故这边再加字段）
	# 仅用在靶向用药提示部分，靶向用药解释还是可以用通用的var_for_regimen
	data["var_for_regimen_ZJZL"] = copy.deepcopy(data["var_for_regimen"])
	data["var_for_regimen_ZJZL"]["level_I"] = process_result_regimen_ZJZL(data["var_for_regimen_ZJZL"]["level_I"])
	data["var_for_regimen_ZJZL"]["level_II"] = process_result_regimen_ZJZL(data["var_for_regimen_ZJZL"]["level_II"])
	# 浙肿胚系致病/疑似致病但无靶向用药的位点解释放在5.1 临床意义明确的变异及靶向用药解析里，逻辑很怪，但还是得做，这边单独做一个列表展示吧
	data["var_germline_nodrug"] = [var for var in var_data if var["clinic_num_g"] in [5, 4] and var["var_origin"] == "germline" and (("evi_sum" in var.keys() and not var["evi_sum"]["evi_split"]) or ("evi_sum" in var.keys() and var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys())))]

	# knb
	data["knb"] = copy.deepcopy(process_knb(jsonDict, config))

	# 子宫内膜癌，其中BRCA1和BRCA2的结果还可以用在HRD上
	data["ec_type"] = {}
	if jsonDict["ec_type"]:
		data["ec_type"] = copy.deepcopy(process_ec_type(jsonDict, config))
		if data["ec_type"]["evi_sum"]["evi_split"] and "Prognostic" in data["ec_type"]["evi_sum"]["evi_split"].keys():
			# 预后 目前只有better和poor，后续会改成三个，列表
			data["ec_type"]["clinical_significance"] = [i["clinical_significance"] for i in data["ec_type"]["evi_sum"]["evi_split"]["Prognostic"]]
			# 证据描述，合并为字符串
			data["ec_type"]["evi_interpretation"] = "".join([i["evi_interpretation"] for i in data["ec_type"]["evi_sum"]["evi_split"]["Prognostic"]])
			# 目前CP40中的临床证据，仅显示“结果显示，”后的内容
			data["ec_type"]["evi_inter_cp40"] = re.split("结果显示，", data["ec_type"]["evi_interpretation"] )[-1] if data["ec_type"]["evi_interpretation"] else ""
	# TP53、POLE、BRCA1、BRCA2基因检测结果，分为I/II类和III类，用于变异和分子分型分开展示的情况
	ec_gene_list = ["TP53", "POLE", "BRCA1", "BRCA2"]
	for gene in ec_gene_list:
		data["ec_type"].update(var_bptm_rule(var_data, gene))

	# 新增BPTM（全血）的情况，该项目var_origin为germline, 不适用体细胞的-20221027
	# BRCA可以用BRCA项目的数据结构，这边仅处理POLE和TP53
	data["ec_type"]["POLE_g"] = g_var_rule_genelist(var_data, ["POLE"])
	data["ec_type"]["TP53_g"] = g_var_rule_genelist(var_data, ["TP53"])
	# 新增结束-20221027

	# MLPA结果
	#data["mlpa"], data["mlpa_image"], data["mlpa_image_del"] = process_mlpa(jsonDict, config)
	data["mlpa"] = mlpa
	data["mlpa_image"] = mlpa_image
	data["mlpa_image_del"] = mlpa_image_del

	# 伴随诊断检测结果
	data["cdx"] = getNCCN_detect(var_data, jsonDict["sample_info"]["tumor_list"], data["mlpa"], config)

	# 检测结果列表
	#data["summary_result"] = getSummary_detect(var_data, data["mlpa"], config)

	# 免疫正负相关
	#data["io"] = {}
	#data["io"]["result"], data["io"]["io_p_summary"], data["io"]["io_n_summary"] = io_detect(var_data)
	#data["io"]["result_116"], data["io"]["io_p_summary_116"], data["io"]["io_n_summary_116"] = io_detect_for_116(var_data)
	#data["io"]["io_p_summary_ZJZL"], data["io"]["io_n_summary_ZJZL"] = io_detect_for_ZJZL(var_data, config)
	data["io"] = get_io_detect(var_data, config)

	# 胚系-林奇综合征：五个基因结果需要分开展示（EPCAM、MLH1、MSH2、MSH6、PMS2）
	# 展示上述基因的3、4、5类变异
	gLS5_gene_list = ["EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	data["gLS5"] = var_lyn5_rule(var_data, gLS5_gene_list)

	# CRC12-MSI需要展示未检测到I/II类变异的基因
	crc12_gene = ["KRAS", "NRAS", "BRAF", "POLE", "PIK3CA", "ERBB2", "ALK", "FGFR3", "NTRK1", "NTRK3", "RET"]
	data["crc12_withoutPathVar_geneList"] = nofoundPath_genelist(var_data, crc12_gene)

	# LC10需展示未检到I/II类变异的基因
	lc10_gene = ["BRAF", "EGFR", "ERBB2", "KRAS", "MET", "ALK", "ROS1", "RET", "NRAS", "PIK3CA"]
	data["lc10_withoutPathVar_geneList"] = nofoundPath_genelist(var_data, lc10_gene) 
	# LC10临检通用模板展示未检测到I/II/III类变异的基因
	data["lc10_withoutVar_geneList"] = nofound_genelist(var_data, lc10_gene)

	# 特殊需求
	data["special"] = {}
	# 1. 复旦中山CP40需判断消化道肿瘤中KNB检出情况
	#data["special"]["FDZS_GA_KNB"] = judge_GA_tumor_KNB(var_data, jsonDict["sample_info"]["tumor_list"])
	# 1. 复旦中山CP40需判断消化道肿瘤中KNB检出情况-更新为只匹配tumor_names_cn-20220927
	#data["special"]["FDZS_GA_KNB"] = judge_GA_tumor_KNB(var_data, [jsonDict["sample_info"]["tumor_names_cn"]])
	# 2. 上海肺科CP40融合分型展示
	#data["special"]["SHFK_SV"] = sv_shfk(var_data)
	# 3. 厦门市一、聊城CP40 I/II类变异分为10基因+30基因展示
	'''
	data["special"]["var_cp40_split"] = {}
	data["special"]["var_cp40_split"]["gene10_level_I_II"], data["special"]["var_cp40_split"]["gene30_level_I_II"], \
	data["special"]["var_cp40_split"]["gene10_level_III"], data["special"]["var_cp40_split"]["gene30_level_III"],\
	data["special"]["var_cp40_split"]["gene10_level_nj_onco_nodrug"], data["special"]["var_cp40_split"]["gene30_level_nj_onco_nodrug"],\
	data["special"]["var_cp40_split"]["gene10_level_nj_III"], data["special"]["var_cp40_split"]["gene30_level_nj_III"] = CP40_split_gene(data)
	# 4. 福建附一CP40 I、II类变异小结
	data["special"]["var_cp40_FJFY"] = {}
	#data["special"]["var_cp40_FJFY"]["level_I_sum"], data["special"]["var_cp40_FJFY"]["level_II_sum"] = CP40_FJFY_summary(data["var_somatic"]["level_I"], data["var_somatic"]["level_II"])
	data["special"]["var_cp40_FJFY"] = CP40_FJFY_summary(data["var_somatic"]["level_I"], data["var_somatic"]["level_II"],data["var_somatic"]["level_III"],data["var_somatic"]["level_onco_nodrug"])
	'''
	# 5. 统计Master体细胞格式（DNA/RNA分开统计）
	data["special"]["Master_level_I"], data["special"]["Master_level_II"], data["special"]["Master_level_onco_nodrug"], data["special"]["Master_level_III"] = Master_sum(var_data)
	# 6. 汇总I/I/III类变异的基因列表
	'''
	data["special"]["SYX_levelI"] = "、".join(sorted(set([var["gene_symbol"] for var in data["var_somatic"]["level_I"]])))
	data["special"]["SYX_levelII"] = "、".join(sorted(set([var["gene_symbol"] for var in data["var_somatic"]["level_II"]])))
	data["special"]["SYX_levelIII"] = "、".join(sorted(set([var["gene_symbol"] for var in data["var_somatic"]["level_III"] + data["var_somatic"]["level_onco_nodrug"]])))
	# 7. 胚系三类变异基因列表，包含BRCA dup
	data["special"]["gene_var3_list"] = "、".join(sorted(set([var["gene_symbol"] for var in (data["var_germline"]["level_3"]+data["mlpa"]["B1_Gain"]+data["mlpa"]["B2_Gain"])])))
	'''
	#print (data["special"]["gene_var3_list"])
	# 8. 返回DNA/RNA MET14共检字段
	data["special"]["judge_mergeMET"] = judge_mergeMET if jsonDict['sample_info']['product_name'] != 'XW1402' else ''
	# 9. 福建附一HRR总结
	'''
	data["special"]["FJFY_HRR_sum"] = HRR_FJFY_summary(data)
	# 10. 301 MP 三类变异，一行展示两个变异
	data["special"]["varS3_301"] = varS3_splitfor301(data["var_somatic"]["level_III"])
	# 10. 301 gHRR 三类变异，一行展示两个变异，还需考虑MLPA dup-2023.09.22
	data["special"]["varS3_301_germline"] = varS3_splitfor301(data["mlpa"]["B1_Gain"]+data["mlpa"]["B2_Gain"]+data["var_germline"]["level_3"])
	# 2023.09.22添加完成
	# 11. 301 MP 变异排序
	data["special"]["var_sort_301"] = {}
		# 11.1 有对照时的“具有临床意义的变异”
	data["special"]["var_sort_301"]["Match_drug"] = sort_var_for301(data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"])
		# 11.2 有对照时的“肿瘤发生发展相关变异”
	data["special"]["var_sort_301"]["Match_withoutdrug"] = sort_var_for301(data["var_somatic"]["level_onco_nodrug"])
		# 11.3 无对照时的“具有临床意义的变异”
	data["special"]["var_sort_301"]["single_drug"] = sort_var_for301(data["var_germline"]["regimen_level_I"] + data["var_somatic"]["level_I"] +\
																	data["var_germline"]["regimen_level_II"] + data["var_somatic"]["level_II"])
		# 11.4 无对照时的“肿瘤发生发展相关变异”
	data["special"]["var_sort_301"]["single_withoutdrug"] = sort_var_for301(data["var_germline_nodrug"] + data["var_somatic"]["level_onco_nodrug"])

	# 12. 301 MP 排序规则2
	data["special"]["var_sort_301_2"] = {}
	# 12.1 有对照时的“体细胞变异”
	data["special"]["var_sort_301_2"]["somatic"] = sort_var_for301_2(data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"] +\
																	data["var_somatic"]["level_onco_nodrug"])
	# 12.2 无对照时的“基因变异”
	data["special"]["var_sort_301_2"]["var"] = sort_var_for301_2(data["var_germline"]["regimen_level_I"] + data["var_somatic"]["level_I"] +\
																data["var_germline"]["regimen_level_II"] + data["var_somatic"]["level_II"]+\
																data["var_germline_nodrug"] + data["var_somatic"]["level_onco_nodrug"])
	
	# 13. 广东省人民gHRR结果小结
	data["special"]["summary_GDRM_gHRR"], data["special"]["summary_GDRM_tHRR"], data["special"]["summary_GDRM_ptHRR"] = get_summary_GDRM_gHRR([var for var in data["var_germline"]["level_5"] if var["bio_category"] == "Snvindel"]\
																+data["mlpa"]["B1_Loss"]\
																+data["mlpa"]["B2_Loss"]\
																+[var for var in data["var_germline"]["level_4"] if var["bio_category"] == "Snvindel"], \
																[var for var in data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"] if var["bio_category"] == "Snvindel"])
	# 14. 厦门市一 116 I/II类变异分为10基因+其他基因展示
	data["special"]["var_pan116_split"] = {}
	data["special"]["var_pan116_split"]["gene10_level_I_II"], \
	data["special"]["var_pan116_split"]["other_level_I_II"], \
	data["special"]["var_pan116_split"]["gene10_level_III"], \
	data["special"]["var_pan116_split"]["other_level_III"], \
	data["special"]["var_pan116_split"]["gene10_germline_45"], \
	data["special"]["var_pan116_split"]["other_germline_45"] = PAN116_split_gene(data)
	# 15. 统计体细胞/来源未明，非多态性变异数-适用郑大一HRR
	data["special"]["varnum_ZDY_HRR"] = varnum_ZDY_HRR(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_III"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_IV"])
	# 16. 上海仁济CP40-重要基因检测结果基因列表和顺序保持不变（I/II/III），变异解读和基因介绍为重要基因I/II类+其他基因I/II类，基因介绍需要去重-2023.09.18
	data["special"]["shrj_cp40"] = {}
	data["special"]["shrj_cp40"]["summary_cdx"],\
	data["special"]["shrj_cp40"]["var_inter"],\
	data["special"]["shrj_cp40"]["gene_inter"],\
	data["special"]["shrj_cp40"]["gene_inter2"] = shrj_cp40_var(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"], config)
	#print (data["special"]["shrj_cp40"]["var_inter"])
	'''
	# xw0417结果
	data['special']['xw0417'] = sort_for_xw0417(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"])
	data['special']['tXW6002'] = getSum_for_tXW6002(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"])
	data['special']['gXW6002'] = getSum_for_gXW6002(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"])
	data['special']['xw5902'] = getSum_for_XW5902(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"]+data["var_somatic"]["level_III"], jsonDict['hd'])
	data['special']['ad3101'] = var_ad3101_summary(data)
	data['special']['tXW6701'] = getSum_for_tXW6701(data, jsonDict['hd'], jsonDict['sample_info']['tumor_list'])
	data['special']['gXW6701'] = getSum_for_gXW6701(data, jsonDict['hd'], jsonDict['sample_info']['tumor_list'])
	data['special']['XW6701'] = getSum_for_XW6701(data, jsonDict['hd'], jsonDict['sample_info']['tumor_list'])
	data['special']['xw6003'] = getSum_for_XW6003(data, jsonDict['hd'])
	data['special']['xw0294'] = getSum_for_tXW6002(data["var_somatic"]["level_I"]+data["var_somatic"]["level_II"]+data["var_somatic"]["level_onco_nodrug"])
	data['special']['gXW7601'] = getSum_for_gXW7601(data)
	data['special']['tXW7601'] = getSum_for_tXW7601(data, jsonDict['hd'], jsonDict['sample_info']['tumor_list'])
	# 临沂肿瘤
	data["special"]["LYZL_116"] = PAN116_LYZL_summary(data["var_somatic"]["level_I"], data["var_somatic"]["level_II"], data["var_somatic"]["level_onco_nodrug"], data["var_germline"]["level_5"], data["var_germline"]["level_4"])

	# BCL2L11基因2号内含子胚系缺失多态性检测结果 
	data["BCL2L11"] = ""
	for var in jsonDict["snvindel"]:
		if var["gene_symbol"] == "BCL2L11" and var["hgvs_c"] == "c.394+1479_394+4381del":
			data["BCL2L11"] = "T"
			break

	# 胃癌分子分型结果
	data["GA_type"] = process_ga_type(jsonDict, var_data)

	# Master HD
	data["gss"] = getgss(jsonDict, var_data, config) 

	# 体细胞I/II类Snvindel变异IGV图-20221020
	# 适用福建附一CP40模板
	data["igv_I_II"] = [var["igvplot"] for var in data["var_somatic"]["level_I"] + data["var_somatic"]["level_II"] if var["bio_category"] == "Snvindel" and "igvplot" in var.keys() and var["igvplot"]]

	#胚系4/5类Snvindel变异IGV图-2022.11.18
	# 适用福建附一HRR模板
	data["igv_4_5"] = [var["igvplot"] for var in data["var_germline"]["level_5"] + data["var_germline"]["level_4"] if var["bio_category"] == "Snvindel" and "igvplot" in var.keys() and var["igvplot"]]

	# 体细胞肿瘤发生发展相关变异snvindel IGV图-2023.06.27
	# 适用武汉协和HRR
	data["igv_onconodrug"] = [var["igvplot"] for var in data["var_somatic"]["level_onco_nodrug"] if var["bio_category"] == "Snvindel" and "igvplot" in var.keys() and var["igvplot"]]
	

	#### BRCA #####################################################################
	data_BRCA = {}
	data_BRCA["snv_s"] = {}
	data_BRCA["snv_m"] = {}
	#data_BRCA["mlpa"], data_BRCA["mlpa_image"], data_BRCA["mlpa_image_del"] = process_mlpa(jsonDict, config)
	data_BRCA["mlpa"] = mlpa
	data_BRCA["mlpa_image"] = mlpa
	data_BRCA["mlpa_image_del"] = mlpa_image_del
	data_BRCA["igvplot"] = {}
	data_BRCA["refer"] = {}
	
	for gene in ["BRCA1", "BRCA2"]:
		# snvindel
		for i in range(1, 6):
			# 单独组织/单独血液，等级判断：胚系：clinic_num_g; 体细胞：clinic_num_s 可能存在L5中存在胚系疑似致病或L4中存胚系致病，注意识别！！
			data_BRCA["snv_s"]["B"+gene[-1]+"_L"+str(i)] = [var for var in snvindel if var["gene_symbol"] == gene and ((var["var_origin"] == "germline" and var["clinic_num_g"] == i) or (var["var_origin"] != "germline" and var["clinic_num_s"] == i))]
			# 仅提取胚系结果，用于配对模板中需要拆分为两份的血液报告
			data_BRCA["snv_m"]["B"+gene[-1]+"_G_L"+str(i)] = [var for var in snvindel if var["gene_symbol"] == gene and var["var_origin"] == "germline" and var["clinic_num_g"] == i]

		# IGV
		data_BRCA["igvplot"]["all"] = [var["igvplot"] for var in snvindel if ((var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5]) or (var["var_origin"] != "germline" and var["clinic_num_s"] in [4, 5])) and var["gene_symbol"] in ["BRCA1","BRCA2"]]

		# 仅提取胚系结果，用于配对模板中需要拆分为两份的血液报告
		data_BRCA["igvplot"]["germline"] = [var["igvplot"] for var in snvindel if var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5] and var["gene_symbol"] in ["BRCA1","BRCA2"]]

		# 获取参考文献，不同模板展示内容不同
		# 目前josn中参考文献是汇总形式，无法直接使用，故脚本中需重新提取
		data_BRCA["refer"] = {}
		# 1. 变异描述PMID汇总
		data_BRCA["refer"]["var"] = []
		data_BRCA["refer"]["var_G"] = []
		# 2. 基因介绍PMID汇总
		data_BRCA["refer"]["gene"] = []
		data_BRCA["refer"]["gene_G"] = []
		# 3. 证据描述PMID汇总
		data_BRCA["refer"]["evi"] = []
		data_BRCA["refer"]["evi_G"] = []
		for var in snvindel:
			## a. BRCA 单样本项目
			# 新增条件，限制为BRCA基因-2023.05.06
			if ((var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5]) or (var["var_origin"] != "germline" and var["clinic_num_s"] in [4, 5])) and \
			var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				data_BRCA["refer"]["var"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
				data_BRCA["refer"]["gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
				data_BRCA["refer"]["evi"].extend(var["evi_sum"]["refer_evi"])
			## b. BRCA 配对样本项目
			if var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5] and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				data_BRCA["refer"]["var_G"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
				data_BRCA["refer"]["gene_G"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
				data_BRCA["refer"]["evi_G"].extend(var["evi_sum"]["refer_evi"])
			
	## MLPA文献漏加了，补充一下-刘炜芬-2023.11.07
	for var in mlpa["B1_Loss"] + mlpa["B2_Loss"]:
		data_BRCA["refer"]["var"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		data_BRCA["refer"]["evi"].extend(var["evi_sum"]["refer_evi"])
		data_BRCA["refer"]["var_G"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene_G"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		data_BRCA["refer"]["evi_G"].extend(var["evi_sum"]["refer_evi"])
	for var in mlpa["B1_Gain"] + mlpa["B2_Gain"]:
		data_BRCA["refer"]["var"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		data_BRCA["refer"]["var_G"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		data_BRCA["refer"]["gene_G"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
	# MLPA文献添加完成-2023.11.07

	# 存放总结性特殊需求
	data_BRCA["special"] = {}	
	# 特殊-广东人民小结
	data_BRCA["special"]["summary_GDRM_G"], data_BRCA["special"]["summary_GDRM_S"] = var_summary_GDRM(data_BRCA)
	# 特殊-贵州肿瘤非4,5类突变
	data_BRCA["special"]["GZZL_level13_G"], data_BRCA["special"]["GZZL_level13_S"] = var_GZZL(data_BRCA)
	# 特殊-上海仁济，结果结论拼接
	data_BRCA["special"]["summary_SHRJ"] = var_summary_SHRJ(data_BRCA)
	# 特殊-福建肿瘤（BRCA1基因致病性变异）
	data_BRCA["special"]["summary_FJZL"] = var_summary_FJZL(data_BRCA)
	# 特殊-山东省立/山东齐鲁非4，5类突变
	data_BRCA["special"]["SD_level1_2_G"], data_BRCA["special"]["SD_level1_2_S"] = var_SD(data_BRCA)
	# 特殊-武汉协和总结（1,2,3类和4,5类）
	data_BRCA["special"]["summary_WHXH_level45"], data_BRCA["special"]["summary_WHXH_level123"], data_BRCA["special"]["summary_WHXH_level123_G"] = var_summary_WHXH(data_BRCA)
	# 特殊-山东肿瘤非4,5类
	data_BRCA["special"]["summary_SDZL_level123_G"], data_BRCA["special"]["summary_SDZL_level123_S"] = var_level123_SDZL(data_BRCA)
	# 特殊-判断体细胞变异最高等级
	data_BRCA["special"]["judge_intumor"] = brca_judge_inApprovalTumor(data_BRCA)
	# 特殊-重庆附一-进院全血总结
	data_BRCA["special"]["sum_g_CQFY"] = var_summary_CQFY(data_BRCA)
	#print (data["special"]["sum_g_CQFY"])
	# 特殊-福建附一-2022.11.11
	data_BRCA["special"]["summary_FJFY"] = BRCA_FJFY_summary(data_BRCA)
	

	### HRR for SHSY ######################################
	data_HRR_SHSY = {}
	# 更新排序规则，按I/II（同等级时按治疗方案最高等级进行排序）、基因、类型、频率（频率在各部分处理的时候已经排序了）
	# 同等级时按治疗方案最高等级进行排序 A > B > C > D > N
	shsy_top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
	shsy_var_data = sorted(snvindel, key=lambda i : (i["clinic_num_s"], shsy_top_level_rule.get(i["top_level"]), i["gene_symbol"]))
	brca_list = ["BRCA1", "BRCA2"]
	hrr_list = ["ATM","ATR","BARD1","BRIP1","CDK12","CHEK1","CHEK2","FANCA","FANCL","MRE11A","NBN","PALB2","PPP2R2A","RAD51B","RAD51C","RAD51D","RAD54L"]
	other_list = ["AR","CDH1","ESR1","HDAC2","HOXB13","PTEN","STK11","TP53","BRAF","ERBB2","KRAS","NRAS","PIK3CA"]

	# 体细胞结果整理
	data_HRR_SHSY["var_s_brca"] = s_var_rule_genelist(shsy_var_data, brca_list)
	data_HRR_SHSY["var_s_hrr"] = s_var_rule_genelist(shsy_var_data, hrr_list)
	data_HRR_SHSY["var_s_other"] = s_var_rule_genelist(shsy_var_data, other_list)

	# 胚系结果整理
	data_HRR_SHSY["var_g_brca"] = g_var_rule_genelist(shsy_var_data, brca_list)
	data_HRR_SHSY["var_g_hrr"] = g_var_rule_genelist(shsy_var_data, hrr_list)
	data_HRR_SHSY["var_g_other"] = g_var_rule_genelist(shsy_var_data, other_list)
	

	
	return data, data_BRCA, data_HRR_SHSY