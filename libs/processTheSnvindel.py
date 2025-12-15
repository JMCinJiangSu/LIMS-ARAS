#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran, typeStran, typeStran_from_inter
from libs.getEvi import varRegimen
import copy
from libs.specialRequest import varInfo_XAJDY, varInter_FJZL, varInfo_SYX
from libs.hgvsp_To_hgvsp_abbr import splitAA
from libs.specialRequest import varInfo_FDZS, varInfo_FJZL
from libs.rule import S_function
from libs.rule import decimal_float, decimal_percen
from libs.getConfig import get_gene_class
import datetime
from decimal import Decimal

'''
Discription
	
	处理snvindel格式。 

'''

def process_snvindel(jsonDict, config):
	snvindel = copy.deepcopy(jsonDict["snvindel"])
	#print  ("######### 3.0.1 step0", datetime.datetime.now())
	function = get_gene_class(config)["function"]
	clinical = get_gene_class(config)["clinical"]
	clinicalNumStran_dict = clinicalNumStran(config)
	functionNumStran_dict = functionNumStran(config)
	for var in snvindel:
		# 致病性/致癌性都有，则胚系匹配致病性，体细胞匹配致癌性；
		# 致病性/致癌性仅有一个，胚系/体细胞仅匹配有数据的那个；
		# 致病性/致癌性均无，胚系/体细胞均默认为3。
		if var["clinical_significance"] != "-" or var["function_classification"] != "-":
			var["clinic_num_g"] = clinicalNumStran_dict.get(var["clinical_significance"], 3) if var["clinical_significance"] != "-" else \
								  clinicalNumStran_dict.get(var["function_classification"], 3) 
			var["clinic_num_s"] = functionNumStran_dict.get(var["function_classification"], 3) if var["function_classification"] != "-" else \
								  functionNumStran_dict.get(var["clinical_significance"], 3)
		else:
			var["clinic_num_g"] = 3
			var["clinic_num_s"] = 3
		# AD3101 新增一个字段，只有致病性解读的返回致病性等级，其他情况按致癌性解读等级输出,2025年3月6日
		#var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 3
		#var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 3
		# 更新ad3101 嵇梦晨 2025年6月18日，修改默认值
		var['clinic_ad3101_g'] = var['clinic_num_g'] if var["clinical_significance"] != "-" and var['function_classification'] == '-' else 0
		var['clinic_ad3101_s'] = var['clinic_num_s'] if var['function_classification'] != '-' else 0

		# freq 兼容百分数格式，均转化为小数格式-适配v4-刘炜芬，2023.11.09
		# freq 兼容字符串格式，均转化为小数格式-嵇梦晨，2024.03.20
		freq_key = ["freq", "freq_sc", "freq_rc", "freq_ctrl", "freq_case"]
		for key in freq_key:
			if key in var.keys():
				var[key] = float(var[key].replace("%", ""))/100 if re.search("%", str(var[key])) else float(var[key]) if isinstance(var[key], str) and var[key] else var[key]
		# 兼容结束-2023.11.09

		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		var["freq_str"] = decimal_percen(var["freq"]) if var["freq"] else ""
		var["freq_sc_str"] = decimal_percen(var["freq_sc"]) if var["freq_sc"] else ""
		# 兼容freq_rc未返回的情况
		var["freq_rc"] = var["freq_rc"] if "freq_rc" in var.keys() and var["freq_rc"] else var["freq_ctrl"] if "freq_ctrl" in var.keys() and var["freq_ctrl"] else ""
		var["freq_rc_str"] = decimal_percen(var["freq_rc"]) if "freq_rc" in var.keys() and var["freq_rc"] else decimal_percen(var["freq_ctrl"]) if var["freq_ctrl"] else ""
		var["freq_ctrl_str"] = decimal_percen(var["freq_ctrl"]) if var["freq_ctrl"] else ""
		var["freq_ss_str"] = decimal_percen(var["freq_ss"]) if var["freq_ss"] else ""
		var["freq_case_str"] = decimal_percen(var["freq_case"]) if var["freq_case"] else ""
		
		# 特殊处理
		# 1. hgvs_p加括号-默认格式， 新增兼容（防止新增变异未填hgvs_p）
		var["hgvs_p"] = var["hgvs_p"] if var["hgvs_p"] else "p.?"
		if not re.search("=", var["hgvs_p"]) and var["hgvs_p"] != "p.?" and not re.search("\(", var["hgvs_p"]):
			var["hgvs_p"] = var["hgvs_p"].replace("p.", "p.(")+")"
		# 2. hgvs_p_abbr加括号-默认格式
		var["hgvs_p_abbr"] = var["hgvs_p_abbr"] if var["hgvs_p_abbr"] else ""
		if var["hgvs_p_abbr"] and not re.search("=", var["hgvs_p_abbr"]) and var["hgvs_p_abbr"] != "p.?" and not re.search("\(", var["hgvs_p_abbr"]):
			var["hgvs_p_abbr"] = var["hgvs_p_abbr"].replace("p.", "p.(")+")"
		# 3. hgvs_p除了同义突变，其余不加括号-适用浙肿
		var["hgvs_p_ZJZL"] = var["hgvs_p"].replace("(", "").replace(")", "") if not re.search("=", var["hgvs_p"]) else var["hgvs_p"]
		# 4. hgvs_p全部不加括号-适用于CP40福建附一
		var["hgvs_p_FJFY"] = var["hgvs_p"].replace("(", "").replace(")", "")
		# 5. 不带小版本号的转录本		
		var["transcript_primary_simple"] = re.split("\.", str(var["transcript_primary"]))[0] if var["transcript_primary"] else ""
		# 6. 变异类型转化为中文，新增变异时变异类型可能会为空，会导致程序报错，这边做下兼容
		#var["type_"] = var["type"] if var["type"] else ""
		var["type_cn"], var["type"] = typeStran_from_inter(var)
		var["type_cn"] = var["type_cn"] if var["type_cn"] else ""
		# 临时加一个AltDepth用于郑大一，等极元返回后再用系统自带的，2023.07.03
		if var["depth"] and var["freq"]:
			var["altdepth"] = str(Decimal(float(var["depth"]) * float(var["freq"])).quantize(Decimal("0"), rounding="ROUND_HALF_UP"))
		else:
			var["altdepth"] = "深度或频率为空，请人工填写！"

		# 变异特殊需求
		# 1. 福建肿瘤变异描述
		var["varInter_FJZL"] = varInter_FJZL(var["hgvs_p"], var["type"], var["gene_region"], var["hgvs_c"])
		# 2. 西安交大一：第X号外显子/内含子XXX突变
		var["varInfo_XAJDY"] = varInfo_XAJDY(var["gene_region"], var["type_cn"]) if var["type_cn"] != "--" else varInfo_XAJDY(var["gene_region"], "突变")
		# 3. 孙逸仙：XX号外显子/内含子XX hgvs_p/hgvs_c XX突变
		var["varInfo_SYX"] = varInfo_SYX(var["gene_region"], var["type_cn"], var["type"], var["hgvs_c"], var["hgvs_p_ZJZL"])
		# 4. 若三字母氨基酸为空，则使用报告脚本转化的结果-2022.09.07
		var["hgvs_p_abbr"] = splitAA(var["hgvs_p"]) if not var["hgvs_p_abbr"] else var["hgvs_p_abbr"]
		# 5. 复旦中山变异描述特殊需求
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], var["gene_region"], var["variant_desc_cn"], config)
		# 6. 浙江肿瘤：变异描述统一下格式，需要加“可能形成功能损伤或失活的蛋白”，可以加在报告里，但为了避免知识库中该条信息格式不统一（如有空格）造成报告展示混乱，这边预处理一下
		var["variant_desc_cn_ZJZL"] = var["variant_desc_cn"].strip()[0:-1] if var["variant_desc_cn"] and var["variant_desc_cn"].strip()[-1] == "。" else var["variant_desc_cn"]
		# 7. 福建肿瘤：返回变异频率相关信息（来源配置表）和治疗方案汇总
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		# 8. 深圳医院：XX号外显子/内含子XX hgvs_p/hgvs_c XX突变 同孙逸仙，但是hgvs_p需要保留括号-2023.08.03
		var["varInfo_SZYY"] = varInfo_SYX(var["gene_region"], var["type_cn"], var["type"], var["hgvs_c"], var["hgvs_p"])
		# 9. XW0413 EGFR分子标签 Ex19del Ex20ins  嵇梦晨
		# EGFR_tag不需要括号及 "p." 2024年6月17日 嵇梦晨
		#var['EGFR_tag'] = 'Ex19 del' if 'EGFR Exon19 del' in var['var_category_names'] else 'Ex20 ins' if 'EGFR Exon20 ins' in var['var_category_names'] \
		#	else var['hgvs_p'].replace('p.', '') if not re.search('=', var['hgvs_p']) and var['hgvs_p'] != 'p.?' else var['hgvs_p']
		# json返回改成null了，加个判断 2024年10月29日 jmc
		if var['var_category_names']:
			var['EGFR_tag'] = 'Ex19 del' if 'EGFR Exon19 del' in var['var_category_names'] else 'Ex20 ins' if 'EGFR Exon20 ins' in var['var_category_names'] else var['hgvs_p_ZJZL'].replace('p.', '') if var['hgvs_p_ZJZL'] != 'p.?' else var['hgvs_p_ZJZL']
		#var['EGFR_tag'] = 'Ex19 del' if 'EGFR Exon19 del' in var['var_category_names'] else 'Ex20 ins' if 'EGFR Exon20 ins' in var['var_category_names'] \
		#	else var['hgvs_p_ZJZL'].replace('p.', '') if var['hgvs_p_ZJZL'] != 'p.?' else var['hgvs_p_ZJZL']
		# 10. XW5101 判断NPM1变异分型
		if jsonDict["sample_info"]["product_name"] == "XW5101":
			if var["gene_symbol"] == "NPM1" and re.search("863", var["hgvs_c"]):
				var["NPM1_type"] = "Type A" if var["alt"] == "TCTG" else "Type B" if var["alt"] == "CATG" else "Type D" if var["alt"] == "CCTG" else ''

	#print  ("######### 3.0.2 step2", datetime.datetime.now())
	# 按频率排序下
	snvindel = sorted(snvindel, key=lambda i:i["freq"], reverse=True)
	#print  ("######### 3.0.3 step3", datetime.datetime.now())
	return snvindel