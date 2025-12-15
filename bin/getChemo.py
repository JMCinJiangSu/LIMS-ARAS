#-*- coding:gbk -*-
import copy
from functools import reduce

'''
Discription
	
	改脚本用来获取化疗信息，并整理为便于报告填充的格式，目前有的格式如下：
	1. 完整描述，直接使用json返回的信息填充就行，全部展示，癌种相关的加粗显示
	2. 以位点为key，简单描述，CP40仅肠癌的时候展示
	3. 以药物为key，简单描述，分癌种（暂未添加）
	4. 完整描述，需要区分癌种，使用json返回信息中is_same_tumor=1的内容，若无则放全部

'''

def getchemo(jsonDict):
	chemo_data = copy.deepcopy(jsonDict["PGx"])
	chemo_result = {}
	# 1. complete适用于MP临检通用
	chemo_result["complete"] = sorted(chemo_data, key=lambda i:(i["gene_symbol"], i["dbsnp"]))
	same_tumor_data = [i for i in chemo_result["complete"] if i["is_same_tumor"] == 1]
	# 2. dbsnp_simple_nosplitdrug适用于CP40-未区分癌种
	chemo_result["dbsnp_simple_nosplitdrug"] = process_dbsnp_nosplitdrug(chemo_data)
	# 2.2 适用于北大人民病理科HRR-区分癌种，若对应癌种无化疗信息，则放全部-2022.09.19新增
	chemo_result["dbsnp_simple_splittumor"] = process_dbsnp_nosplitdrug(same_tumor_data) if same_tumor_data else process_dbsnp_nosplitdrug(chemo_result["complete"])
	# 3. HRR,完整版描述，但需要区分癌种，若对应癌种无化疗信息，则放全部chemo_result["complete"] 【is_same_tumor==1 为同癌种，==0为非同癌种】
	chemo_result["complete_split_tumor"] = reduce_chemo(same_tumor_data) if same_tumor_data else reduce_chemo(chemo_result["complete"])
	# 4. 116化疗，去重过的
	chemo_result["reduce_116"] = process_116(chemo_data)
	# 5. 以药物为主，且分癌种（适用HRR北大人民检验科-旧等模板）
	chemo_result["drug_splittumor"] = process_drug_splittumor(chemo_data)
	return chemo_result

def process_dbsnp_nosplitdrug(chemo_data):
	# 适用于CP40
	effict_drug_dict = {}
	for i in chemo_data:
		key = (i["gene_symbol"], i["dbsnp"], i["genotype"], i["evi_level"])
		effict_drug_dict.setdefault(key, [])
		effict_drug_dict[key].append((i["drug_name_cn"], i["impact_type_cn"]))
	dbsnp_list = set([(i["gene_symbol"], i["dbsnp"], i["genotype"], i["evi_level"]) for i in chemo_data])

	result = []
	for info in dbsnp_list:
		tmpdict = {}
		drug_inter = []
		if info in effict_drug_dict.keys():
			for i in effict_drug_dict[info]:
				tmpdict.setdefault(i[1], [])
				tmpdict[i[1]].append(i[0])
		for k, v in tmpdict.items():
			drug_inter.append("、".join(v)+str(k))
		result.append({
			"gene_symbol" : info[0],
			"dbsnp" : info[1],
			"genotype" : info[2],
			"evi_level" : info[3],
			"inter" : drug_inter
		})

	return sorted(result, key = lambda i:(i["gene_symbol"], i["dbsnp"]))

def process_116(chemo_data):
	# 116化疗部分要去重，让改知识库不改，非要代码去重，那我只能再写个函数单独处理。
	# 字典还要重构，原始的里面还有药物名称，乱去重的话可能会造成其他报告不准确的风险。
	chemo_116 = [{"gene_symbol" : i["gene_symbol"], "dbsnp" : i["dbsnp"], "genotype" : i["genotype"], "clin_anno_cn" : i["clin_anno_cn"], "evi_level" : i["evi_level"]} for i in chemo_data]
	chemo_116 = reduce(lambda x, y:x if y in x else x + [y], [[],]+chemo_116)
	return sorted(chemo_116, key=lambda i:(i["gene_symbol"], i["dbsnp"]))

def process_drug_splittumor(chemo_data):
	# 以药物为主进行数据处理，适用HRR北大人民检验科，若同癌种的无化疗结果，则放泛癌种
	split_tumor_data = [i for i in chemo_data if i["is_same_tumor"] == 1]
	split_tumor_data = split_tumor_data if split_tumor_data else chemo_data
	result = []
	result_tmp = {}
	for i in split_tumor_data:
		result_tmp.setdefault(i["drug_name_cn"], [])
		result_tmp[i["drug_name_cn"]].append(i)
	
	for k,v in result_tmp.items():
		v_sort = sorted(v, key=lambda i:(i["evi_level"], i["gene_symbol"]))
		result.append({
			"drug_name_cn" : k,
			"info" : v_sort
		})
	result = sorted(result, key=lambda i:i["drug_name_cn"])
	return result

# 去重，适用HRR
def reduce_chemo(chemo_data):
	result = []
	chemo_list = [{"gene_symbol" : i["gene_symbol"], "dbsnp" : i["dbsnp"], "genotype" : i["genotype"], "clin_anno_cn" : i["clin_anno_cn"], "evi_level" : i["evi_level"]} for i in chemo_data]
	chemo_list = reduce(lambda x, y:x if y in x else x + [y], [[],]+chemo_list)
	return sorted(chemo_list, key=lambda i:(i["gene_symbol"], i["dbsnp"]))