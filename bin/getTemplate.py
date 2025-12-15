#-*- coding:gbk -*-
import re
from libs.getRequirExcel import getRequirenment
from libs.getProdName_alias import alias_name

'''
Discription
	
	匹配报告模板配置表，获取模板名称，用于报告填充。 

'''

def MatchReport(jsonDict, config):
	requir_dict, report_merge_dict = getRequirenment(config)
	i = jsonDict.get("sample_info")
	# 获取产品别名，转化为主要产品名，便于匹配模板
	alias_name_dict = alias_name(config)
	i["prod_names"] = alias_name_dict[i["prod_names"]] if i["prod_names"] in alias_name_dict.keys() else i["prod_names"]
	#匹配报告模板
	# 1. 临检
	 # 新增药企配置-匹配放前面-2023.11.17
	 # 1.0 先匹配"临检:定制:医院-产品-项目"
	 # 1.1 先匹配"临检:定制:医院-产品"
	 # 1.2 无定制的话汇总匹配"临检:通用-简版-产品"， 否则匹配"临检:通用-完整-产品"
	# 2. 进院
	 # 2.1 先匹配定制"进院:定制:医院-科室-产品"，无的话匹配"进院:定制:医院-产品"
	 # 2.2 无定制的话匹配"进院:通用:产品"
	report_name = ""
	if i["report_module_type"] == "rummage":
		# 20220627-这边做个兼容：临检样本需要匹配带后缀送检单位名称，而手动上传的订单缺失该项内容，若需要用手动上传订单出临检报告，会报错。
		# 这边兼容成缺少带后缀送检单位名称时，直接出临检完整版报告
		i["origin_company"] = i["origin_company"] if i["origin_company"] else i["company"]
		# 兼容完成-20220627
		# 新增药企配置-2023.11.17，单位名称先用company吧，现在还是手动上传的订单，不清楚对接LIMS后会是什么样
		if (i["company"], i["prod_names"], i["product_name"]) in requir_dict["rummage"]["CustomEdition"].keys():
			report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["product_name"]))
		else:
		# 新增药企配置代码更新结束-2023.11.17
			if (i["origin_company"], i["prod_names"]) in requir_dict["rummage"]["CustomEdition"].keys():
				report_name = requir_dict["rummage"]["CustomEdition"].get((i["origin_company"], i["prod_names"]))
			else:
				report_name = requir_dict["rummage"]["Universal_simple"].get(i["prod_names"], "") if re.search("汇总", i["origin_company"]) else requir_dict["rummage"]["Universal_complete"].get(i["prod_names"], "")
	elif i["report_module_type"] == "hospital":
		if (i["company"], i["hosp_depart"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
		elif (i["company"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
		else:
			report_name = requir_dict["hospital"]["Universal"].get(i["prod_names"], "")
	
	merge_template = report_merge_dict.get(report_name)

	return report_name, merge_template