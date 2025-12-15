#-*- coding:gbk -*-
import os
import re
import xlrd

'''
Discription
	
	获取config文件夹中的报告模板配置文件，转化为便于其他程序调用的格式。 
	{
		临检 : {
			定制 : {
				*新增 (医院，产品1，项目) ：模板名1,
				(医院, 产品1) : 模板名1
				},
			通用-简版 : {
				产品1 : 模板名1
				},
			通用-完整 : {
				产品1 : 模板名1
				}
			},
		进院 : {
			定制 : {
				(医院, 科室, 产品1) : 模板名1,
				(医院, 产品1) : 模板名1
 				},
			通用 : {
				产品1 : 模板名1
				}
		}
	
'''

#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#requirenment_path = os.path.join(BASE_DIR, "config/report_requirenment.xlsx")

key_stran = {
	"产品名称" : "prod_name",
	"检测业务类型" : "business_type",
	"模板类型" : "report_type",
	"单位名称" : "company",
	"科室" : "hosp_depart",
	"模板名" : "report_name",
	"状态" : "status",
	"模板开发者" : "developer",
	"添加人" : "auditors",
	"添加时间" : "auditors_time",
	"备注" : "note",
	"更新记录" : "update",
	"临检" : "rummage",
	"进院" : "hospital",
	"定制" : "CustomEdition",
	"通用" : "Universal",
	"通用-简版" : "Universal_simple",
	"通用-完整" : "Universal_complete",
	"基础模板拼接" : "merge_template",
	"基础模板-动态" : "temp_dynamic",
	"基础模板-固定" : "temp_fixed",
	"药企项目名称" : "product_name"
	}

def getRequirenment(config):
	#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	requirenment_path = os.path.join(config, "report_requirenment.xlsx")
	xls = xlrd.open_workbook(requirenment_path)
	requirenment_sheet = xls.sheet_by_name("report_requirenment")
	key = requirenment_sheet.row_values(0)
	Data = []
	for num in range(1, requirenment_sheet.nrows):
		rows = requirenment_sheet.row_values(num)
		tmpdict = {}
		for i in range(len(key)):
			tmpdict[key_stran.get(key[i])] = key_stran[rows[i]] if rows[i] in key_stran.keys() else rows[i]
		Data.append(tmpdict)

	#配置信息字典结构	
	requir_dict = {
		"rummage" : {
			"CustomEdition" : {},
			"Universal_simple" : {},
			"Universal_complete" : {}
			},
		"hospital" : {
			"CustomEdition" : {},
			"Universal" : {}
			}
		}

	# 科室目前仅在进院样本中需要做区分，临检和药企就不考虑了
	for i in Data:
		if str(int(i["status"])) == "0":
			if i["company"]:
				if i["hosp_depart"]:
					requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["hosp_depart"], i["prod_name"]), "")
					requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["hosp_depart"], i["prod_name"])] = i["report_name"]
				else:
					# 新增药企配置-2023.11.17
					if i["product_name"]:
						requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"], i["product_name"]), "")
						requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"], i["product_name"])] = i["report_name"]
					else:
					# 新增结束-2023.11.17
						#print (i)
						requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"]), "")
						requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"])] = i["report_name"]
			else:	
				requir_dict[i["business_type"]][i["report_type"]].setdefault(i["prod_name"], "")
				requir_dict[i["business_type"]][i["report_type"]][i["prod_name"]] = i["report_name"]

	# 是否需要拼接基础模板对照表
	report_merge_dict = {}
	for i in Data:
		if str(int(i["status"])) == "0":
			if i["report_name"] not in report_merge_dict.keys():
				report_merge_dict.setdefault(i["report_name"],{})
			report_merge_dict[i["report_name"]] = {
				"merge_template" : str(int(i["merge_template"])),
				"temp_dynamic" : re.split(";", i["temp_dynamic"].replace("\n","").replace(" ","")) if i["temp_dynamic"] else "",
				"temp_fixed" : re.split(";", i["temp_fixed"].replace("\n","").replace(" ","")) if i["temp_fixed"] else ""
			}

	return requir_dict, report_merge_dict
