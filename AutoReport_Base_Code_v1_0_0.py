#-*- coding:gbk -*-
import re
import os
import json
import datetime
import argparse
import time
import copy
import lxml
from docxtpl import DocxTemplate, InlineImage, R
from docx.shared import Mm
from bin import getSampleInfo
from bin import getQC
from bin import getTemplate
from bin import getApprovalDrug
#from bin import getClinicTrial
from bin import getPDL1
from bin import getMSI
from bin import getTMB
from bin import getVar
from bin import getRNAexp
#from bin import getChemo
from bin import getRefence
from bin import getGEP
from bin import getHRD
from bin import getTME
from bin import getApprovalRegimen
from bin import baseTemplate
from libs.getProdName_alias import alias_name
from bin import getImage
import customize_filters
from libs import processMRD
from bin import getRelatedVar

### 防止药企环境的报错，先把不用的模块注释掉-2023.11.17 ###

def get_data(json_name, outfile, config, report_template, outjson, image):
	print ("开始生成填充数据！", datetime.datetime.now())
	data = {}
	# 解析json文件
	with open(outfile+"/"+json_name+".json", "r", encoding='utf-8') as file_json:
		jsonDict = json.load(file_json)
		# 系统将rummage改为clinical了，为了报告脚本代码改动最小化，在这边做转化
		jsonDict["sample_info"]["report_module_type"] = "rummage" if jsonDict["sample_info"]["report_module_type"] == "clinical" else jsonDict["sample_info"]["report_module_type"]
		# 将产品别名进行转换
		alias_name_dict = alias_name(config)
		jsonDict["sample_info"]["prod_names"] = alias_name_dict[jsonDict["sample_info"]["prod_names"]] if jsonDict["sample_info"]["prod_names"] in alias_name_dict.keys() else jsonDict["sample_info"]["prod_names"]
	# 模板选择
	report_name, merge_template = getTemplate.MatchReport(jsonDict, config)
	print (report_name)
	#print (report_name)
	#print (merge_template)
	data["sample"] = getSampleInfo.getSample(jsonDict)
	data["sample"]["report_name"] = report_name
	data["qc"], data["lib_quality_control"] = getQC.getJsonQC(jsonDict)
	data["drug"] = getApprovalDrug.getDrug(jsonDict)
	data["therapeutic_regimen"] = getApprovalRegimen.getRegimen(jsonDict)
	#data["clinic_trial"] = getClinicTrial.getClinic(jsonDict)
	data["pdl1"] = getPDL1.getPDL1(jsonDict)
	data["msi"] = getMSI.getMSI(jsonDict, config)
	data["tmb"] = getTMB.getTMB(jsonDict, config)
	data["var"], data["var_brca"], data["var_hrr_shsy"] = getVar.getVar(jsonDict, config, report_name)
	data["rna_exp"] = getRNAexp.getRNA_exp(jsonDict)
	#data["chemo"] = getChemo.getchemo(jsonDict)
	data["gep"] = getGEP.getGEP(jsonDict)
	data["hrd"] = getHRD.getHRD(jsonDict, data["var"]["ec_type"]["BRCA1_level12"]+data["var"]["ec_type"]["BRCA2_level12"], config)
	data["tme"], data["tme_score"] = getTME.getTME(jsonDict)
	data["hd"] = customize_filters.xw1402_hd(jsonDict["hd"]) if jsonDict['sample_info']['product_name'] in ['XW1402', 'XW1405', 'XW1404'] else jsonDict['hd']
	# 处理MRD数据
	data["mrd"], data['mrd_last'] = processMRD.process_MRD(jsonDict)
	data["mrd_related"] = getRelatedVar.get_related_var(json_name, outfile, config, report_name)
	# 用于浙肿BRCA判定是否有检出CNV信号
	data["judge_CNV"] = "T" if jsonDict["cnv"] else ""
	# 参考文献
	data["refer"] = {}
	data["refer"]["fixed"] = getRefence.getfixed_refer(report_name, data["sample"]["tumor_list"], data["var"]["knb"], data["msi"], config)
	data["refer"]["dynamic"] = getRefence.getdynamic_refer(jsonDict, data["var"], data["msi"], data["hrd"], data["var_brca"])

	# 参数为T时，填充用数据转化为json输出，便于开发
	if outjson == "T":
		dataJson = json.dumps(data, ensure_ascii = False)
		with open(outfile+"/"+json_name+"_to_word.json", "w", encoding = "utf-8") as outFile:
			outFile.write(dataJson)

	# 加一个分页符-2023.05.17
	data["page_break"] = R("\f")
	# 模板填充
	if report_name:
		path = os.path.join(report_template, "template_main", report_name)
		tpl = DocxTemplate(path)
		# 图片填充
		data["image"] = getImage.render_image(tpl, data, jsonDict, report_name, image, config)
		print ("填充数据生成完毕！", datetime.datetime.now())
		print ("开始填充报告：", datetime.datetime.now())
		# 拼接模板-20221216
		data["subdoc"] = baseTemplate.BaseReport(data, tpl, merge_template, report_template, json_name, outfile)
		tpl.render(data)
		tpl.save(outfile+"/"+json_name+".docx")
		print ("报告填充完成！", datetime.datetime.now())
		

		# 加个自动更新域功能-仅供需要更新目录的模板使用，待开发-20220922
		#judge_update = "yes"
		#if judge_update:
			#namespace = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
			#element_updatefields = lxml.etree.SubElement(tpl.settings.element, namespace+"updateFields")
			#element_updatefields.set(namespace+"val", "true")
			#element_updatefields.set(namespace+"val","true")
			#tpl.save(outfile+"/"+json_name+".docx")

	else:
		print ("未匹配到报告模板")

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--json_name", dest = "json_name", required = True)
	parser.add_argument("-o", "--outfile", dest = "outfile", required = True)
	parser.add_argument("-c", "--config", dest = "config", required = True)
	parser.add_argument("-r", "--report_template", dest = "report_template", required = True)
	parser.add_argument("-j", "--outjson", dest = "outjson", required = True)
	parser.add_argument("-i", "--image", dest = "image", required = True)
	arg = parser.parse_args()

	return arg

if __name__ == '__main__':
	args = parse_args()
	get_data(json_name=args.json_name, outfile=args.outfile, config=args.config, report_template=args.report_template, outjson=args.outjson, image = args.image)