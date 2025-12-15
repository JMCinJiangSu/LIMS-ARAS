#-*- coding:gbk -*-
import re
import copy
from functools import reduce
from libs.rule import decimal_float, decimal_percen

'''
Discription
	
	获取json文件中的qc，转化为报告模板方便填充的格式。 

'''

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

# 用于判断测序仪，Miseq/NextSeq、华大-2023.07.03
# 规则：1）五位数的为Miseq；2）S开头的为华大测序仪；3）A开头的为CN500
# 规则更新：外部测序的编号可能会填成“SJMXXX”，属于CN500，程序会判定为华大，这边加个判定-2023.11.01
# 开头为S+一个数字的判定为华大，其他返回"huada/cn500"，标个颜色
def judge_Sequencer(flowcell_lane):
	flowcell = re.split("_", flowcell_lane)[0] if flowcell_lane else ""
	if len(flowcell) == 5:
		return "miseq"
	elif re.match("S", flowcell):
	#	return "huada"
		if re.match("S[0-9]", flowcell):
			return "huada"
		else:
			return ""
	elif re.match("A", flowcell):
		return "cn500"
	else:
		return ""

def QCStran_dict(qcDict):
	QC_result = {}
	for k, v in qcDict.items():
		if k != "qc_type" and v:
			QC_result[k+"_num"] = float(v[:-1]) if re.search("cleandata_size|dv200|tumor_content|tumor_content_macrodissection_performed", k) else float(v) if is_number(v) else v
			#QC_result[k] = "{:.2%}".format(float(v)) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration", k) else v if re.search("cleandata_size", k) else "{:.2f}".format(float(v)) if is_number(v) else v
			QC_result[k] = decimal_percen(v) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration|registration_certificate", k) else \
						   v if re.search("cleandata_size", k) else \
						   decimal_float(v) if is_number(v) else \
						   v

		else:
			QC_result[k+"_num"] = 0
			QC_result[k] = ""
		# 新增测序仪字段-2023.07.30
		if k == "flowcell_lane":
			QC_result["Sequencer"] = judge_Sequencer(v)
		# 新增结束-2023.07.03
	return QC_result

def QCStran_list(qclist):
	QC_result = {}
	for i in qclist:
		QC_result.setdefault(i["qc_type"], {})
		for k, v in i.items():
			if k != "qc_type" and v:
				QC_result[i["qc_type"]][k+"_num"] = float(v[:-1]) if re.search("cleandata_size|dv200|tumor_content|tumor_content_macrodissection_performed", k) else float(v) if is_number(v) else v
				#QC_result[i["qc_type"]][k] = "{:.2%}".format(float(v)) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration", k) else v if re.search("cleandata_size", k) else "{:.2f}".format(float(v)) if is_number(v) else v
				QC_result[i["qc_type"]][k] = decimal_percen(v) if re.search("q30|q20|ratio|uni20", k) and not re.search("concentration", k) else \
											 v if re.search("cleandata_size", k) else \
											 decimal_float(v) if is_number(v) else \
											 v
			else:
				QC_result[i["qc_type"]][k+"_num"] = 0
				QC_result[i["qc_type"]][k] = ""	
			# 新增测序仪字段-2023.07.30
			if k == "flowcell_lane":
				QC_result["Sequencer"] = judge_Sequencer(v)
			# 新增结束-2023.07.03
	return QC_result
	
# 主函数
def getJsonQC(jsonDict):
	qc = copy.deepcopy(jsonDict["qc"])
	# 用于XW0248 判断CNV质控是否合格
	#qc['dna_data_qc'][0]['cnv_noise'] = jsonDict['cnv'][0]['cnv_noise'] if jsonDict['cnv'] else ''
	
	qc_items = [k for k in qc]
	data = {}
	for item in qc_items:
		# qc_gradient单独处理
		if item == "qc_gradient" and qc["qc_gradient"]:
			data["qc_gradient"] = {}
			for i in qc["qc_gradient"]:
				data["qc_gradient"][i["qc_source"]+"_"+i["gradient_num"]] = "{:.2%}".format(float(i["gradient_ratio"])) 
		else:
			# 若质控内容只有一个且为字典格式时，使用该转换方式
			if type(qc[item]).__name__=="dict":
				data[item] = QCStran_dict(qc[item])
			# 若质控内容为列表
			else:
				# 列表长度为1时，转化为字段，使用字典的转换方式
				if len(qc[item])==1:
					data[item] = QCStran_dict(qc[item][0])
				# 列表中含有多个元素时，使用列表的转换方式
				elif len(qc[item]) >= 2:
					# 先进行去重，防止出现重复的质控数据导致报错
					qc_data = reduce(lambda x, y:x if y in x else x + [y], [[],]+qc[item])
					# 去重后还报错的原因可能是数据审核时没有把无关数据F掉！！
					if len(qc_data) ==1:
						data[item] = QCStran_dict(qc_data[0])
					else:
						data[item] = QCStran_list(qc[item])
	
	# 新增湿实验质控-2022.08.10
	qc_lib = copy.deepcopy(jsonDict["lib_quality_control"]) if "lib_quality_control" in jsonDict.keys() and jsonDict["lib_quality_control"] else {}
	qc_lib_items = [k for k in qc_lib]
	lib_data = {}
	for item in qc_lib_items:
		if type(qc_lib[item]).__name__=="dict":
			lib_data[item] = QCStran_dict(qc_lib[item][0])
		else:
			if len(qc_lib[item]) == 1:
				lib_data[item] = QCStran_dict(qc_lib[item][0])
			elif len(qc_lib[item]) >= 2:
				qc_data = reduce(lambda x, y:x if y in x else x + [y], [[],]+qc_lib[item])
				if len(qc_data) == 1:
					lib_data[item] = QCStran_dict(qc_data[0])
				else:
					lib_data[item] = QCStran_list(qc_lib[item])

	return data, lib_data