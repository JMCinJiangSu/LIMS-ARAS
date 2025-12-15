#-*- coding:gbk -*-
import datetime
import copy
import re

'''
Discription
	
	获取json文件中的sample info，转化为报告模板方便填充的格式。 

'''

def getSample(jsonDict):
	data = copy.deepcopy(jsonDict["sample_info"])
	data["report_date"] = str(datetime.date.today())
	for k, v in data.items():
		if not v:
			data[k] = ""
	
	#新增日期格式：XX年XX月XX日
	receive_date = re.split("-", data["receive_data"])
	data["receive_date_special_1"] = receive_date[0] + "年" + receive_date[1] + "月" + receive_date[2] + "日" if receive_date and len(receive_date) >= 3 else ""
	report_date = re.split("-", data["report_date"])
	data["report_date_special_1"] = report_date[0] + "年" + report_date[1] + "月" + report_date[2] + "日" if report_date and len(report_date) >= 3 else ""

	#新增生信分析时间
	json_name_list = re.split("_", data["json_batch_name"])
	json_date = json_name_list[0] if json_name_list else ""
	data["json_date"] = json_date[:4] + "年" + json_date[4:6] + "月" + json_date[6:] + "日" if json_date and len(json_date) == 8 else ""

	#判断临检样本是厦门接收还是上海接收的，以首字母是否带S进行判断
	data["locate"] = "SH" if re.match("S", data["sample_id"]) else "XM"

	# 浙江妇保，寄送备注处理
	data["ZJFB_mark_dict"] = {}
	if data["mark"]:
		for i in re.split("；|;", data["mark"]):
			tmp = re.split("：|:", i)
			if len(tmp) == 2:
				data["ZJFB_mark_dict"][tmp[0]] = ""
				data["ZJFB_mark_dict"][tmp[0]] = tmp[1]
	
	# 肿瘤细胞含量新增数值字段，便于比对
	data["tumor_content_num"] = data["tumor_content"].replace("%", "") if "tumor_content" in data.keys() and data["tumor_content"] else ""

	# 特殊处理：临检项目，年龄为0的放空，性别为未知的放空-20220923
	if jsonDict["sample_info"]["report_module_type"] != "hospital":
		data["age"] = "" if data["age"] and data["age"] == "0" else data["age"]
		data["gender"] = "" if data["gender"] and data["gender"] == "未知" else data["gender"]

	# 香港大学深圳医院，寄送备注处理-2022.10.10
	data["XGDX_mark_dict"] = {}
	if data["mark"]:
		for i in re.split("，", data["mark"]):
			tmp = re.split("：", i)
			if tmp:
				data["XGDX_mark_dict"][tmp[0]] = tmp[1] if len(tmp) >= 2 else ""
	
	# 采集日期，组织采集日期还是原来的字段：gather_data
	# 新增一个血液采集日期，这边做兼容-2022.11.29
#	data["blood_collection_date"] = data["blood_collection_date"] if "blood_collection_date" in data.keys() and data["blood_collection_date"] else ""
	data["blood_collection_date"] = re.split(" ", data["blood_collection_date"])[0] if "blood_collection_date" in data.keys() and data["blood_collection_date"] else ""

	# 现在日期返回格式有的会带具体时间，这边加个兼容处理-刘炜芬-20231025
	date_key = ["receive_data", "tissue_collection_date", "blood_collection_date", "gather_data", "submission_date", "application_date", "sampling_time", "analysis_date", "section_date", "blood_date_received", "tissue_date_received"]
	for key in date_key:
		if key in data.keys():
			data[key] = re.split(" ", data[key])[0] if data[key] else ""
	# 添加完成-20231025

	# 新增兼容-组织采集日期，目前大部分模板都是使用	gather_data字段，新增了一个tissue_collection_date-刘炜芬，2023.11.14
	# 若存在gather_data，则gather_data使用该字段内容，若不存在gather_data，则使用tissue_collection_date
	data["gather_data"] = data["gather_data"] if "gather_data" in data.keys() and data["gather_data"] else \
						  data["tissue_collection_date"] if "tissue_collection_date" in data.keys() and data["tissue_collection_date"] else \
						  ""
	# 兼容完成-2023.11.14
	
	return data