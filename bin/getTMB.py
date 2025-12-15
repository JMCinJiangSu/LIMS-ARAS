#-*- coding:gbk -*-
import copy
from libs import listResultToDict
from libs.getEvi import varRegimen

'''
Discription
	
	该脚本用来获取TMB结果。
	提取内容：
	1. 结果
	2. 图片
	3. 治疗方案
	用于TMB特殊处理内容 

'''

def getTMB(jsonDict, config):
	tmb_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["tmb"]))
	if 'TMB_value' in tmb_dict.keys() and tmb_dict['TMB_value']:
		tmb_dict['TMB_value'] = round(float(tmb_dict['TMB_value']), 2)
	if "evi_sum" in tmb_dict.keys() and tmb_dict["evi_sum"]:
		tmb_dict["evi_sum"] = varRegimen(jsonDict, tmb_dict["evi_sum"], config, tmb_dict)

	return tmb_dict