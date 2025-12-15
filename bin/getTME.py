#-*- coding:gbk -*-
import copy
from libs import listResultToDict
from libs.rule import decimal_float, is_number

'''
Discription
	
	该脚本用来获取TME结果。
	提取内容：
	1. 结果
	2. 图片 

'''

def getTME(jsonDict):
	tme_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["tme_type"])) if "tme_type" in jsonDict.keys() else {}
	tme_score = copy.deepcopy(jsonDict["tme_score"]) if "tme_score" in jsonDict.keys() else {}
	if tme_score:
		for k,v in tme_score.items():
			#tme_score[k] = "{:.2f}".format(float(v))
			#tme_score[k] = decimal_float(v)
			# tme_score返回中新增了不是数值型的字段（sample_id等），加个兼容-2023.09.20
			tme_score[k] = decimal_float(v) if is_number(v) else v
			# 兼容完成-2023.09.20

	return tme_dict, tme_score