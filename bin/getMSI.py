#-*- coding:gbk -*-
import copy
from libs import listResultToDict
from libs.getEvi import varRegimen

'''
Discription
	
	该脚本用来获取MSI结果。
	提取内容：
	1. 结果
	2. 图片
	3. 治疗方案
	用于MSI特殊处理内容 

'''

def getMSI(jsonDict, config):
	msi_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["msi"]))
	# MSI_num json无返回，通过脚本计算
	if msi_dict:
		msi_dict["msi_num_cp40"] = int(float(msi_dict["msi_score"]) * 55 + 0.5) if msi_dict["msi_score"] or msi_dict["msi_score"]==0 else ""
	if "evi_sum" in msi_dict.keys() and msi_dict["evi_sum"]:
		msi_dict["evi_sum"] = varRegimen(jsonDict, msi_dict["evi_sum"], config, msi_dict)
	return msi_dict