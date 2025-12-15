#-*- coding:gbk -*-
import copy
from libs import listResultToDict
'''
Discription
	
	该脚本用来获取PD-L1结果。
	提取内容：
	1. 结果，阴性/阳性（肺癌TPS/其他癌种CPS）
	2. 图片
	3. 治疗方案（未确定需不需要提取）
	用于PD-L1特殊处理内容 

'''

def getPDL1(jsonDict):
	pdl1_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["pdl1"]))
	if pdl1_dict:
		pdl1_dict["value"] = str(int(float(pdl1_dict["value"])*100))+"%" if pdl1_dict["result"] == "阳性" and pdl1_dict["type"] == "TPS" else int(pdl1_dict["value"]) if pdl1_dict["result"] == "阳性" and pdl1_dict["type"] == "CPS" else pdl1_dict["value"]

	return pdl1_dict