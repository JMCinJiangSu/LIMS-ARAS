#-*- coding:gbk -*-
import copy

'''
Discription
	
	该脚本用来获取GEP结果。
	提取内容：
	1. 结果
	2. 图片 

'''

def getGEP(jsonDict):

	gep_dict = copy.deepcopy(jsonDict["gep"])
	
	return gep_dict