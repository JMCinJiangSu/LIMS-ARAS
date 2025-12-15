#-*- coding:gbk -*-

'''
Discription
	
	将MSI、PD-L1、TMB等应该为字典却输出为列表的信息转回字典，若有多个，则默认选择列表第一个元素 
	
'''

def ListToDict(json_result):
	return json_result if type(json_result).__name__=="dict" else json_result[0] if json_result else {}