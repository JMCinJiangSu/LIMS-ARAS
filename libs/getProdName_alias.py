#-*- coding:gbk -*-
import os
#import xlrd
import pandas as pd 

'''
Discription
	
	获取产品别名配置信息，转化为主要产品名 
	
'''

#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#requirenment_path = os.path.join(BASE_DIR, "config/report_requirenment.xlsx")

def alias_name(config):
	df = pd.read_excel(os.path.join(config, "report_requirenment.xlsx"), sheet_name="prod_alias_name", header=0)
	#print (df)
	#print (df.set_index(df.columns[0]).to_dict()[df.columns[1]])
	return df.set_index(df.columns[0]).to_dict()[df.columns[1]]

#def alias_name(config):
#	requirenment_path = os.path.join(config, "report_requirenment.xlsx")
#	xls = xlrd.open_workbook(requirenment_path)
#	requirenment_sheet = xls.sheet_by_name("prod_alias_name")
#	key = requirenment_sheet.row_values(0)
#	Data = {}
#	for num in range(1, requirenment_sheet.nrows):
#		rows = requirenment_sheet.row_values(num)
#		Data[rows[0]] = rows[1]
#	return Data