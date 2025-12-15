#-*- coding:gbk -*-
import copy
from libs import listResultToDict

'''
Discription
	
	该脚本用来获取RNA表达结果，并根据模板需求进行分列处理。 

'''

def getRNA_exp(jsonDict):
	rna_exp_data = copy.deepcopy(jsonDict["rna_exp"])

	# 新增一个兼容-适配v4-2023.11.06
	# tpm为空时，v3返回“0.0”，v4返回变成了空值，反馈过说是要改报告脚本
	for i in rna_exp_data:
		if not i["tpm"]:
			i["tpm"] = "0.0"
	# 兼容完成-2023.11.06
	
	rna_exp = {}
	if rna_exp_data:
		# 临检通用版展示的是四列
		rna_exp["column_4"] = process_data(rna_exp_data, 4) 
		# 浙肿展示五列
		rna_exp["column_5"] = process_data(rna_exp_data, 5)
		# 浙肿需要展示GEP相关基因的表达
		rna_exp["gep"] = gep_tpm(rna_exp_data)
		# xw5101 展示25个基因的tpm
		rna_exp['xw5101'] = process_data(exp_xw5101(rna_exp_data), 2)

	return rna_exp

# num指的是模板中要按几列展示，基因+TPM为一列
def process_data(rna_exp_data, num):
	rna_exp = []
	for i in range(0, len(rna_exp_data)-num, num):
		tmp_dict = {}
		for j in range(1, num + 1):
			tmp_dict["gene"+str(j)] = rna_exp_data[i+j-1]["gene_symbol"]
			tmp_dict["tpm"+str(j)] = rna_exp_data[i+j-1]["tpm"]
		rna_exp.append(tmp_dict)
	
	rest_gene = len(rna_exp_data) % num
	rest_tmp_dict = {}
	for j in range(1, num+1):
		rest_tmp_dict["gene"+str(j)] = ""
		rest_tmp_dict["tpm"+str(j)] = ""

	rest_num = 1
	last_row_num = len(rna_exp_data)-rest_gene if rest_gene != 0 else len(rna_exp_data)-rest_gene-num
	for j in range(last_row_num, len(rna_exp_data)):
		rest_tmp_dict["gene"+str(rest_num)] = rna_exp_data[j]["gene_symbol"]
		rest_tmp_dict["tpm"+str(rest_num)] = rna_exp_data[j]["tpm"]
		rest_num += 1
	rna_exp.append(rest_tmp_dict)

	return rna_exp

# GEP 相关基因的表达量-适用浙肿Master
def gep_tpm(rna_exp_data):
	gene_list = ["CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT"]
	gep_exp = {}
	for i in rna_exp_data:
		if i["gene_symbol"] in gene_list:
			gep_exp.setdefault(i["gene_symbol"].replace("-",""), "")
			gep_exp[i["gene_symbol"].replace("-","")] = i["tpm"]

	return gep_exp

# XW5101 v3 RNAseq只报25个基因
def exp_xw5101(rna_exp_data):
	gene_list = ['BCL2', 'BCL2L1', 'CD14', 'CD34', 'CD38', 'CDK2', 'CDK4', 'CDK6', 'CDK7', 'CDK9', 'FLT3', 'HLA-A', 'HLA-B', 'HLA-C', 'HOXA9', 'ITGAM', 'MCL1', 'MEIS1', 'MEN1', 'MNDA', 'MYC', 'PBX3', 'PRMT5', 'PTPRC', 'TP53']
	exp_xw5101 = []
	for i in rna_exp_data:
		if i['gene_symbol'] in gene_list:
			exp_xw5101.append(i)
	
	return exp_xw5101