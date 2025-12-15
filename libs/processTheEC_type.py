#-*- coding:gbk -*-
import re
from libs.getEvi import varRegimen
from libs import listResultToDict
import copy

'''
Discription
	
	处理子宫内膜癌分子分型格式。 

'''

def process_ec_type(jsonDict, config):
	ec_type_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["ec_type"]))
	ec_type_dict["evi_sum"] = varRegimen(jsonDict, ec_type_dict["evi_sum"], config, ec_type_dict)
	
	return ec_type_dict