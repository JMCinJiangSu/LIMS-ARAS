#-*- coding:gbk -*-
from libs import listResultToDict
import copy
def getHD(jsonDict):
    hd_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["hd"]))
