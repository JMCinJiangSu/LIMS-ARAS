#-*- coding:gbk -*-
import re

'''
Discription

	获取描述中的PMID，并转化为模板可填充的格式。 
	
'''

# 获取描述中的PMID号
def getPMID_from_inter(inter):
	pmid_list = []
	mat = re.compile(r"PMID.\s?\d+")
	for i in mat.findall(str(inter)):
		if re.search(":|: |：|： ", i):
			pmid = (re.split(":|: |：|： ", i))[1].replace(" ", "")
		else:
			pmid = (re.split("PMID", i))[1]
		pmid_list.append(pmid)
	
	return pmid_list

# 获取json中的参考文献集，除了固定文献，其余都从里面获取，没有获取到的为空就行
def getRef_from_json(jsonDict):
	ref_dict = {}
	for ref in jsonDict["refer"]:
		# pmid为空的先不处理-2022.09.22
		if ref["pmid"]:
			pmid = re.split(":", ref["pmid"])[-1].replace("]","")
			ref["authors"] = ref["authors"] if ref["authors"] else ""
			ref["date"] = ref["date"] if ref["date"] else ""
			ref["title"] = ref["title"] if ref["title"] else ""
			ref["journal"] = ref["journal"] if ref["journal"] else ""
			ref["vol"] = ref["vol"] if ref["vol"] else ""
			ref["pmid"] = ref["pmid"] if ref["pmid"] else ""
			ref_dict[pmid] = " ".join([ref["authors"], ref["date"], ref["title"], ref["journal"], ref["vol"], ref["pmid"]])
	
	return ref_dict

# 获取描述中的参考文献
def getRef_from_inter(jsonDict, inter):
	ref_dict = getRef_from_json(jsonDict)
	pmid_list = getPMID_from_inter(inter)
	refer_list = [ref_dict[i] for i in pmid_list if i in ref_dict.keys()]

	return refer_list