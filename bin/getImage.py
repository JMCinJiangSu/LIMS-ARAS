#-*- coding:gbk -*-
from docxtpl import InlineImage
from docx.shared import Mm
import os
import re
from libs.getImageSize import image_file

'''
Discription
	
	该脚本用来处理填充的图片，主要处理使用插入方法（InlineImage）的图片

'''
def render_image(tpl, data, jsonDict, report_name, image, config):
	result = {}
	# 插入图片
	# 1. BRCA IGV图
	result["brca_all"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["igvplot"]["all"] if i != None and os.path.exists(i)]
#	print (result["brca_all"])
	result["brca_g"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["igvplot"]["germline"] if i != None and os.path.exists(i)]
	# 2. MSI图
	if data["msi"] and "img_path" in data["msi"].keys() and data["msi"]["img_path"] and os.path.exists(data["msi"]["img_path"]):
		result["msi"] = InlineImage(tpl, data["msi"]["img_path"], width=Mm(90))
	# 3. TME图
	if data["qc"] and "rna_data_qc" in data["qc"].keys() and data["qc"]["rna_data_qc"] and "tmeplot" in data["qc"]["rna_data_qc"].keys() and \
		data["qc"]["rna_data_qc"]["tmeplot"] and os.path.exists(data["qc"]["rna_data_qc"]["tmeplot"]):
		result["tme"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["tmeplot"], width=Mm(110))
	# 4. GEP图 先判断gep.img_path，有的话直接展示，无则判断qc.rna_data_qc.gepplot
	if data["gep"] and "img_path" in data["gep"].keys() and data["gep"]["img_path"] and os.path.exists(data["gep"]["img_path"]):
		result["gep"] = InlineImage(tpl, data["gep"]["img_path"], width=Mm(90))
	elif data["qc"] and "rna_data_qc" in data["qc"].keys() and data["qc"]["rna_data_qc"] and "gepplot" in data["qc"]["rna_data_qc"].keys() and \
		data["qc"]["rna_data_qc"]["gepplot"] and os.path.exists(data["qc"]["rna_data_qc"]["gepplot"]):
		result["gep"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["gepplot"], width=Mm(90))
	# 5. 新增I/II类变异Snvindel IGV图-2022.10.20
	result["igv_I_II"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_I_II"] if os.path.exists(i)]
	# 新增适用福建附一CP40的IGV图，长11cm，宽5cm-2023.10.20
	result["igv_I_II_FJFY"] = [InlineImage(tpl, i, width=Mm(110), height=Mm(50)) for i in data["var"]["igv_I_II"] if os.path.exists(i)]
	# 6. 新增4/5类变异Snvindel IGV图-2022.11.18
	result["igv_4_5"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_4_5"] if os.path.exists(i)]
	# 7. TMB图
	if data["tmb"] and "img_path" in data["tmb"].keys() and os.path.exists(data["tmb"]["img_path"]):
		result["tmb"] = InlineImage(tpl, data["tmb"]["img_path"], width=Mm(100))
	# 8. PD-L1图
	if data["pdl1"] and "file_pdl1" in data["pdl1"].keys() and os.path.exists(data["pdl1"]["file_pdl1"]):
		result["pdl1"] = InlineImage(tpl, data["pdl1"]["file_pdl1"], width=Mm(80))
	# 9. CNV图
	if "cnv_file_path" in jsonDict.keys() and jsonDict["cnv_file_path"] and "abs_path" in jsonDict["cnv_file_path"].keys() and \
		jsonDict["cnv_file_path"]["abs_path"] and os.path.exists(jsonDict["cnv_file_path"]["abs_path"]):
		result["cnv"] = InlineImage(tpl, jsonDict["cnv_file_path"]["abs_path"], width=Mm(160))
	# 10. MLPA图（仅del）
	result["mlpa_image_del"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["mlpa_image_del"] if os.path.exists(i)]
	# 11. MRD图
	if data['mrd_last'] and 'mrd_img_path' in data['mrd_last'].keys() and os.path.exists(data['mrd_last']['mrd_img_path']):
		result['mrd'] = InlineImage(tpl, data['mrd_last']['mrd_img_path'], width=Mm(180))

	# 新增肿瘤发生发展相关Snindel IGV图-2023.06.27
	result["igv_onconodrug"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_onconodrug"] if os.path.exists(i)]
	# 新增结束-2023.06.27

	# 获取配置表中图片指定尺寸，默认为100
	image_sise_dict = image_file(config)
	# 模板中固定的图片，直接拼接会显示不出来，这边把固定图片按插入的形式进行填充
	image_list = [i for i in os.listdir(image)]
	for i in image_list:
		image_name = re.split("\.", i)[0]
		width = image_sise_dict.get(image_name, 100) 
		result["fixed_"+image_name] = InlineImage(tpl, image+"/"+i, width=Mm(width))


	# 插入或替换图片
	# 7. TMB图
#	if data["tmb"] and "img_path" in data["tmb"].keys() and os.path.exists(data["tmb"]["img_path"]):
#		# 报告打开报错，TMB图插入方式修改-仅临检通用Master样本（2022.09.27新增复旦中山）
#		if re.search("Master|90|87", data["sample"]["prod_names"]) and re.search("rummage|FDZS|SDSL", report_name):
#			tpl.replace_pic("test_TMB.jpg", data["tmb"]["img_path"])
#			print ("replace")
#		else:
#			result["tmb"] = InlineImage(tpl, data["tmb"]["img_path"], width=Mm(100))
#			print (inline)

	# 替换图片
	# 8. PD-L1图 图片太多可能会报错，PDL1用替换
#	if data["pdl1"] and "file_pdl1" in data["pdl1"].keys() and os.path.exists(data["pdl1"]["file_pdl1"]):
#		tpl.replace_pic("test_PDL1.jpg", data["pdl1"]["file_pdl1"])
	# 9. CNV图
	# 曲线贴图-这个原因待查找！！！
	# 使用插入图片的方法，本地测试正常、线上运行打开word会报错，可能的原因有1）word格式问题，2）模块版本差异；排查太耗时了，暂时放弃
	# 改用replace方法，插入图片名跟本地不符合，代码找不到，排查耗时，先不找了。暂时使用之前正常替换图片的图片模板，成功运行，先用着吧，有空再来排查
#	if re.search("Pan116（组织）|LC76（组织）|CRC25（组织）|GA18（组织）|TC21（组织）", data["sample"]["prod_names"]) and re.search("rummage", report_name):
#		if "cnv_file_path" in jsonDict.keys() and jsonDict["cnv_file_path"] and "abs_path" in jsonDict["cnv_file_path"].keys() and jsonDict["cnv_file_path"]["abs_path"] and os.path.exists(jsonDict["cnv_file_path"]["abs_path"]):
#			tpl.replace_pic("test_MSI.png", jsonDict["cnv_file_path"]["abs_path"])
	
	return result
