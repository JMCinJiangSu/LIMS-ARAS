#-*- coding:gbk -*-
import os
import re
import jinja2
from docxtpl import DocxTemplate
import datetime
import io
#import threading
import customize_filters

jinja_env = jinja2.Environment()

'''
Discription
	
	处理基础模板，分为固定模板和动态填充模板

'''

def BaseReport(data, tpl, merge_template, report_template, json_name, outfile):
	dynamic_dir = os.path.join(report_template, "template_dynamic")
	fixed_dir = os.path.join(report_template, "template_fixed")
	data["base_report"] = {}

	result = {}

	## 动态填充模板
	# 需要拼接，merge_template == 0
	if merge_template["merge_template"] == "0":
		temp_dynamic_list = []
		# 动态模板，指定列表，为no时则不处理
		if "no" not in merge_template["temp_dynamic"]:
			# 汇总动态模板和对应的样式
			if merge_template["temp_dynamic"]:
				# 格式：module_name.type，若只有module_name的话，type赋为默认值default
				temp_dynamic_list = [{"module_name" : re.split("\.", i)[0], "type" : re.split("\.", i)[1] if len(re.split("\.", i)) == 2 else "default"} for i in merge_template["temp_dynamic"]]
			# 动态模板，默认，需要遍历文件夹，模板样式也为默认的
			else:
				temp_dynamic_list = [{"module_name" : re.split("\.", i)[0], "type" : ""} for i in os.listdir(dynamic_dir)]
		
		for temp in temp_dynamic_list:
			# 将动态模板和对应样式ID返回到data数据结构中，用于指定样式选择
			data["base_report"][temp["module_name"]+"_type"] = temp["type"]
#			print (data["base_report"])
			# 动态模板填充
			tpl_name = DocxTemplate(os.path.join(dynamic_dir, temp["module_name"]+".docx"))
			tpl_name.render(data, jinja_env)
			# 使用IO BytesIO写到内存中，不再输出到磁盘中
			file_stream = io.BytesIO()
			tpl_name.save(file_stream)
			result[temp["module_name"]] = tpl.new_subdoc(file_stream)
#			tpl_name.save("./tmp/"+temp["module_name"]+"subdoc.docx")
			#print (temp, datetime.datetime.now())
		

		'''
		# 多线程-测试结果：速度没有提升
		def newsub(temp):
			# 将动态模板和对应样式ID返回到data数据结构中，用于指定样式选择
			data["base_report"][temp["module_name"]+"_type"] = temp["type"]
			# 动态模板填充
			tpl_name = DocxTemplate(os.path.join(dynamic_dir, temp["module_name"]+".docx"))
			tpl_name.render(data, jinja_env)
			# 使用IO BytesIO写到内存中，不再输出实体文件
			file_stream = io.BytesIO()
			tpl_name.save(file_stream)
			result[temp["module_name"]] = tpl.new_subdoc(file_stream)
			print (temp, datetime.datetime.now())
		
		for temp in temp_dynamic_list:
			t = threading.Thread(target = newsub, args = (temp, ))
			t.start()
			print (threading.active_count())
		# 测试多线程
		'''

		# 固定列表需要评估下是否要做指定列表，按默认的处理运行速度差别如何？
		# 固定模板，指定列表，为no时不处理
		if merge_template["temp_fixed"]:
			temp_fixed_list = merge_template["temp_fixed"] if "no" not in merge_template["temp_fixed"] else []
		else:
			temp_fixed_list = [re.split("\.", i)[0] for i in os.listdir(fixed_dir)]
		
		for temp in temp_fixed_list:
			result[temp] = tpl.new_subdoc(fixed_dir+"/"+temp+".docx")

	return result