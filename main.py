#-*- coding:gbk -*-
import argparse
import os
import re
#import subprocess as sp
import AutoReport_Base_Code_v1_0_0


BASE_DIR = (AutoReport_Base_Code_v1_0_0.__file__).replace("AutoReport_Base_Code_v1_0_0.py", "")


config = os.path.join(BASE_DIR, "config")

report_template = os.path.join(BASE_DIR, "report_template")

base_script = os.path.join(BASE_DIR, "AutoReport_Base_Code_v1_0_0.py")

outjson = "F"

image = os.path.join(BASE_DIR, "image")



def main(json_name, outfile):
	
	AutoReport_Base_Code_v1_0_0.get_data(json_name, outfile, config, report_template, outjson, image)

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--json_name", dest = "json_name", required = True)
	parser.add_argument("-o", "--outfile", dest = "outfile", required = True)
	arg = parser.parse_args()

	return arg

def script_main_auto_output_report_file(json_name=None, outfile_dir=None):

	out_report_file, word_name = main(json_name=json_name, outfile=outfile_dir)
	return out_report_file, word_name

if __name__ == '__main__':
	args = parse_args()
	main(json_name=args.json_name, outfile=args.outfile)