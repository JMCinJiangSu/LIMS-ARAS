
import re
import os
import json
from bin import getVar

def get_related_var(json_name, outfile, config, report_name):
    data = []
    # 解析json文件
    with open(os.path.join(outfile, json_name + ".json"), "r", encoding='utf-8') as file_json:
        jsonDict = json.load(file_json)
        current_filepath = os.path.join(outfile, json_name + ".json")
        # 系统将rummage改为clinical了，为了报告脚本代码改动最小化，在这边做转化
        jsonDict["sample_info"]["report_module_type"] = "rummage" if jsonDict["sample_info"]["report_module_type"] == "clinical" else jsonDict["sample_info"]["report_module_type"]
        tmp_dict, _, _ = getVar.getVar(jsonDict, config, report_name)
        current_var = get_var(tmp_dict)
        data.append({
            "mrd_date" : jsonDict["sample_info"]["receive_data"].split()[0] if jsonDict["sample_info"]["receive_data"] else "",
            "var" : current_var if current_var else [],
            "json_name" : json_name,
            'var_count' : len(current_var) if current_var else 0
        })
    
    # 获取关联的json
    if "relation_order_code" in jsonDict["sample_info"].keys() and jsonDict["sample_info"]["relation_order_code"]:
        for name in jsonDict["sample_info"]["relation_order_code"]:
            related_json_path = os.path.join(os.path.dirname(current_filepath), f"{name}.json")
            with open(related_json_path, "r", encoding='utf-8') as file_json:
                related_json = json.load(file_json)
                mrd = related_json["mrd"]
                mrd_date = [item["receive_date"].split()[0] for item in mrd] if mrd else []
                tmp_dict, _, _ = getVar.getVar(related_json, config, report_name)
                result = get_var(tmp_dict)
                data.append({
                    "mrd_date" : mrd_date[-1] if mrd_date else "",
                    "var" : result if result else [],
                    "json_name" : name,
                    'var_count' : len(result) if result else 0
                })
    data = sorted(data, key=lambda x: x["mrd_date"]) if data else []
    if data:
        for i in range(len(data)):
            data[i]["num"] = i + 1

    return data

def get_var(tmp_dict):
    var_list = tmp_dict["var_somatic"]["level_I"] + tmp_dict["var_somatic"]["level_II"] + tmp_dict["var_somatic"]["level_onco_nodrug"] + tmp_dict["var_somatic"]["level_III"]
    # MRD应该只看Snvindel的结果
    var_list = [var for var in var_list if var["bio_category"] == "Snvindel"]
    result = []
    if var_list:
        for var in var_list:
            result.append(
                {
                    'gene_symbol' : var['gene_symbol'],
                    'gene_region' : var['gene_region'],
                    'hgvs_c' : var['hgvs_c'],
                    'hgvs_p' : var['hgvs_p'],
                    'hgvs_p_abbr' : var['hgvs_p_abbr'],
                    'transcript' : var['transcript_primary'],
                    'freq' : var['freq_str']
                }
            )
    return result