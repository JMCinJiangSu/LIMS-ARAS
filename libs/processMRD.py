#-*- coding:gbk -*-
import copy
from libs.rule import decimal_float

def process_MRD(jsonDict):
    mrd = copy.deepcopy(jsonDict['mrd']) if "mrd" in jsonDict.keys() else []
    mrd_result = []
    if mrd:
        for item in mrd:
            item['mrd_status'] = '阳性' if item['mrd_status'] == 'Positive' else '阴性'
            item['receive_date'] = item['receive_date'].split()[0].replace('-', '.')
            item['ctdna_conc'] = decimal_float(item['ctdna_conc'])
            mrd_result.append(item)
    
    mrd_result = sorted(mrd_result, key=lambda i :i['receive_date'])

    
    if mrd_result:
        for i in range(len(mrd_result)):
            mrd_result[i]['num'] = i + 1
    
    latest_result = mrd_result[-1] if mrd_result else ''

    return mrd_result, latest_result
