# -*- coding: utf-8 -*-
"""
SampleInfo.py - 带详细日志的 SQL Server 数据提取脚本
输出：LIMSResult.json + LIMS结果-YYYY-MM-DD.xlsx + SampleInfo.log
"""
from typing import Any, Dict
import pyodbc
import traceback
from timeit import default_timer
import os
from collections import ChainMap
from openpyxl import Workbook
from datetime import date, datetime
import json
import logging


# === 全局配置 ===
base_dir = os.getcwd()

host_sql = '192.168.2.21'
usr_sql = 'amoydxreader'
passwd_sql = 'Adx@361027'
db = 'ELIMSDB'


# === 日志配置 ===
log_file = os.path.join(base_dir, "SampleInfo.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    encoding="utf-8"
)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, datetime):
            return o.isoformat()
        return super().default(o)


class UsingMssql(object):
    """上下文管理的 SQL Server 连接类（带日志）"""

    def __init__(self, commit=False, log_time=True, log_label='SQL执行耗时'):
        self._log_time = log_time
        self._commit = commit
        self._log_label = log_label

    def __enter__(self):
        try:
            if self._log_time:
                self._start = default_timer()

            conn_str = (
                "DRIVER={ODBC Driver 17 for SQL Server};"
                f"SERVER={host_sql};"
                f"DATABASE={db};"
                f"UID={usr_sql};"
                f"PWD={passwd_sql};"
                "TrustServerCertificate=yes;"
            )
            self._conn = pyodbc.connect(conn_str, timeout=10)
            self._cursor = self._conn.cursor()
            logging.info(f"✅ 成功连接到数据库：{host_sql} / {db}")
        except Exception as e:
            logging.error("❌ 数据库连接失败: %s", str(e))
            logging.error(traceback.format_exc())
            raise
        return self

    def __exit__(self, *exc_info):
        try:
            if self._commit:
                self._conn.commit()
            self._cursor.close()
            self._conn.close()

            if self._log_time:
                diff = default_timer() - self._start
                logging.info("-- %s: %.6f 秒", self._log_label, diff)
        except Exception as e:
            logging.error("关闭数据库连接时出错: %s", str(e))

    @property
    def cursor(self):
        return self._cursor

'''
def mydecode_dic(dic):
    if not dic or not isinstance(dic, dict):
        return dic
    res = dict()
    for item in dic.items():
        try:
            key = item[0].encode('latin1').decode('gbk') if isinstance(item[0], str) else item[0]
            val = item[1].encode('latin1').decode('gbk') if isinstance(item[1], str) else item[1]
            res[key] = val
        except Exception:
            logging.warning("⚠️ 解码异常字段: %s", str(item))
    return res
'''

def mydecode_dic(dic: Dict[Any, Any]) -> Dict[Any, Any]:
    """
    解码字典中的中文字符串（从latin1错误解码恢复为GBK/UTF-8）
    
    Args:
        dic: 需要解码的字典
        
    Returns:
        解码后的字典
    """
    if not dic or not isinstance(dic, dict):
        return dic
    
    res = {}
    for key, value in dic.items():
        try:
            # 解码key
            decoded_key = key
            if isinstance(key, str):
                decoded_key = decode_chinese_string(key)
            
            # 解码value
            decoded_value = value
            if isinstance(value, str):
                decoded_value = decode_chinese_string(value)
            elif isinstance(value, (list, tuple)):
                # 处理嵌套的列表/元组
                decoded_value = type(value)(decode_chinese_string(item) if isinstance(item, str) else item 
                                          for item in value)
            elif isinstance(value, dict):
                # 递归处理嵌套字典
                decoded_value = mydecode_dic(value)
            
            res[decoded_key] = decoded_value
            
        except Exception as e:
            logging.warning("⚠️ 解码异常字段: key=%s, value=%s, error=%s", 
                          str(key), str(value), str(e))
            # 出错时保留原始值
            res[key] = value
    
    return res

def decode_chinese_string(s: str, encodings: list = None) -> str:
    """
    解码中文字符串，尝试多种编码方式
    
    Args:
        s: 需要解码的字符串
        encodings: 尝试的编码列表，默认['gbk', 'utf-8', 'gb2312']
        
    Returns:
        解码后的字符串
    """
    if encodings is None:
        encodings = ['gbk', 'utf-8', 'gb2312']
    
    # 如果字符串看起来已经是正常中文，直接返回
    if is_chinese_string(s):
        return s
    
    # 方法1: latin1 -> 目标编码 (修复错误解码)
    for encoding in encodings:
        try:
            return s.encode('latin1').decode(encoding)
        except (UnicodeEncodeError, UnicodeDecodeError):
            continue
    
    # 方法2: 尝试直接转换（针对ASCII码形式的字符串）
    if is_ascii_code_string(s):
        try:
            return decode_ascii_codes(s, encodings[0])  # 默认使用第一个编码
        except Exception:
            pass
    
    # 方法3: 尝试其他可能的编码组合
    for encoding in encodings:
        try:
            return s.encode(encoding).decode(encoding)  # 尝试直接解码
        except (UnicodeEncodeError, UnicodeDecodeError):
            continue
    
    # 所有方法都失败，返回原始字符串
    logging.debug("无法解码字符串: %s", s)
    return s

def is_chinese_string(s: str) -> bool:
    """检查字符串是否包含中文字符"""
    import re
    chinese_pattern = re.compile(r'[\u4e00-\u9fff]')
    return bool(chinese_pattern.search(s))

def is_ascii_code_string(s: str) -> bool:
    """检查字符串是否由ASCII码数字组成"""
    return all(c in '0123456789 ' for c in s.strip())

def decode_ascii_codes(ascii_str: str, encoding: str = 'gbk') -> str:
    """将ASCII码字符串转换为中文"""
    ascii_list = ascii_str.strip().split()
    bytes_data = bytes(int(x) for x in ascii_list if 0 <= int(x) <= 255)
    return bytes_data.decode(encoding)


class ReportData:
    def __init__(self, sample_id):
        self.sample_id = sample_id

    def get_sample_info(self):
        data_js = dict()
        sample_id = self.sample_id

        sql = f"""
        SELECT FilterCode, SSZCode, VisitName, Birthday, CreateDate,
               PathologyDiagnose, SiteofBiopsyCollection, SurgerySource, Extend1,
               CONVERT(VARCHAR(10), SampleCollectionTime, 23) AS BloodCollectionTime,
               CONVERT(VARCHAR(10), BodySampleCollectionTime, 23) AS TissueCollectionTime,
               CONVERT(VARCHAR(10), CreateDate, 23) AS CreateDate,
               CONVERT(VARCHAR(10), PreparationTime, 23) AS PreparationTime
        FROM [ProjectManger].[ProjectInfo]
        WHERE MainBarcode = '{sample_id}'
        """

        try:
            with UsingMssql(log_time=False) as um:
                um.cursor.execute(sql)
                desc = [col[0] for col in um.cursor.description]
                rows = um.cursor.fetchall()
                if not rows:
                    logging.warning(f"样本 {sample_id} 未在 ProjectInfo 表中找到记录。")
                    return {}
                data_js = [dict(zip(desc, rows[0]))][0]
        except Exception as e:
            logging.error(f"查询 ProjectInfo 失败 ({sample_id}): {e}")
            logging.error(traceback.format_exc())
            return {}

        sql = f"""
        SELECT HospitalName, SexName, SampleTypeName, CancerType, CuttingPart, CONVERT(VARCHAR(10), SampleReceiveDate, 23) AS ReceiveDate
        FROM [PerWork].[ApplySampleInfo]
        WHERE Barcode = '{sample_id}'
        """

        try:
            with UsingMssql(log_time=False) as um:
                um.cursor.execute(sql)
                descrip = [col[0] for col in um.cursor.description]
                tmp = um.cursor.fetchall()
                if tmp:
                    tmp_dict = [dict(zip(descrip, tmp[0]))][0]
                    data_js = dict(ChainMap(data_js, tmp_dict))
        except Exception as e:
            logging.error(f"查询 ApplySampleInfo 失败 ({sample_id}): {e}")
            logging.error(traceback.format_exc())

        data_js['Birthday'] = str(data_js.get('Birthday', ''))[:4] if data_js.get('Birthday') else ''
        return mydecode_dic(data_js)

    def lib_qc(self):
        sample_id = self.sample_id
        sql = f"""
        SELECT DISTINCT DataKeyFiled, DataValue
        FROM WF.WFRuningRecordResultParse
        WHERE BusinessParam1 = '{sample_id}'
        AND DataKeyFiled IN (
            '肿瘤细胞含量', '是否富集', 'DNA浓度(ng/μL)', 'DNA总量(ng)',
            'RNA浓度(ng/μL)', 'RNA总量(ng)', 'RNA DV200(%)',
            '文库总量(ng)', '捕获文库总量(ng)', '文库浓度(ng/μL)', '捕获文库浓度(ng/μL)'
        )
        """
        try:
            with UsingMssql(log_time=False) as um:
                um.cursor.execute(sql)
                rows = um.cursor.fetchall()
                if not rows:
                    logging.info(f"样本 {sample_id} 在 WFRuningRecordResultParse 无 QC 数据。")
                    return {}
                data_qc = dict(rows)
            return mydecode_dic(data_qc)
        except Exception as e:
            logging.error(f"查询 QC 数据失败 ({sample_id}): {e}")
            logging.error(traceback.format_exc())
            return {}


def gen_report(sample_id):
    logging.info(f"▶️ 开始处理样本 {sample_id}")
    data = ReportData(sample_id)
    info = data.get_sample_info()
    library_qc = data.lib_qc()
    merged = {**info, **library_qc}
    logging.info(f"✅ 样本 {sample_id} 数据提取完成 ({len(merged)} 字段)")
    return merged


def Samplelist(f):
    with open(f, 'r', encoding='utf-8') as file:
        lines = file.read().split('\n')
    return [x.strip() for x in lines if x.strip()]


def dict_to_excel(result):
    wb = Workbook()
    ws = wb.active
    ws.title = 'LIMS结果'
    filename = f'LIMS结果-{datetime.now().strftime("%Y-%m-%d %H-%M-%S")}.xlsx'

    if not result:
        logging.warning("⚠️ 没有数据可写入 Excel。")
        return

    # 写入顶部说明
    gen_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    ws.append([f"生成时间：{gen_time}"])
    ws.append([f"样本总数：{len(result)}"])
    ws.append([])  # 空行分隔

    # 获取第一个样本的字段名
    first_sample = next(iter(result))
    colnames = ['sampleid'] + list(result[first_sample].keys())
    ws.append(colnames)

    # 写入所有样本数据
    for sampleid, info in result.items():
        row_content = [sampleid] + [info.get(col, '') for col in result[first_sample].keys()]
        ws.append(row_content)

    # 自动调整列宽（可选，美观）
    for col in ws.columns:
        max_length = max(len(str(cell.value)) if cell.value else 0 for cell in col)
        adjusted_width = (max_length + 2)
        ws.column_dimensions[col[0].column_letter].width = adjusted_width

    wb.save(filename)
    logging.info(f"📄 Excel 文件已生成（含时间与样本说明）：{filename}")



if __name__ == '__main__':
    logging.info("=== 启动 LIMS 报告数据提取任务 ===")
    result = dict()
    try:
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                if 'sampleid.txt' not in file:
                    continue
                path = os.path.join(root, file)
                sample_list = Samplelist(path)
                for sid in sample_list:
                    info = gen_report(sid)
                    result[sid] = info

        #datajson = json.dumps(result, ensure_ascii=False, cls=DateTimeEncoder)
        #with open('LIMSResult.json', "w", encoding="utf-8") as outFile:
        #    outFile.write(datajson)
        #logging.info("✅ JSON 文件生成完成：LIMSResult.json")

        if result:
            dict_to_excel(result)
    except Exception as e:
        logging.error("脚本执行过程中出现异常: %s", str(e))
        logging.error(traceback.format_exc())
    finally:
        logging.info("=== LIMS 报告数据提取任务结束 ===\n")
