import os
import json
import re

def rename_json_files(folder_path):
    """
    批量重命名文件夹中的JSON文件，使用文件内 sample_info.order_id 作为新文件名
    
    :param folder_path: 包含JSON文件的文件夹路径
    """
    # 检查文件夹是否存在
    if not os.path.exists(folder_path):
        print(f"错误：文件夹 '{folder_path}' 不存在")
        return
    
    # 计数器
    renamed_count = 0
    error_count = 0
    
    # 遍历文件夹中的所有文件
    for filename in os.listdir(folder_path):
        if not filename.lower().endswith('.json'):
            continue
            
        file_path = os.path.join(folder_path, filename)
        
        try:
            # 读取JSON文件
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # 提取order_id
            order_id = data['sample_info']['order_id']
            prods_name = data['sample_info']['prod_names']
            prods_name = re.sub(r'[^\w\-]+', '_', prods_name)  # 替换非法字符
            
            # 构建新文件名
            new_filename = f"{order_id}_{prods_name}.json"
            new_path = os.path.join(folder_path, new_filename)
            
            # 避免覆盖现有文件
            if os.path.exists(new_path):
                print(f"警告：文件 '{new_filename}' 已存在，跳过重命名 '{filename}'")
                error_count += 1
                continue
            
            # 重命名文件
            os.rename(file_path, new_path)
            print(f"重命名: {filename} -> {new_filename}")
            renamed_count += 1
            
        except Exception as e:
            print(f"处理文件 '{filename}' 时出错: {str(e)}")
            error_count += 1
    
    # 输出摘要
    print(f"\n操作完成: 成功重命名 {renamed_count} 个文件, 遇到 {error_count} 个错误")

if __name__ == "__main__":
    # 设置你的文件夹路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_folder = os.path.join(script_dir, "test")
    
    # 执行重命名
    rename_json_files(target_folder)