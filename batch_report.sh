#!/bin/bash

BASE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
echo "当前脚本目录: $BASE_DIR"
INPUT_DIR="${BASE_DIR}/test"
OUTPUT_DIR="${BASE_DIR}/test"
PYTHON_SCRIPT="${BASE_DIR}/main.py"

# 创建输出目录
#mkdir -p "$OUTPUT_DIR"

echo "开始批量处理JSON文件..."
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
python3 "${BASE_DIR}/json_rename.py"
# 计数器
success_count=0
fail_count=0

# 遍历所有JSON文件
for json_file in "$INPUT_DIR"/*.json; do
    if [ -f "$json_file" ]; then
        filename=$(basename "$json_file")
        sample_name="${filename%.json}"  # 移除.json扩展名
        echo "正在处理: $filename (样本名: $sample_name)"
        
        # 生成报告
        python3 "$PYTHON_SCRIPT" -s "$sample_name" -o "$OUTPUT_DIR"
        
        # 检查命令是否执行成功
        if [ $? -eq 0 ]; then
            echo "✓ $filename 处理成功 (样本名: $sample_name)"
            ((success_count++))
        else
            echo "✗ $filename 处理失败"
            ((fail_count++))
        fi
    fi
done

echo "批量处理完成!"
echo "成功: $success_count 个文件"
echo "失败: $fail_count 个文件"

# 显示使用的样本名列表
echo ""
echo "处理的样本名列表 (不含.json后缀):"
for json_file in "$INPUT_DIR"/*.json; do
    if [ -f "$json_file" ]; then
        filename=$(basename "$json_file")
        sample_name="${filename%.json}"
        echo "  - $sample_name"
    fi
done