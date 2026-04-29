#!/usr/bin/env python3
"""格式化考试资料，生成整洁的排版"""

import re

# 读取原始内容
with open('/tmp/exam_full.txt', 'r', encoding='utf-8', errors='ignore') as f:
    lines = f.readlines()

# 处理每一行，添加适当的格式
formatted_lines = []
current_section = ""

for line in lines:
    line = line.strip()
    if not line:
        # 空行保持空行（段落的分隔）
        formatted_lines.append("")
        continue
    
    # 检测章节标题
    if line.startswith('<') and line.endswith('>'):
        section_name = line.strip('<>')
        current_section = section_name
        formatted_lines.append("")
        formatted_lines.append(f"## {section_name}")
        formatted_lines.append("")
        continue
    
    # 检测数字编号的条目（如 "1. xxx" 或 "1、xxx"）
    if re.match(r'^\d+[.、]', line):
        formatted_lines.append(line)
        continue
    
    # 普通文本，保持原样
    formatted_lines.append(line)

# 写入格式化后的文件
output_path = '/root/.openclaw/workspace/exam_materials/眼科学硕士中期考核复习资料_整理版.md'

with open(output_path, 'w', encoding='utf-8') as f:
    f.write("# 眼科学硕士中期考核 - 完整复习资料\n\n")
    f.write("**整理日期**: 2026-04-15\n\n")
    f.write("---\n\n")
    
    # 添加目录
    f.write("## 目录\n\n")
    f.write("1. 名词解释\n")
    f.write("2. 问答题\n")
    f.write("3. 病例分析\n")
    f.write("4. 历年真题\n\n")
    f.write("---\n\n")
    
    # 写入正文，每个段落之间加空行
    for i, line in enumerate(formatted_lines):
        # 在章节标题后添加空行
        if line.startswith("## "):
            f.write("\n" + line + "\n\n")
        # 在普通段落之间添加空行
        elif line and not line.startswith("##"):
            f.write(line + "\n\n")
        else:
            f.write(line + "\n")

print(f"格式化完成: {output_path}")
print(f"总行数: {len(formatted_lines)}")