#!/usr/bin/env python3
"""
快速清理参考文献格式，使其更符合Vancouver风格。
"""

import re
import sys

def clean_reference(line, ref_num):
    """快速清理单条参考文献"""
    # 移除编号
    line = re.sub(r'^\d+\.\s+', '', line).strip()
    
    # 常见清理操作
    # 1. 确保卷期格式：卷(期)
    line = re.sub(r'(\d+) \(', r'\1(', line)
    
    # 2. 确保页码格式：起始页-结束页
    line = re.sub(r'(\d+)\-(\d+)', r'\1-\2', line)
    
    # 3. 确保年份后使用分号（期刊文章）
    # 查找模式：期刊名. 年份;卷(期):页码
    if re.search(r'\.\s+\d{4}[,\s]', line):
        line = re.sub(r'\.\s+(\d{4})[,\s]', r'. \1;', line, count=1)
    
    # 4. 修复双分号问题
    line = line.replace(';;', ';')
    
    # 5. 清理多余空格
    line = re.sub(r'\s+', ' ', line).strip()
    
    # 6. 确保句点结尾
    if not line.endswith('.'):
        line += '.'
    
    # 7. 特殊处理：WHO报告格式
    if 'World Health Organization' in line and 'Geneva:' in line:
        # 已经是正确的书籍/报告格式
        pass
    
    return f"{ref_num}. {line}"

def main():
    if len(sys.argv) < 2:
        print("使用方法: python quick_clean_references.py <输入文件>")
        return
    
    input_file = sys.argv[1]
    
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    in_references = False
    formatted_lines = []
    ref_count = 0
    
    for line in lines:
        stripped = line.strip()
        
        # 检测参考文献部分开始
        if stripped.startswith('# References'):
            in_references = True
            formatted_lines.append(line)
            continue
        
        # 跳过格式说明行但保留它
        if in_references and stripped.startswith('*Note:'):
            formatted_lines.append('\n')  # 替换为空白行
            continue
        
        # 检测参考文献部分结束
        if in_references and stripped.startswith('---') and 'Word count' in line:
            in_references = False
            formatted_lines.append(line)
            continue
        
        if in_references and stripped and not stripped.startswith('*'):
            # 是参考文献行
            ref_count += 1
            cleaned = clean_reference(stripped, ref_count)
            formatted_lines.append(cleaned + '\n')
        else:
            formatted_lines.append(line)
    
    # 写回文件
    with open(input_file, 'w', encoding='utf-8') as f:
        f.writelines(formatted_lines)
    
    print(f"已快速清理 {ref_count} 条参考文献")
    print(f"文件已保存: {input_file}")
    print(f"注意：这是快速清理，建议人工检查关键条目")

if __name__ == '__main__':
    main()