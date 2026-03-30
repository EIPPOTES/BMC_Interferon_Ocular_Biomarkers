#!/usr/bin/env python3
"""
最终节编号修复脚本
确保所有3.x节编号连续且逻辑正确
"""

import re

def read_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        return f.read()

def write_file(filepath, content):
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

def renumber_sections(content):
    """重新编号所有3.x节"""
    print("重新编号所有3.x节...")
    
    # 定义正确的节顺序和标题
    correct_order = [
        "3.1 Participant Characteristics",
        "3.2 Macular Structural Changes",
        "3.3 Optic Disc Structural Changes",
        "3.4 Correlation with Depression Severity",
        "3.5 Multiple Comparison Correction Summary",
        "3.6 Diagnostic Performance",
        "3.7 Multivariate Analysis",
        "3.8 Machine Learning Optimization of Diagnostic Performance",
        "3.9 Subgroup Analysis by Depression Severity",
        "3.10 Sensitivity Analyses for Age Differences",
        "3.11 Subgroup Analysis by Sex and Age",
        "3.12 Analysis of Refractive Error as a Confounding Factor",
        "3.13 Inter-eye Consistency",
        "3.14 Summary of Key Findings",
        "3.15 Supplementary Analyses"
    ]
    
    # 查找所有现有的3.x节
    lines = content.split('\n')
    section_indices = []
    
    for i, line in enumerate(lines):
        match = re.match(r'^## (3\.\d+)(?!\d)(.*)$', line)
        if match:
            section_indices.append((i, match.group(1), match.group(2).strip()))
    
    print(f"找到 {len(section_indices)} 个3.x节")
    
    # 如果节数量不匹配，可能需要调整
    if len(section_indices) != len(correct_order):
        print(f"⚠ 节数量不匹配: 现有{len(section_indices)}个, 应有{len(correct_order)}个")
    
    # 重新编号现有的节
    for idx, (line_idx, old_num, title) in enumerate(section_indices):
        if idx < len(correct_order):
            new_title = correct_order[idx]
            old_line = lines[line_idx]
            new_line = f"## {new_title}"
            
            # 只更改编号，保留原标题内容
            # 提取原标题内容（去掉编号部分）
            if title:
                # 如果原标题与正确标题不同，使用原标题内容
                # 但使用新的编号
                match = re.match(r'3\.\d+\s+(.*)', f"3.0 {title}")
                if match:
                    title_content = match.group(1)
                    new_line = f"## {correct_order[idx].split(' ', 1)[0]} {title_content}"
            
            lines[line_idx] = new_line
            print(f"  {old_num} → {correct_order[idx].split(' ', 1)[0]}")
    
    # 修复重复的"Inter-eye Consistency"节
    # 查找重复的3.13节
    intereye_count = 0
    for i, line in enumerate(lines):
        if "Inter-eye Consistency" in line and "## 3." in line:
            intereye_count += 1
            if intereye_count > 1:
                # 删除重复的节
                print(f"  删除重复的Inter-eye Consistency节 (第{i+1}行)")
                # 找到节内容结束
                for j in range(i+1, len(lines)):
                    if re.match(r'^## [23]\.', lines[j]):
                        del lines[i:j]
                        break
    
    return '\n'.join(lines)

def main():
    """主函数"""
    paper_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/01_Manuscript/OCT_MDD_Manuscript_463eyes_ML_with_Tables_ROC_Updated_20260315.md"
    
    content = read_file(paper_path)
    
    # 重新编号节
    content = renumber_sections(content)
    
    # 保存
    backup_path = paper_path.replace('.md', '.backup_before_final_section_fix.md')
    write_file(backup_path, read_file(paper_path))
    
    write_file(paper_path, content)
    
    print("\n✅ 节编号修复完成!")
    print(f"备份文件: {backup_path}")

if __name__ == "__main__":
    main()
