#!/usr/bin/env python3
"""
检查论文中引用的表格哪些已经生成，哪些缺失
"""

import os
import re

def extract_table_references(paper_path):
    """从论文中提取表格引用"""
    with open(paper_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 查找表格引用
    table_pattern = r'Table\s+([0-9S]+)'
    matches = re.findall(table_pattern, content, re.IGNORECASE)
    
    # 去重并排序
    tables = []
    for match in matches:
        # 处理S1, S2等补充表格
        if match.startswith('S'):
            table_ref = f"Supplementary Table {match}"
        else:
            table_ref = f"Table {match}"
        
        if table_ref not in tables:
            tables.append(table_ref)
    
    # 按数字排序（主表格在前，补充表格在后）
    def sort_key(ref):
        if 'Supplementary' in ref:
            num = ref.replace('Supplementary Table S', '').strip()
            return (1, int(num) if num.isdigit() else 999)
        else:
            num = ref.replace('Table ', '').strip()
            return (0, int(num) if num.isdigit() else 999)
    
    tables.sort(key=sort_key)
    return tables

def check_existing_tables(tables_dir, manuscript_dir):
    """检查现有表格文件"""
    existing_files = {}
    
    # 检查03_Tables目录
    for root, dirs, files in os.walk(tables_dir):
        for file in files:
            if file.endswith('.xlsx'):
                existing_files[file] = os.path.join(root, file)
    
    # 检查01_Manuscript目录
    for root, dirs, files in os.walk(manuscript_dir):
        for file in files:
            if file.endswith('.xlsx'):
                existing_files[file] = os.path.join(root, file)
    
    return existing_files

def map_table_to_files(table_ref, existing_files):
    """将表格引用映射到可能的文件"""
    table_num = table_ref.split()[-1]
    
    # 可能的文件名模式
    possible_patterns = []
    
    if 'Supplementary' in table_ref:
        # 补充表格
        possible_patterns = [
            f"Supplementary_Table_{table_num}.xlsx",
            f"TableS{table_num}.xlsx",
            f"补充表格_{table_num}.xlsx",
            f"Supplementary_{table_num}.xlsx"
        ]
    else:
        # 主表格
        possible_patterns = [
            f"Table{table_num}_*.xlsx",
            f"Table{table_num}.xlsx",
            f"Table_{table_num}_*.xlsx",
            f"Table_{table_num}.xlsx"
        ]
    
    # 尝试匹配
    matched_files = []
    for pattern in possible_patterns:
        # 简化模式匹配
        pattern_lower = pattern.lower()
        for file in existing_files.keys():
            file_lower = file.lower()
            if pattern_lower.replace('*', '') in file_lower:
                if file not in matched_files:
                    matched_files.append(file)
    
    return matched_files

def main():
    """主函数"""
    # 路径
    paper_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/01_Manuscript/OCT_MDD_Manuscript_463eyes_ML_with_Tables_ROC_Updated_20260315.md"
    tables_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/03_Tables"
    manuscript_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/01_Manuscript"
    
    print("检查论文表格完整性...")
    print("=" * 60)
    
    # 提取表格引用
    table_refs = extract_table_references(paper_path)
    print(f"论文中引用的表格: {len(table_refs)} 个")
    for i, ref in enumerate(table_refs, 1):
        print(f"  {i}. {ref}")
    
    # 检查现有文件
    existing_files = check_existing_tables(tables_dir, manuscript_dir)
    print(f"\n现有表格文件: {len(existing_files)} 个")
    for file in sorted(existing_files.keys()):
        print(f"  {file}")
    
    # 映射表格到文件
    print(f"\n表格匹配情况:")
    print("-" * 60)
    
    missing_tables = []
    for table_ref in table_refs:
        matched_files = map_table_to_files(table_ref, existing_files)
        
        if matched_files:
            print(f"✓ {table_ref}:")
            for file in matched_files:
                print(f"    → {file}")
        else:
            print(f"✗ {table_ref}: 未找到对应文件")
            missing_tables.append(table_ref)
    
    # 分析缺失的表格
    print(f"\n缺失表格分析:")
    print("-" * 60)
    
    if missing_tables:
        print(f"总共缺失 {len(missing_tables)} 个表格:")
        for table_ref in missing_tables:
            table_num = table_ref.split()[-1]
            print(f"  {table_ref}")
            
            # 根据表格编号推测内容
            if table_ref == "Table 4":
                print("    推测内容: Correlation analysis between PHQ-9 scores and OCT parameters")
                print("    可能对应文件: 相关性分析_OCT_vs_PHQ9_20260315.xlsx")
            elif table_ref == "Table 6":
                print("    推测内容: Multiple linear regression analysis")
                print("    可能对应文件: 多变量回归_线性模型结果_20260315.xlsx 或 多变量回归_混合效应模型结果_20260315.xlsx")
            elif table_ref == "Table 7":
                print("    推测内容: MDD patients stratified by PHQ-9 scores")
                print("    可能对应文件: 未找到明确的PHQ-9分层分析文件")
            elif "Supplementary" in table_ref:
                print(f"    推测内容: 补充表格 {table_num}")
    else:
        print("✓ 所有表格都有对应文件")
    
    # 检查表格文件是否在论文目录中（便于整合）
    print(f"\n表格文件位置检查:")
    print("-" * 60)
    
    table_files_in_manuscript = []
    for file, path in existing_files.items():
        if '01_Manuscript' in path:
            table_files_in_manuscript.append(file)
    
    if table_files_in_manuscript:
        print(f"在论文目录中的表格文件 ({len(table_files_in_manuscript)} 个):")
        for file in table_files_in_manuscript:
            print(f"  ✓ {file}")
    else:
        print("没有表格文件在论文目录中")
    
    # 建议
    print(f"\n建议:")
    print("-" * 60)
    
    if missing_tables:
        print("1. 为以下缺失表格创建独立的Excel文件:")
        for table_ref in missing_tables:
            print(f"   - {table_ref}")
        
        print("\n2. 将所有表格文件复制到论文目录中便于投稿:")
        print("   mv 03_Tables/*Table*.xlsx 01_Manuscript/")
        
        print("\n3. 确保论文中的表格引用与文件命名一致")
    else:
        print("✓ 所有表格都已生成")
        print("\n1. 确保表格文件已正确嵌入论文Word文档")
        print("2. 检查表格格式是否符合期刊要求")
        print("3. 验证表格数据与论文正文的一致性")

if __name__ == "__main__":
    main()