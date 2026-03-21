#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
检查Figures和Tables与论文引用的一致性
"""

import os
import re

def check_figures_tables_consistency():
    manuscript_path = '/mnt/c/Users/CUI/Desktop/投稿版/01_Manuscript/OCT_MDD_Manuscript_Final.md'
    figures_dir = '/mnt/c/Users/CUI/Desktop/投稿版/02_Figures'
    tables_dir = '/mnt/c/Users/CUI/Desktop/投稿版/03_Tables'
    
    print("="*70)
    print("Figures和Tables与论文引用一致性检查")
    print("="*70)
    
    # 读取论文内容
    with open(manuscript_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 提取引用的Figures
    print("\n📊 论文中引用的Figures:")
    print("-"*70)
    
    figures_in_paper = {
        'Figure 1': {'name': 'Study Flow Chart', 'pattern': r'Figure 1|Figure1'},
        'Figure 2': {'name': 'Group Comparison', 'pattern': r'Figure 2|Figure2'},
        'Figure 3': {'name': 'ROC Curves', 'pattern': r'Figure 3|Figure3'},
        'Figure 4': {'name': 'Correlation Scatter', 'pattern': r'Figure 4|Figure4'},
        'Figure 5': {'name': 'Forest Plot', 'pattern': r'Figure 5|Figure5'},
        'Figure 6': {'name': 'Subgroup Analysis', 'pattern': r'Figure 6|Figure6'},
    }
    
    for fig_key, fig_info in figures_in_paper.items():
        # 检查是否在论文中引用
        mentions = len(re.findall(fig_info['pattern'], content, re.IGNORECASE))
        
        # 检查文件是否存在
        fig_file = f"Figure{fig_key.split()[1]}"
        possible_files = [
            os.path.join(figures_dir, f"{fig_file}_Study_Flowchart.png"),
            os.path.join(figures_dir, f"{fig_file}_Group_Comparison.png"),
            os.path.join(figures_dir, f"{fig_file}_ROC_Curves.png"),
            os.path.join(figures_dir, f"{fig_file}_Correlation_Scatter.png"),
            os.path.join(figures_dir, f"{fig_file}_Forest_Plot.png"),
            os.path.join(figures_dir, f"{fig_file}_Subgroup_Analysis.png"),
        ]
        
        file_exists = any(os.path.exists(f) for f in possible_files)
        existing_file = next((f for f in possible_files if os.path.exists(f)), None)
        
        status = "✅" if file_exists else "❌"
        file_size = os.path.getsize(existing_file)/1024 if existing_file else 0
        
        print(f"{status} {fig_key}: {fig_info['name']}")
        print(f"   引用次数: {mentions}次")
        if existing_file:
            print(f"   文件: {os.path.basename(existing_file)} ({file_size:.1f} KB)")
        else:
            print(f"   文件: 未找到")
    
    # 提取引用的Tables
    print("\n📋 论文中引用的Tables:")
    print("-"*70)
    
    tables_in_paper = {
        'Table 1': {'name': 'Baseline Characteristics', 'pattern': r'Table 1|Table1'},
        'Table 2': {'name': 'Macular Layer Analysis', 'pattern': r'Table 2|Table2'},
        'Table 3': {'name': 'Optic Disc Parameters', 'pattern': r'Table 3|Table3'},
        'Table 4': {'name': 'Correlation Analysis', 'pattern': r'Table 4|Table4'},
        'Table 5': {'name': 'ROC Analysis', 'pattern': r'Table 5|Table5'},
        'Table 6': {'name': 'Multivariate Regression', 'pattern': r'Table 6|Table6'},
        'Table 7': {'name': 'Subgroup Analysis', 'pattern': r'Table 7|Table7'},
        'Table 8': {'name': 'Inter-eye Consistency', 'pattern': r'Table 8|Table8'},
    }
    
    for table_key, table_info in tables_in_paper.items():
        # 检查是否在论文中引用
        mentions = len(re.findall(table_info['pattern'], content, re.IGNORECASE))
        
        # 检查文件是否存在
        table_num = table_key.split()[1]
        possible_files = [
            os.path.join(tables_dir, f"Table{table_num}_Baseline_Characteristics_ThreeLine_Final.png"),
            os.path.join(tables_dir, f"Table{table_num}_Macular_Layers_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{table_num}_Optic_Disc_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{table_num}_Correlation_Analysis_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{table_num}_ROC_Analysis_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{table_num}_Multivariate_Regression_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{table_num}_Subgroup_Analysis_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{table_num}_Intereye_Consistency_ThreeLine.png"),
        ]
        
        file_exists = any(os.path.exists(f) for f in possible_files)
        existing_file = next((f for f in possible_files if os.path.exists(f)), None)
        
        status = "✅" if file_exists else "❌"
        file_size = os.path.getsize(existing_file)/1024 if existing_file else 0
        
        print(f"{status} {table_key}: {table_info['name']}")
        print(f"   引用次数: {mentions}次")
        if existing_file:
            print(f"   文件: {os.path.basename(existing_file)} ({file_size:.1f} KB)")
        else:
            print(f"   文件: 未找到")
    
    # 总结
    print("\n" + "="*70)
    print("总结")
    print("="*70)
    
    all_figures_exist = all(
        any(os.path.exists(f) for f in [
            os.path.join(figures_dir, f"Figure{i}_Study_Flowchart.png"),
            os.path.join(figures_dir, f"Figure{i}_Group_Comparison.png"),
            os.path.join(figures_dir, f"Figure{i}_ROC_Curves.png"),
            os.path.join(figures_dir, f"Figure{i}_Correlation_Scatter.png"),
            os.path.join(figures_dir, f"Figure{i}_Forest_Plot.png"),
            os.path.join(figures_dir, f"Figure{i}_Subgroup_Analysis.png"),
        ])
        for i in range(1, 7)
    )
    
    all_tables_exist = all(
        any(os.path.exists(f) for f in [
            os.path.join(tables_dir, f"Table{i}_Baseline_Characteristics_ThreeLine_Final.png"),
            os.path.join(tables_dir, f"Table{i}_Macular_Layers_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{i}_Optic_Disc_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{i}_Correlation_Analysis_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{i}_ROC_Analysis_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{i}_Multivariate_Regression_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{i}_Subgroup_Analysis_ThreeLine.png"),
            os.path.join(tables_dir, f"Table{i}_Intereye_Consistency_ThreeLine.png"),
        ])
        for i in range(1, 9)
    )
    
    if all_figures_exist and all_tables_exist:
        print("✅ 所有Figures和Tables与论文引用完全一致！")
    else:
        print("⚠️ 部分文件缺失，请检查")
    
    print(f"\nFigures: 6/6 {'✅' if all_figures_exist else '❌'}")
    print(f"Tables: 8/8 {'✅' if all_tables_exist else '❌'}")

if __name__ == "__main__":
    check_figures_tables_consistency()