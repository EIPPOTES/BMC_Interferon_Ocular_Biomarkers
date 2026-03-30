#!/usr/bin/env python3
"""
查找最新版统计分析结果
"""

import os
import pandas as pd
from datetime import datetime
import re

# 目录路径
base_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改"
tables_dir = os.path.join(base_dir, "03_Tables")
unified_dir = os.path.join(base_dir, "08_Unified_Analysis_20260315", "unified_analysis_results_20260315")

print("=" * 80)
print("查找最新版统计分析结果")
print(f"搜索目录: {tables_dir}")
print("=" * 80)

# 定义分析类型和关键词映射
analysis_types = {
    "描述性统计分析": ["描述性", "Descriptive", "descriptive"],
    "组间比较": ["组间比较", "Group_Comparison", "Group Comparison", "MDD vs Control"],
    "相关性分析": ["相关性", "Correlation", "correlation"],
    "ROC分析": ["ROC", "AUC"],
    "多变量回归分析": ["多变量回归", "Multivariate", "multivariate", "回归分析"],
    "亚组分析": ["亚组分析", "Subgroup", "subgroup"]
}

# 收集所有Excel文件
excel_files = []
for root, dirs, files in os.walk(tables_dir):
    for file in files:
        if file.endswith('.xlsx'):
            excel_files.append(os.path.join(root, file))

print(f"发现 {len(excel_files)} 个Excel文件")

# 按分析类型分类文件
categorized_files = {atype: [] for atype in analysis_types.keys()}
categorized_files["其他"] = []

for file_path in excel_files:
    filename = os.path.basename(file_path)
    matched = False
    
    for atype, keywords in analysis_types.items():
        for keyword in keywords:
            if keyword in filename:
                # 尝试获取文件修改时间
                mtime = os.path.getmtime(file_path)
                mtime_str = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M')
                
                # 获取文件大小
                size_kb = os.path.getsize(file_path) / 1024
                
                categorized_files[atype].append({
                    'file': filename,
                    'path': file_path,
                    'mtime': mtime,
                    'mtime_str': mtime_str,
                    'size_kb': size_kb,
                    'full_path': file_path
                })
                matched = True
                break
        if matched:
            break
    
    if not matched:
        mtime = os.path.getmtime(file_path)
        mtime_str = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M')
        size_kb = os.path.getsize(file_path) / 1024
        categorized_files["其他"].append({
            'file': filename,
            'path': file_path,
            'mtime': mtime,
            'mtime_str': mtime_str,
            'size_kb': size_kb,
            'full_path': file_path
        })

# 检查统一分析结果
unified_file = os.path.join(unified_dir, "统一分析结果汇总.xlsx")
if os.path.exists(unified_file):
    mtime = os.path.getmtime(unified_file)
    mtime_str = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M')
    size_kb = os.path.getsize(unified_file) / 1024
    
    # 统一分析包含多个分析类型
    unified_entry = {
        'file': "统一分析结果汇总.xlsx",
        'path': unified_file,
        'mtime': mtime,
        'mtime_str': mtime_str,
        'size_kb': size_kb,
        'full_path': unified_file,
        'note': "包含多变量回归、亚组分析、ROC分析、相关性分析"
    }
    
    # 添加到相应类别
    for atype in ["多变量回归分析", "亚组分析", "ROC分析", "相关性分析"]:
        if atype in categorized_files:
            categorized_files[atype].append(unified_entry.copy())

# 找出每个类别的最新文件
latest_files = {}
for atype, files in categorized_files.items():
    if files:
        # 按修改时间排序
        sorted_files = sorted(files, key=lambda x: x['mtime'], reverse=True)
        latest = sorted_files[0]
        
        # 检查是否有多个最新文件
        latest_mtime = latest['mtime']
        same_time_files = [f for f in sorted_files if abs(f['mtime'] - latest_mtime) < 60]  # 60秒内视为同时
        
        latest_files[atype] = {
            'latest': latest,
            'all_latest': same_time_files,
            'total_count': len(files)
        }

# 打印结果
print("\n" + "=" * 80)
print("最新版统计分析结果汇总")
print("=" * 80)

for atype in analysis_types.keys():
    if atype in latest_files:
        info = latest_files[atype]
        latest = info['latest']
        print(f"\n📊 {atype}:")
        print(f"   最新文件: {latest['file']}")
        print(f"   修改时间: {latest['mtime_str']}")
        print(f"   文件大小: {latest['size_kb']:.1f} KB")
        print(f"   文件位置: {latest['path']}")
        
        if 'note' in latest:
            print(f"   备注: {latest['note']}")
        
        # 如果有多个文件同时最新
        if len(info['all_latest']) > 1:
            print(f"   同时最新的文件还有:")
            for f in info['all_latest'][1:]:
                print(f"     • {f['file']} ({f['mtime_str']})")
        
        print(f"   该类别共有 {info['total_count']} 个文件")
        
        # 显示文件内容预览
        try:
            if os.path.exists(latest['full_path']):
                # 检查文件是否可读
                try:
                    xl = pd.ExcelFile(latest['full_path'])
                    sheets = xl.sheet_names
                    print(f"   包含工作表: {sheets}")
                    
                    # 读取第一个工作表的前几行
                    if sheets:
                        df = xl.parse(sheets[0])
                        print(f"   数据形状: {df.shape[0]} 行 × {df.shape[1]} 列")
                        
                        if df.shape[0] > 0 and df.shape[1] <= 10:
                            print(f"   前3行数据:")
                            print(df.head(3).to_string(index=False))
                        elif df.shape[1] > 10:
                            print(f"   列名示例: {list(df.columns)[:5]}...")
                except Exception as e:
                    print(f"   读取文件内容时出错: {str(e)[:50]}")
        except Exception as e:
            print(f"   访问文件时出错: {str(e)[:50]}")

# 打印未分类的文件
if categorized_files["其他"]:
    print(f"\n📁 未分类文件 ({len(categorized_files['其他'])} 个):")
    for file_info in categorized_files["其他"][:10]:  # 只显示前10个
        print(f"   • {file_info['file']} ({file_info['mtime_str']}, {file_info['size_kb']:.1f} KB)")

print("\n" + "=" * 80)
print("分析完成")
print("=" * 80)