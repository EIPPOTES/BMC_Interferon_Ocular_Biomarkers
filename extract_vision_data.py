#!/usr/bin/env python3
"""
从原始数据提取视力基线数据
"""
import pandas as pd
import numpy as np

# 读取原始数据
df = pd.read_excel('/mnt/c/Users/CUI/Desktop/final/04_Data/data_499eyes_20260315.xlsx')

print("="*60)
print("原始数据列名检查")
print("="*60)
print(f"总列数: {len(df.columns)}")
print()

# 查找视力相关列
vision_keywords = ['vision', 'visual', 'acuity', '视力', 'VA', 'BCVA', 'logMAR', 'LogMAR', 'logmar']
found_vision_cols = []

for col in df.columns:
    for keyword in vision_keywords:
        if keyword.lower() in str(col).lower():
            found_vision_cols.append(col)
            break

print("发现的视力相关列:")
if found_vision_cols:
    for col in found_vision_cols:
        print(f"  - {col}")
else:
    print("  未发现视力相关列，显示所有列名:")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i}. {col}")

print()
print("="*60)

# 如果有视力列，提取统计信息
if found_vision_cols:
    for vision_col in found_vision_cols:
        print(f"\n视力列 '{vision_col}' 统计:")
        print(f"  数据类型: {df[vision_col].dtype}")
        print(f"  非空值数: {df[vision_col].notna().sum()}")
        print(f"  空值数: {df[vision_col].isna().sum()}")
        print(f"  描述统计:")
        print(df[vision_col].describe())
        
        # 按组别统计
        if 'Group' in df.columns or 'group' in df.columns:
            group_col = 'Group' if 'Group' in df.columns else 'group'
            print(f"\n  按组别统计 ({group_col}):")
            for group in df[group_col].unique():
                group_data = df[df[group_col] == group][vision_col].dropna()
                print(f"    {group}: n={len(group_data)}, Mean={group_data.mean():.2f}, SD={group_data.std():.2f}")
