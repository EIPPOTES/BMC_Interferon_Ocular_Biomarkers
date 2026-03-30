#!/usr/bin/env python3
"""
检查视力数据完整性
"""
import pandas as pd
import os

# 读取描述性统计详细文件
df = pd.read_excel('/mnt/c/Users/CUI/Desktop/final/03_Tables/Descriptive_Statistics_Detailed_463eyes_20260315.xlsx')

print("="*60)
print("视力数据检查报告")
print("="*60)
print()

print("数据文件列名:")
print(df.columns.tolist())
print()

# 检查是否有视力相关列
vision_keywords = ['vision', 'visual', 'acuity', '视力', 'VA', 'BCVA', 'logMAR']
found_vision = []

for col in df.columns:
    for keyword in vision_keywords:
        if keyword.lower() in str(col).lower():
            found_vision.append(col)
            break

print("="*60)
if found_vision:
    print(f"✅ 发现视力相关数据列: {found_vision}")
else:
    print("❌ 未发现视力相关数据列")
print("="*60)
print()

# 显示所有Category
print("数据类别(Category):")
categories = df['Category'].unique()
for i, cat in enumerate(categories, 1):
    print(f"{i}. {cat}")
