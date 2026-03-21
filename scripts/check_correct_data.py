#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
检查正确的原始数据文件：02_OCT数据_完整整合.xlsx
"""

import pandas as pd
import os

def check_correct_data_file():
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    
    print("="*70)
    print("检查正确的原始数据文件")
    print("="*70)
    
    if not os.path.exists(data_file):
        print(f"❌ 文件不存在: {data_file}")
        return
    
    print(f"\n📁 文件: {data_file}")
    print(f"   大小: {os.path.getsize(data_file)/1024:.1f} KB")
    
    # 读取所有工作表
    xl = pd.ExcelFile(data_file)
    sheet_names = xl.sheet_names
    
    print(f"\n📊 工作表数量: {len(sheet_names)}")
    print(f"   工作表名称: {sheet_names}")
    
    # 检查每个工作表
    for sheet in sheet_names:
        print(f"\n{'='*70}")
        print(f"工作表: {sheet}")
        print(f"{'='*70}")
        
        df = pd.read_excel(data_file, sheet_name=sheet)
        
        print(f"   形状: {df.shape} (行×列)")
        print(f"   列名: {list(df.columns[:10])}{'...' if len(df.columns) > 10 else ''}")
        
        # 检查分组列
        if '分组' in df.columns:
            print(f"\n   分组分布:")
            print(df['分组'].value_counts())
        
        # 检查年龄列
        if '年龄' in df.columns:
            age = df['年龄'].dropna()
            print(f"\n   年龄统计:")
            print(f"     均值: {age.mean():.1f} ± {age.std():.1f}岁")
            print(f"     范围: {age.min():.0f} - {age.max():.0f}岁")
            print(f"     样本量: {len(age)}")
        
        # 显示前几行
        print(f"\n   前3行数据:")
        print(df.head(3))

if __name__ == "__main__":
    check_correct_data_file()