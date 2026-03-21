#!/usr/bin/env python3
"""
分析data.xlsx文件结构，为修复Figure 5做准备
"""

import pandas as pd
import numpy as np
import os

def analyze_data_file():
    data_file = '../data/raw/data.xlsx'
    
    if not os.path.exists(data_file):
        print(f"错误: 数据文件不存在: {data_file}")
        return False
    
    try:
        print(f"读取数据文件: {data_file}")
        df = pd.read_excel(data_file)
        
        print(f"\n数据基本信息:")
        print(f"  行数: {len(df)}")
        print(f"  列数: {len(df.columns)}")
        
        print(f"\n列名列表:")
        for i, col in enumerate(df.columns):
            print(f"  {i+1:2d}. {col}")
        
        print(f"\n数据前几行:")
        print(df.head())
        
        print(f"\n数据类型:")
        print(df.dtypes.head(20))
        
        # 查找分组列
        group_cols = [col for col in df.columns if 'group' in col.lower() or '分组' in col]
        print(f"\n分组相关列: {group_cols}")
        
        # 查找OCT相关列
        oct_cols = [col for col in df.columns if 'OCT' in col or 'retina' in col.lower() or 
                   'rnfl' in col.lower() or 'gcl' in col.lower() or 'thickness' in col.lower()]
        print(f"\nOCT相关列 (前20个):")
        for col in oct_cols[:20]:
            print(f"  {col}")
        
        # 如果有分组列，显示分布
        if group_cols:
            group_col = group_cols[0]
            print(f"\n分组分布:")
            print(df[group_col].value_counts())
            
            # 显示每组样本量
            print(f"\n每组样本量:")
            for group in df[group_col].unique():
                n = len(df[df[group_col] == group])
                print(f"  {group}: {n}行")
        
        return True
        
    except Exception as e:
        print(f"读取数据时出错: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    analyze_data_file()