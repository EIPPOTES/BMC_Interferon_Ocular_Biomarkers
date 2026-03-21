#!/usr/bin/env python3
# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
完整检查Table 3-4,6-8的数据一致性
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
import os

def load_data():
    """加载数据"""
    data_file = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    df = pd.read_excel(data_file)
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    return df, df_patient, mdd, control

def check_table3(df, mdd, control):
    """Table 3 - 视盘参数"""
    print("\n" + "="*70)
    print("Table 3: Optic Disc Parameters (视盘参数)")
    print("="*70)
    
    optic_params = [
        ('RNFL_Total', 'Total RNFL Thickness'),
        ('Cup Area', 'Cup Area'),
        ('Rim Area', 'Rim Area'),
        ('C/D Area Ratio', 'Cup-to-Disc Ratio'),
    ]
    
    print(f"\n{'参数':<30} {'MDD':<20} {'Control':<20} {'P-value':<10}")
    print("-"*80)
    
    for col, name in optic_params:
        if col in df.columns:
            mdd_data = mdd[col].dropna()
            ctrl_data = control[col].dropna()
            stat, p = mannwhitneyu(mdd_data, ctrl_data, alternative='two-sided')
            
            print(f"{name:<30} {mdd_data.mean():.2f}±{mdd_data.std():.2f}  "
                  f"{ctrl_data.mean():.2f}±{ctrl_data.std():.2f}  {p:.6f}")

def check_table4(df, mdd, control):
    """Table 4 - 相关性分析"""
    print("\n" + "="*70)
    print("Table 4: Correlation Analysis (PHQ-9与OCT参数相关性)")
    print("="*70)
    
    # 仅MDD组
    mdd_with_phq9 = mdd[mdd['PHQ-9'].notna()]
    
    corr_params = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('RNFL_Total', 'Total RNFL Thickness'),
    ]
    
    print(f"\n{'参数':<35} {'r':<10} {'P-value':<10} {'n':<10}")
    print("-"*70)
    
    for col, name in corr_params:
        if col in mdd_with_phq9.columns:
            x = mdd_with_phq9['PHQ-9']
            y = mdd_with_phq9[col]
            
            # 移除缺失值
            mask = x.notna() & y.notna()
            x_clean = x[mask]
            y_clean = y[mask]
            
            if len(x_clean) > 2:
                r, p = spearmanr(x_clean, y_clean)
                print(f"{name:<35} {r:>8.3f}  {p:>8.3f}  {len(x_clean):>8}")

def check_table6(df, mdd, control):
    """Table 6 - 多变量回归"""
    print("\n" + "="*70)
    print("Table 6: Multivariate Regression (多变量回归)")
    print("="*70)
    
    # 简单的线性回归分析
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler
    
    # 准备数据
    df_reg = df_patient[['年龄', '性别', 'Retina_平均厚度', '分组']].dropna()
    df_reg['性别_编码'] = (df_reg['性别'] == '男').astype(int)
    df_reg['分组_编码'] = (df_reg['分组'] == '抑郁症').astype(int)
    
    X = df_reg[['年龄', '性别_编码', '分组_编码']]
    y = df_reg['Retina_平均厚度']
    
    model = LinearRegression()
    model.fit(X, y)
    
    print(f"\n因变量: Mean Macular Thickness")
    print(f"\n{'自变量':<20} {'β系数':<15} {'解释':<30}")
    print("-"*70)
    
    variables = ['Age', 'Male Gender', 'MDD Status']
    for var, coef in zip(variables, model.coef_):
        direction = "增加" if coef < 0 else "减少"
        print(f"{var:<20} {coef:>10.3f}     MDD与黄斑厚度{direction}相关")
    
    print(f"\nR² = {model.score(X, y):.3f}")

def check_table7(mdd):
    """Table 7 - 亚组分析"""
    print("\n" + "="*70)
    print("Table 7: Subgroup Analysis by PHQ-9 Severity")
    print("="*70)
    
    # 按PHQ-9分层
    def classify_phq9(score):
        if pd.isna(score):
            return 'Unknown'
        elif score < 5:
            return 'Minimal (0-4)'
        elif score < 10:
            return 'Mild (5-9)'
        elif score < 15:
            return 'Moderate (10-14)'
        else:
            return 'Severe (≥15)'
    
    mdd_analysis = mdd.copy()
    mdd_analysis['PHQ9_Group'] = mdd_analysis['PHQ-9'].apply(classify_phq9)
    
    print(f"\n{'PHQ-9分组':<20} {'n':<10} {'Mean Macular Thickness':<25}")
    print("-"*60)
    
    for group in ['Minimal (0-4)', 'Mild (5-9)', 'Moderate (10-14)', 'Severe (≥15)']:
        subset = mdd_analysis[mdd_analysis['PHQ9_Group'] == group]
        if len(subset) > 0:
            thickness = subset['Retina_平均厚度'].dropna()
            print(f"{group:<20} {len(subset):<10} {thickness.mean():.2f}±{thickness.std():.2f}")

def check_table8(df, mdd, control):
    """Table 8 - 双眼一致性"""
    print("\n" + "="*70)
    print("Table 8: Inter-eye Consistency (双眼一致性)")
    print("="*70)
    
    # 找出有双眼数据的参与者
    eye_counts = df.groupby('Patient_ID')['Eye'].nunique()
    bilateral_patients = eye_counts[eye_counts == 2].index
    
    df_bilateral = df[df['Patient_ID'].isin(bilateral_patients)]
    
    # 计算左右眼相关性
    results = []
    for patient in bilateral_patients[:20]:  # 取前20个作为示例
        patient_data = df_bilateral[df_bilateral['Patient_ID'] == patient]
        if len(patient_data) == 2:
            left = patient_data[patient_data['Eye'] == 'L']['Retina_平均厚度'].values
            right = patient_data[patient_data['Eye'] == 'R']['Retina_平均厚度'].values
            if len(left) > 0 and len(right) > 0:
                results.append((left[0], right[0]))
    
    if results:
        left_vals = [r[0] for r in results]
        right_vals = [r[1] for r in results]
        r, p = spearmanr(left_vals, right_vals)
        
        print(f"\n双眼黄斑平均厚度相关性:")
        print(f"  Pearson r = {r:.3f}")
        print(f"  P-value = {p:.3f}")
        print(f"  n = {len(results)} participants")
        print(f"  Mean absolute difference = {np.mean(np.abs(np.array(left_vals) - np.array(right_vals))):.2f} μm")

def main():
    print("="*70)
    print("完整检查 Table 3-4,6-8 数据一致性")
    print("="*70)
    
    global df_patient  # 用于table6
    df, df_patient, mdd, control = load_data()
    
    print(f"\n数据文件: 02_OCT数据_完整整合.xlsx")
    print(f"样本量: MDD={len(mdd)}人, Control={len(control)}人")
    
    # 检查各表格
    check_table3(df, mdd, control)
    check_table4(df, mdd, control)
    check_table6(df, mdd, control)
    check_table7(mdd)
    check_table8(df, mdd, control)
    
    print("\n" + "="*70)
    print("✅ Table 3-4,6-8 检查完成")
    print("="*70)

if __name__ == "__main__":
    main()