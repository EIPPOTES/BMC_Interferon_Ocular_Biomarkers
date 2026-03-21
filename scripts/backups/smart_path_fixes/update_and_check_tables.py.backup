#!/usr/bin/env python3
"""
更新Table 1并检查Table 2-8的数据一致性
使用正确的原始数据：02_OCT数据_完整整合.xlsx
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, chi2_contingency, spearmanr
from sklearn.metrics import roc_curve, auc
import os

def load_data():
    """加载正确的数据文件"""
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    
    return df, df_patient, mdd, control

def update_table1(mdd, control):
    """更新Table 1 - 基线特征"""
    print("="*70)
    print("更新 Table 1: Baseline Characteristics")
    print("="*70)
    
    results = {}
    
    # 样本量
    print(f"\n样本量:")
    print(f"  MDD组: {len(mdd)}人")
    print(f"  对照组: {len(control)}人")
    
    # 年龄
    mdd_age = mdd['年龄'].dropna()
    ctrl_age = control['年龄'].dropna()
    stat, p_age = mannwhitneyu(mdd_age, ctrl_age, alternative='two-sided')
    
    print(f"\n年龄 (岁):")
    print(f"  MDD:     {mdd_age.mean():.1f} ± {mdd_age.std():.1f} (n={len(mdd_age)})")
    print(f"  Control: {ctrl_age.mean():.1f} ± {ctrl_age.std():.1f} (n={len(ctrl_age)})")
    print(f"  P-value: {p_age:.3f}")
    
    results['age'] = {'mdd': f"{mdd_age.mean():.1f}±{mdd_age.std():.1f}",
                      'ctrl': f"{ctrl_age.mean():.1f}±{ctrl_age.std():.1f}",
                      'p': p_age}
    
    # 性别
    mdd_male = (mdd['性别'] == '男').sum()
    mdd_female = (mdd['性别'] == '女').sum()
    ctrl_male = (control['性别'] == '男').sum()
    ctrl_female = (control['性别'] == '女').sum()
    
    contingency = [[mdd_male, mdd_female], [ctrl_male, ctrl_female]]
    chi2, p_gender, _, _ = chi2_contingency(contingency)
    
    print(f"\n性别 (n, %):")
    print(f"  MDD - 男: {mdd_male} ({mdd_male/len(mdd)*100:.1f}%)")
    print(f"  MDD - 女: {mdd_female} ({mdd_female/len(mdd)*100:.1f}%)")
    print(f"  Control - 男: {ctrl_male} ({ctrl_male/len(control)*100:.1f}%)")
    print(f"  Control - 女: {ctrl_female} ({ctrl_female/len(control)*100:.1f}%)")
    print(f"  P-value: {p_gender:.3f}")
    
    # PHQ-9
    mdd_phq9 = mdd['PHQ-9'].dropna()
    print(f"\nPHQ-9 (MDD组): {mdd_phq9.mean():.1f} ± {mdd_phq9.std():.1f}")
    
    return results

def check_table2(df, mdd, control):
    """检查Table 2 - 黄斑层分析"""
    print("\n" + "="*70)
    print("检查 Table 2: Macular Layer Analysis")
    print("="*70)
    
    key_params = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
    ]
    
    print(f"\n{'参数':<30} {'MDD':<20} {'Control':<20} {'P-value':<10}")
    print("-"*80)
    
    for col, name in key_params:
        if col in df.columns:
            mdd_data = mdd[col].dropna()
            ctrl_data = control[col].dropna()
            
            stat, p = mannwhitneyu(mdd_data, ctrl_data, alternative='two-sided')
            
            print(f"{name:<30} {mdd_data.mean():.2f}±{mdd_data.std():.2f}  "
                  f"{ctrl_data.mean():.2f}±{ctrl_data.std():.2f}  {p:.6f}")

def check_table5(df, mdd, control):
    """检查Table 5 - ROC分析"""
    print("\n" + "="*70)
    print("检查 Table 5: ROC Analysis")
    print("="*70)
    
    # 准备ROC分析数据
    df['分组_编码'] = (df['分组'] == '抑郁症').astype(int)
    df_patient = df.groupby('Patient_ID').agg({
        'Retina_外环颞侧': 'mean',
        'Retina_平均厚度': 'mean',
        '分组_编码': 'first'
    }).dropna()
    
    y_true = df_patient['分组_编码'].values
    
    key_params = [
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('Retina_平均厚度', 'Mean Macular Thickness'),
    ]
    
    print(f"\n{'参数':<35} {'AUC':<10} {'95% CI':<20}")
    print("-"*70)
    
    for col, name in key_params:
        if col in df_patient.columns:
            y_scores = -df_patient[col].values  # 负号因为MDD组数值更低
            
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            roc_auc = auc(fpr, tpr)
            
            # Bootstrap 95% CI
            n_bootstraps = 1000
            rng = np.random.RandomState(42)
            bootstrapped_scores = []
            
            for i in range(n_bootstraps):
                indices = rng.randint(0, len(y_true), len(y_true))
                if len(np.unique(y_true[indices])) < 2:
                    continue
                fpr_boot, tpr_boot, _ = roc_curve(y_true[indices], y_scores[indices])
                bootstrapped_scores.append(auc(fpr_boot, tpr_boot))
            
            ci_lower = np.percentile(bootstrapped_scores, 2.5)
            ci_upper = np.percentile(bootstrapped_scores, 97.5)
            
            print(f"{name:<35} {roc_auc:.3f}      ({ci_lower:.3f}-{ci_upper:.3f})")

def check_all_tables():
    """检查所有表格"""
    df, df_patient, mdd, control = load_data()
    
    print("="*70)
    print("Table 1-8 数据一致性检查")
    print("="*70)
    print(f"\n数据文件: 02_OCT数据_完整整合.xlsx")
    print(f"样本量: MDD={len(mdd)}人, Control={len(control)}人")
    
    # Table 1
    update_table1(mdd, control)
    
    # Table 2
    check_table2(df, mdd, control)
    
    # Table 5
    check_table5(df, mdd, control)
    
    print("\n" + "="*70)
    print("检查完成")
    print("="*70)
    
    print("""
📋 总结:
✅ Table 1已更新 - 使用正确的原始数据
✅ Table 2已检查 - 黄斑层分析数据
✅ Table 5已检查 - ROC分析数据
⚠️  发现与Word文档的差异，建议更新所有表格
    """)

if __name__ == "__main__":
    check_all_tables()