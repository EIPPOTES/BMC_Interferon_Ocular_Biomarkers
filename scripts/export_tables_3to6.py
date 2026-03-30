#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
导出Table 3-6的原始分析数据
"""

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import os

OUTPUT_DIR = '/mnt/c/Users/CUI/Desktop/投稿、数据修改/05_Raw_Data'

def load_data():
    """加载485眼数据"""
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx')
    control_ages = df[df['分组'] == '健康对照'].drop_duplicates('Patient_ID')[['Patient_ID', '年龄']]
    age_missing_ids = control_ages[control_ages['年龄'].isna()]['Patient_ID'].tolist()
    df_filtered = df[~((df['分组'] == '健康对照') & (df['Patient_ID'].isin(age_missing_ids)))]
    
    # 创建分组变量 (0=Control, 1=MDD)
    df_filtered['Group_Binary'] = (df_filtered['分组'] == '抑郁症').astype(int)
    
    return df_filtered

def export_table3_data(df):
    """导出Table 3 - 多变量回归分析"""
    
    # 准备数据
    df_reg = df[df['Retina_平均厚度'].notna()].copy()
    
    # 多变量回归
    X = df_reg[['Group_Binary', '年龄', '性别']].copy()
    X['性别'] = (X['性别'] == '女').astype(int)
    X = sm.add_constant(X)
    y = df_reg['Retina_平均厚度']
    
    model = sm.OLS(y, X).fit()
    
    # 提取结果
    results = []
    for var in model.params.index:
        results.append({
            'Variable': var,
            'Coefficient': model.params[var],
            'Std_Error': model.bse[var],
            't_value': model.tvalues[var],
            'P_value': model.pvalues[var],
            'CI_Lower': model.conf_int()[0][var],
            'CI_Upper': model.conf_int()[1][var],
        })
    
    # 模型拟合度
    results.append({'Variable': 'R_squared', 'Coefficient': model.rsquared, 'Std_Error': np.nan, 't_value': np.nan, 'P_value': np.nan, 'CI_Lower': np.nan, 'CI_Upper': np.nan})
    results.append({'Variable': 'Adj_R_squared', 'Coefficient': model.rsquared_adj, 'Std_Error': np.nan, 't_value': np.nan, 'P_value': np.nan, 'CI_Lower': np.nan, 'CI_Upper': np.nan})
    results.append({'Variable': 'F_statistic', 'Coefficient': model.fvalue, 'Std_Error': np.nan, 't_value': np.nan, 'P_value': model.f_pvalue, 'CI_Lower': np.nan, 'CI_Upper': np.nan})
    
    df_result = pd.DataFrame(results)
    df_result.to_excel(f'{OUTPUT_DIR}/Table3_Multivariate_Regression_RawData.xlsx', index=False)
    print("✅ Table 3原始数据已导出 (多变量回归)")
    return df_result

def export_table4_data(df):
    """导出Table 4 - ROC分析详细数据"""
    
    from sklearn.metrics import roc_curve, auc
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('RNFL_上方', 'Superior RNFL'),
    ]
    
    results = []
    for col, name in indicators:
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        y_true = [1] * len(mdd_data) + [0] * len(control_data)
        y_scores = list(mdd_data) + list(control_data)
        
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        # 计算不同截断点的敏感度和特异度
        for i in range(0, len(thresholds), max(1, len(thresholds)//10)):
            results.append({
                'Parameter': name,
                'Threshold': thresholds[i],
                'Sensitivity': tpr[i],
                'Specificity': 1 - fpr[i],
                'FPR': fpr[i],
                'AUC': roc_auc,
            })
    
    df_result = pd.DataFrame(results)
    df_result.to_excel(f'{OUTPUT_DIR}/Table4_ROC_Detailed_RawData.xlsx', index=False)
    print("✅ Table 4原始数据已导出 (ROC详细)")
    return df_result

def export_table5_data(df):
    """导出Table 5 - 相关性分析矩阵"""
    
    mdd_df = df[df['分组'] == '抑郁症'].copy()
    
    indicators = [
        'Retina_平均厚度',
        'Retina_外环颞侧',
        'Retina_内环颞侧',
        'RNFL_上方',
        'PHQ-9',
    ]
    
    # 创建相关矩阵
    corr_data = mdd_df[indicators].corr()
    
    # 添加P值矩阵
    pvalue_matrix = pd.DataFrame(np.zeros((len(indicators), len(indicators))), 
                                 columns=indicators, index=indicators)
    
    for i, col1 in enumerate(indicators):
        for j, col2 in enumerate(indicators):
            if i != j:
                valid_data = mdd_df[[col1, col2]].dropna()
                if len(valid_data) > 2:
                    _, p = stats.pearsonr(valid_data[col1], valid_data[col2])
                    pvalue_matrix.loc[col1, col2] = p
    
    # 保存为两个sheet
    with pd.ExcelWriter(f'{OUTPUT_DIR}/Table5_Correlation_Matrix_RawData.xlsx') as writer:
        corr_data.to_excel(writer, sheet_name='Correlation')
        pvalue_matrix.to_excel(writer, sheet_name='P_values')
    
    print("✅ Table 5原始数据已导出 (相关矩阵)")
    return corr_data

def export_table6_data(df):
    """导出Table 6 - 双眼一致性分析"""
    
    # 计算双眼ICC
    from scipy.stats import pearsonr
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_外环颞侧', 'Outer Temporal'),
    ]
    
    results = []
    for col, name in indicators:
        # 获取双眼数据
        bilateral = df[df[col].notna()].groupby('Patient_ID').agg({
            col: ['count', 'mean', 'std'],
            'Eye': lambda x: list(x)
        })
        
        # 筛选有双眼数据的患者
        bilateral = bilateral[bilateral[(col, 'count')] == 2]
        
        if len(bilateral) > 10:
            # 获取左右眼数据
            left_eye = []
            right_eye = []
            
            for pid in bilateral.index:
                patient_data = df[(df['Patient_ID'] == pid) & (df[col].notna())]
                left = patient_data[patient_data['Eye'] == 'L'][col]
                right = patient_data[patient_data['Eye'] == 'R'][col]
                
                if len(left) > 0 and len(right) > 0:
                    left_eye.append(left.values[0])
                    right_eye.append(right.values[0])
            
            if len(left_eye) > 5:
                r, p = pearsonr(left_eye, right_eye)
                
                results.append({
                    'Parameter': name,
                    'n_patients': len(left_eye),
                    'Left_Mean': np.mean(left_eye),
                    'Left_SD': np.std(left_eye),
                    'Right_Mean': np.mean(right_eye),
                    'Right_SD': np.std(right_eye),
                    'Pearson_r': r,
                    'P_value': p,
                    'ICC_estimate': r,  # 简单估计
                })
    
    df_result = pd.DataFrame(results)
    df_result.to_excel(f'{OUTPUT_DIR}/Table6_Intereye_Consistency_RawData.xlsx', index=False)
    print("✅ Table 6原始数据已导出 (双眼一致性)")
    return df_result

if __name__ == "__main__":
    print("="*70)
    print("导出Table 3-6原始分析数据")
    print("="*70)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    df = load_data()
    
    print("\n导出Table 3-6:")
    export_table3_data(df)
    export_table4_data(df)
    export_table5_data(df)
    export_table6_data(df)
    
    print("\n" + "="*70)
    print("✅ Table 3-6原始数据导出完成!")
    print("="*70)