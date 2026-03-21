# 路径配置
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
#!/usr/bin/env python3
"""
提取所有Table和Figure的原始分析数据到Excel
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_curve, auc
import os

OUTPUT_DIR = '/mnt/c/Users/CUI/Desktop/投稿、数据修改/05_Raw_Data'

def load_data():
    """加载485眼数据"""
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx')
    control_ages = df[df['分组'] == '健康对照'].drop_duplicates('Patient_ID')[['Patient_ID', '年龄']]
    age_missing_ids = control_ages[control_ages['年龄'].isna()]['Patient_ID'].tolist()
    return df[~((df['分组'] == '健康对照') & (df['Patient_ID'].isin(age_missing_ids)))]

def calculate_cohens_d(x, y):
    """计算Cohen's d"""
    nx, ny = len(x), len(y)
    pooled_std = np.sqrt(((nx-1)*x.var() + (ny-1)*y.var()) / (nx + ny - 2))
    return (x.mean() - y.mean()) / pooled_std

def export_table1_data(df):
    """导出Table 1原始数据"""
    
    # 基线特征
    mdd_patients = df[df['分组'] == '抑郁症'].drop_duplicates('Patient_ID')
    control_patients = df[df['分组'] == '健康对照'].drop_duplicates('Patient_ID')
    
    data = {
        'Characteristic': [],
        'MDD_Mean': [],
        'MDD_SD': [],
        'Control_Mean': [],
        'Control_SD': [],
        'P_value': [],
    }
    
    # 年龄
    data['Characteristic'].append('Age (years)')
    data['MDD_Mean'].append(mdd_patients['年龄'].mean())
    data['MDD_SD'].append(mdd_patients['年龄'].std())
    data['Control_Mean'].append(control_patients['年龄'].mean())
    data['Control_SD'].append(control_patients['年龄'].std())
    _, p = stats.mannwhitneyu(mdd_patients['年龄'], control_patients['年龄'])
    data['P_value'].append(p)
    
    # 性别
    data['Characteristic'].append('Female (%)')
    data['MDD_Mean'].append((mdd_patients['性别'] == '女').mean() * 100)
    data['MDD_SD'].append(np.nan)
    data['Control_Mean'].append((control_patients['性别'] == '女').mean() * 100)
    data['Control_SD'].append(np.nan)
    data['P_value'].append(np.nan)
    
    # PHQ-9
    data['Characteristic'].append('PHQ-9 Score')
    data['MDD_Mean'].append(mdd_patients['PHQ-9'].mean())
    data['MDD_SD'].append(mdd_patients['PHQ-9'].std())
    data['Control_Mean'].append(np.nan)
    data['Control_SD'].append(np.nan)
    data['P_value'].append(np.nan)
    
    df_result = pd.DataFrame(data)
    df_result.to_excel(f'{OUTPUT_DIR}/Table1_Baseline_Characteristics_RawData.xlsx', index=False)
    print("✅ Table 1原始数据已导出")
    return df_result

def export_table2_data(df):
    """导出Table 2原始数据"""
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('Retina_内环颞侧', 'Inner Temporal Thickness'),
        ('Retina_外环上方', 'Outer Superior Thickness'),
        ('Retina_外环鼻侧', 'Outer Nasal Thickness'),
        ('Retina_外环下方', 'Outer Inferior Thickness'),
        ('RNFL_上方', 'Superior RNFL'),
        ('RNFL_颞侧', 'Temporal RNFL'),
        ('RNFL_鼻侧', 'Nasal RNFL'),
        ('RNFL_下方', 'Inferior RNFL'),
        ('Cup Area', 'Cup Area'),
        ('Rim Volume', 'Rim Volume'),
        ('C/D Area Ratio', 'C/D Area Ratio'),
    ]
    
    results = []
    for col, name in indicators:
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        if len(mdd_data) < 10 or len(control_data) < 10:
            continue
        
        stat, pvalue = stats.mannwhitneyu(mdd_data, control_data, alternative='two-sided')
        d = calculate_cohens_d(mdd_data, control_data)
        
        results.append({
            'Parameter': name,
            'MDD_n': len(mdd_data),
            'MDD_Mean': mdd_data.mean(),
            'MDD_SD': mdd_data.std(),
            'MDD_Median': mdd_data.median(),
            'MDD_Q1': mdd_data.quantile(0.25),
            'MDD_Q3': mdd_data.quantile(0.75),
            'Control_n': len(control_data),
            'Control_Mean': control_data.mean(),
            'Control_SD': control_data.std(),
            'Control_Median': control_data.median(),
            'Control_Q1': control_data.quantile(0.25),
            'Control_Q3': control_data.quantile(0.75),
            'P_value': pvalue,
            'Cohens_d': d,
        })
    
    df_result = pd.DataFrame(results)
    df_result.to_excel(f'{OUTPUT_DIR}/Table2_OCT_Parameters_RawData.xlsx', index=False)
    print("✅ Table 2原始数据已导出")
    return df_result

def export_figure3_data(df):
    """导出Figure 3原始数据 (ROC分析)"""
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('RNFL_上方', 'Superior RNFL'),
        ('RNFL_颞侧', 'Temporal RNFL'),
        ('Cup Area', 'Cup Area'),
        ('C/D Area Ratio', 'C/D Area Ratio'),
    ]
    
    results = []
    for col, name in indicators:
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        if len(mdd_data) < 10 or len(control_data) < 10:
            continue
        
        y_true = [1] * len(mdd_data) + [0] * len(control_data)
        y_scores = list(mdd_data) + list(control_data)
        
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        # 找到最佳截断点
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        
        results.append({
            'Parameter': name,
            'AUC': roc_auc,
            'Optimal_Threshold': optimal_threshold,
            'Sensitivity': tpr[optimal_idx],
            'Specificity': 1 - fpr[optimal_idx],
            'FPR_at_optimal': fpr[optimal_idx],
            'TPR_at_optimal': tpr[optimal_idx],
        })
    
    df_result = pd.DataFrame(results)
    df_result.to_excel(f'{OUTPUT_DIR}/Figure3_ROC_Analysis_RawData.xlsx', index=False)
    print("✅ Figure 3原始数据已导出")
    return df_result

def export_figure5_data(df):
    """导出Figure 5原始数据 (效应量)"""
    
    indicators = [
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_中心厚度', 'Central Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_内环鼻侧', 'Inner Nasal'),
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('RNFL_上方', 'Superior RNFL'),
        ('RNFL_颞侧', 'Temporal RNFL'),
        ('RNFL_鼻侧', 'Nasal RNFL'),
        ('Cup Area', 'Cup Area'),
        ('Rim Volume', 'Rim Volume'),
        ('C/D Area Ratio', 'C/D Area Ratio'),
    ]
    
    results = []
    for col, name in indicators:
        mdd_data = df[(df['分组'] == '抑郁症') & (df[col].notna())][col]
        control_data = df[(df['分组'] == '健康对照') & (df[col].notna())][col]
        
        if len(mdd_data) < 10 or len(control_data) < 10:
            continue
        
        d = calculate_cohens_d(mdd_data, control_data)
        
        # 95% CI
        nx, ny = len(mdd_data), len(control_data)
        se = np.sqrt((nx + ny) / (nx * ny) + d**2 / (2 * (nx + ny)))
        ci_lower = d - 1.96 * se
        ci_upper = d + 1.96 * se
        
        results.append({
            'Parameter': name,
            'MDD_n': len(mdd_data),
            'MDD_Mean': mdd_data.mean(),
            'MDD_SD': mdd_data.std(),
            'Control_n': len(control_data),
            'Control_Mean': control_data.mean(),
            'Control_SD': control_data.std(),
            'Mean_Difference': mdd_data.mean() - control_data.mean(),
            'Cohens_d': d,
            'CI_Lower': ci_lower,
            'CI_Upper': ci_upper,
            'Interpretation': 'Small' if abs(d) < 0.5 else 'Medium' if abs(d) < 0.8 else 'Large',
        })
    
    df_result = pd.DataFrame(results)
    df_result.to_excel(f'{OUTPUT_DIR}/Figure5_Effect_Sizes_RawData.xlsx', index=False)
    print("✅ Figure 5原始数据已导出")
    return df_result

def export_correlation_data(df):
    """导出相关性分析原始数据 (Figure 4)"""
    
    mdd_df = df[df['分组'] == '抑郁症'].copy()
    
    indicators = [
        'Retina_平均厚度',
        'Retina_外环颞侧',
        'Retina_内环颞侧',
        'RNFL_上方',
        'RNFL_颞侧',
    ]
    
    results = []
    for col in indicators:
        valid_data = mdd_df[(mdd_df[col].notna()) & (mdd_df['PHQ-9'].notna())]
        
        if len(valid_data) < 10:
            continue
        
        x = valid_data['PHQ-9']
        y = valid_data[col]
        
        r, pvalue = stats.pearsonr(x, y)
        
        # 线性回归
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        
        results.append({
            'Parameter': col,
            'n': len(valid_data),
            'Pearson_r': r,
            'P_value': pvalue,
            'R_squared': r_value**2,
            'Slope': slope,
            'Intercept': intercept,
            'Std_Error': std_err,
        })
    
    df_result = pd.DataFrame(results)
    df_result.to_excel(f'{OUTPUT_DIR}/Figure4_Correlation_RawData.xlsx', index=False)
    print("✅ Figure 4原始数据已导出")
    return df_result

if __name__ == "__main__":
    print("="*70)
    print("导出所有Table和Figure的原始分析数据")
    print("="*70)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    df = load_data()
    
    print("\n导出Table数据:")
    export_table1_data(df)
    export_table2_data(df)
    
    print("\n导出Figure数据:")
    export_figure3_data(df)
    export_correlation_data(df)
    export_figure5_data(df)
    
    print("\n" + "="*70)
    print("✅ 所有原始数据已导出!")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)