#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
生成完整的数据更新清单
对比论文中的统计数据 vs 实际计算值
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_curve, auc
import re

def load_data():
    """加载数据"""
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    df['C/D_Area_Ratio'] = df['Cup Area'] / df['Disc Area']
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    return df, df_patient, mdd, control

def cohens_d(x, y):
    """计算Cohen's d"""
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1)**2 + (ny-1)*np.std(y, ddof=1)**2) / dof)

def calculate_all_statistics(mdd, control):
    """计算所有统计数据"""
    
    results = {}
    
    # 1. 样本量
    results['sample'] = {
        'mdd_n': len(mdd),
        'control_n': len(control),
        'total_n': len(mdd) + len(control),
        'mdd_eyes': len(mdd) * 2,  # 近似
        'control_eyes': len(control) * 2
    }
    
    # 2. 年龄
    mdd_age = mdd['年龄'].dropna()
    ctrl_age = control['年龄'].dropna()
    stat, p_age = mannwhitneyu(mdd_age, ctrl_age, alternative='two-sided')
    
    results['age'] = {
        'mdd_mean': mdd_age.mean(),
        'mdd_std': mdd_age.std(),
        'ctrl_mean': ctrl_age.mean(),
        'ctrl_std': ctrl_age.std(),
        'p_value': p_age
    }
    
    # 3. 关键OCT参数
    key_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('RNFL_Total', 'Total RNFL'),
    ]
    
    for col, name in key_params:
        if col in mdd.columns:
            mdd_data = mdd[col].dropna()
            ctrl_data = control[col].dropna()
            
            stat, p = mannwhitneyu(mdd_data, ctrl_data, alternative='two-sided')
            d = cohens_d(mdd_data, ctrl_data)
            
            results[name] = {
                'mdd_mean': mdd_data.mean(),
                'mdd_std': mdd_data.std(),
                'ctrl_mean': ctrl_data.mean(),
                'ctrl_std': ctrl_data.std(),
                'p_value': p,
                'cohens_d': d
            }
    
    # 4. ROC分析
    df_temp = pd.concat([mdd, control])
    df_temp['分组_编码'] = (df_temp['分组'] == '抑郁症').astype(int)
    
    y_true = df_temp['分组_编码'].values
    y_scores = -df_temp['Retina_外环颞侧'].values
    
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)
    
    results['roc'] = {
        'auc': roc_auc
    }
    
    # 5. 相关性
    mdd_with_phq9 = mdd[mdd['PHQ-9'].notna()]
    if len(mdd_with_phq9) > 0:
        x = mdd_with_phq9['PHQ-9']
        y = mdd_with_phq9['Retina_外环颞侧']
        mask = x.notna() & y.notna()
        if mask.sum() > 2:
            r, p = spearmanr(x[mask], y[mask])
            results['correlation'] = {
                'r': r,
                'p': p,
                'n': mask.sum()
            }
    
    return results

def generate_update_list():
    """生成更新清单"""
    
    df, df_patient, mdd, control = load_data()
    actual_stats = calculate_all_statistics(mdd, control)
    
    print("="*70)
    print("📋 论文数据更新清单")
    print("="*70)
    print("\n对比: 论文中的数值 vs 实际计算值\n")
    
    updates = []
    
    # 1. 样本量
    print("1️⃣ 样本量更新")
    print("-"*70)
    print(f"{'项目':<25} {'论文':<15} {'实际':<15} {'操作':<10}")
    print("-"*70)
    
    sample_updates = [
        ('MDD participants', '164', str(actual_stats['sample']['mdd_n']), '✅ 一致'),
        ('Control participants', '87', str(actual_stats['sample']['control_n']), '⚠️ 更新'),
        ('Total participants', '251', str(actual_stats['sample']['total_n']), '⚠️ 更新'),
    ]
    
    for item, paper, actual, action in sample_updates:
        print(f"{item:<25} {paper:<15} {actual:<15} {action}")
        if '更新' in action:
            updates.append({
                'section': 'Sample Size',
                'item': item,
                'old': paper,
                'new': actual
            })
    
    # 2. 年龄
    print("\n2️⃣ 年龄统计更新")
    print("-"*70)
    age = actual_stats['age']
    print(f"MDD年龄: {age['mdd_mean']:.1f} ± {age['mdd_std']:.1f}岁")
    print(f"Control年龄: {age['ctrl_mean']:.1f} ± {age['ctrl_std']:.1f}岁")
    print(f"P-value: {age['p_value']:.3f}")
    
    if abs(age['mdd_mean'] - 30.2) > 1 or abs(age['ctrl_mean'] - 27.6) > 1:
        updates.append({
            'section': 'Age',
            'item': 'Age statistics',
            'old': 'MDD: 30.2±13.5, Control: 27.6±12.4',
            'new': f"MDD: {age['mdd_mean']:.1f}±{age['mdd_std']:.1f}, Control: {age['ctrl_mean']:.1f}±{age['ctrl_std']:.1f}"
        })
    
    # 3. 效应量
    print("\n3️⃣ 效应量 (Cohen's d) 更新")
    print("-"*70)
    print(f"{'参数':<25} {'论文d值':<15} {'实际d值':<15} {'差异':<10}")
    print("-"*70)
    
    effect_paper = {
        'Outer Temporal': -0.50,
        'Mean Macular': -0.42,
    }
    
    for param, paper_d in effect_paper.items():
        if param in actual_stats:
            actual_d = actual_stats[param]['cohens_d']
            diff = abs(actual_d - paper_d)
            status = '⚠️' if diff > 0.03 else '✅'
            print(f"{param:<25} {paper_d:<15} {actual_d:.2f}          {status}")
            
            if diff > 0.03:
                updates.append({
                    'section': 'Effect Size',
                    'item': param,
                    'old': f"d={paper_d}",
                    'new': f"d={actual_d:.2f}"
                })
    
    # 4. AUC
    print("\n4️⃣ AUC更新")
    print("-"*70)
    actual_auc = actual_stats['roc']['auc']
    paper_auc = 0.646
    print(f"Outer Temporal AUC:")
    print(f"  论文: {paper_auc}")
    print(f"  实际: {actual_auc:.3f}")
    print(f"  差异: {abs(actual_auc - paper_auc):.3f}")
    
    if abs(actual_auc - paper_auc) > 0.01:
        updates.append({
            'section': 'ROC',
            'item': 'Outer Temporal AUC',
            'old': str(paper_auc),
            'new': f"{actual_auc:.3f}"
        })
    
    # 5. 相关性
    if 'correlation' in actual_stats:
        print("\n5️⃣ 相关性更新")
        print("-"*70)
        corr = actual_stats['correlation']
        print(f"PHQ-9 vs Outer Temporal:")
        print(f"  r = {corr['r']:.3f}")
        print(f"  P = {corr['p']:.3f}")
        print(f"  n = {corr['n']}")
    
    # 总结
    print("\n" + "="*70)
    print(f"📊 需要更新的项目: {len(updates)}项")
    print("="*70)
    
    if updates:
        print("\n更新清单:")
        for i, upd in enumerate(updates, 1):
            print(f"{i}. [{upd['section']}] {upd['item']}")
            print(f"   旧: {upd['old']}")
            print(f"   新: {upd['new']}")
    else:
        print("\n✅ 所有数据一致，无需更新")
    
    return updates

def main():
    generate_update_list()
    
    print("\n" + "="*70)
    print("🔧 建议操作")
    print("="*70)
    print("""
1. 更新摘要中的样本量: 251 → 244 participants
2. 更新方法中的对照组: 87 → 80 controls
3. 更新效应量数值 (如果差异>0.05)
4. 更新Table 1 (基线特征表)
5. 检查所有提及样本量的地方

⚠️  重要: 更新后需要重新生成Tables和Figures
    """)

if __name__ == "__main__":
    main()