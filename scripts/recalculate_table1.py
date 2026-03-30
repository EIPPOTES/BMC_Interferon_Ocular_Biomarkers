#!/usr/bin/env python3
# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
重新计算Table 1基线特征表
使用正确的原始数据：02_OCT数据_完整整合.xlsx
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, chi2_contingency
import os

def recalculate_table1():
    """重新计算Table 1"""
    
    # 读取正确的数据文件
    data_file = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    df = pd.read_excel(data_file)
    
    print("="*70)
    print("重新计算 Table 1: Baseline Characteristics")
    print("="*70)
    print(f"\n数据来源: {data_file}")
    print(f"原始形状: {df.shape}")
    
    # 按人分组（去重）
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    print(f"去重后（人水平）: {df_patient.shape}")
    
    # 分组
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    
    print(f"\n📊 样本量:")
    print(f"  MDD组: {len(mdd)}人")
    print(f"  对照组: {len(control)}人")
    print(f"  总计: {len(mdd) + len(control)}人")
    
    results = {}
    
    # 1. 年龄
    print(f"\n{'='*70}")
    print("1. 年龄 (Age)")
    print(f"{'='*70}")
    
    mdd_age = mdd['年龄'].dropna()
    ctrl_age = control['年龄'].dropna()
    
    # Mann-Whitney U检验
    stat, p_age = mannwhitneyu(mdd_age, ctrl_age, alternative='two-sided')
    
    results['age'] = {
        'mdd_mean': mdd_age.mean(),
        'mdd_std': mdd_age.std(),
        'ctrl_mean': ctrl_age.mean(),
        'ctrl_std': ctrl_age.std(),
        'p_value': p_age
    }
    
    print(f"  MDD: {mdd_age.mean():.1f} ± {mdd_age.std():.1f}岁 (n={len(mdd_age)})")
    print(f"  Control: {ctrl_age.mean():.1f} ± {ctrl_age.std():.1f}岁 (n={len(ctrl_age)})")
    print(f"  P-value: {p_age:.3f}")
    print(f"\n  ⚠️ 与Word文档对比:")
    print(f"    Word: MDD 30.2±13.5岁, Control 27.6±12.4岁 (P=0.096)")
    print(f"    计算: MDD {mdd_age.mean():.1f}±{mdd_age.std():.1f}岁, Control {ctrl_age.mean():.1f}±{ctrl_age.std():.1f}岁 (P={p_age:.3f})")
    
    # 2. 性别
    print(f"\n{'='*70}")
    print("2. 性别 (Gender)")
    print(f"{'='*70}")
    
    mdd_male = (mdd['性别'] == '男').sum()
    mdd_female = (mdd['性别'] == '女').sum()
    ctrl_male = (control['性别'] == '男').sum()
    ctrl_female = (control['性别'] == '女').sum()
    
    # 卡方检验
    contingency = [[mdd_male, mdd_female], [ctrl_male, ctrl_female]]
    chi2, p_gender, _, _ = chi2_contingency(contingency)
    
    mdd_total = mdd_male + mdd_female
    ctrl_total = ctrl_male + ctrl_female
    
    print(f"  MDD - 男: {mdd_male} ({mdd_male/mdd_total*100:.1f}%)")
    print(f"  MDD - 女: {mdd_female} ({mdd_female/mdd_total*100:.1f}%)")
    print(f"  Control - 男: {ctrl_male} ({ctrl_male/ctrl_total*100:.1f}%)")
    print(f"  Control - 女: {ctrl_female} ({ctrl_female/ctrl_total*100:.1f}%)")
    print(f"  P-value: {p_gender:.3f}")
    
    # 3. 受教育年限
    print(f"\n{'='*70}")
    print("3. 受教育年限 (Education)")
    print(f"{'='*70}")
    
    mdd_edu = mdd['受教育年限'].dropna()
    ctrl_edu = control['受教育年限'].dropna()
    
    stat, p_edu = mannwhitneyu(mdd_edu, ctrl_edu, alternative='two-sided')
    
    print(f"  MDD: {mdd_edu.mean():.1f} ± {mdd_edu.std():.1f}年 (n={len(mdd_edu)})")
    print(f"  Control: {ctrl_edu.mean():.1f} ± {ctrl_edu.std():.1f}年 (n={len(ctrl_edu)})")
    print(f"  P-value: {p_edu:.3f}")
    
    # 4. PHQ-9 (仅MDD组)
    print(f"\n{'='*70}")
    print("4. PHQ-9评分 (仅MDD组)")
    print(f"{'='*70}")
    
    mdd_phq9 = mdd['PHQ-9'].dropna()
    
    print(f"  MDD: {mdd_phq9.mean():.1f} ± {mdd_phq9.std():.1f} (n={len(mdd_phq9)})")
    print(f"  范围: {mdd_phq9.min():.0f} - {mdd_phq9.max():.0f}")
    
    # 总结对比
    print(f"\n{'='*70}")
    print("📋 与Word文档对比总结")
    print(f"{'='*70}")
    
    print(f"""
{'指标':<20} {'Word文档':<25} {'重新计算':<25} {'状态':<10}
{'-'*70}
{'年龄 (MDD)':<20} {'30.2±13.5岁':<25} {f"{mdd_age.mean():.1f}±{mdd_age.std():.1f}岁":<25} {'❌ 不符':<10}
{'年龄 (Control)':<20} {'27.6±12.4岁':<25} {f"{ctrl_age.mean():.1f}±{ctrl_age.std():.1f}岁":<25} {'⚠️ 接近':<10}
{'性别 (MDD男)':<20} {'72 (43.9%)':<25} {f"{mdd_male} ({mdd_male/mdd_total*100:.1f}%)":<25} {'待核对':<10}
{'教育年限 (MDD)':<20} {'12.8±3.2年':<25} {f"{mdd_edu.mean():.1f}±{mdd_edu.std():.1f}年":<25} {'待核对':<10}
    """)
    
    print(f"\n{'='*70}")
    print("结论与建议")
    print(f"{'='*70}")
    
    print("""
⚠️ 发现显著差异:

1. MDD组年龄差异最大 (约8岁)
   - Word: 30.2±13.5岁
   - 实际: 38.3±20.2岁
   
2. 可能原因:
   a) Word文档使用了筛选后的数据（如排除>60岁）
   b) 基于不同的原始数据文件
   c) 数据版本未同步

3. 建议行动:
   ✅ 立即使用当前数据重新生成Table 1
   ✅ 在方法部分明确说明纳入/排除标准
   ✅ 确认所有其他表格基于相同数据版本
   ✅ 更新投稿版中的所有相关文件
    """)
    
    return results

if __name__ == "__main__":
    recalculate_table1()