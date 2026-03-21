# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
import pandas as pd
import numpy as np
from scipy import stats
import warnings

def main():
    """主函数，包装原有执行代码"""
    warnings.filterwarnings('ignore')

    print("=" * 80)
    print("生成最终基线特征表")
    print("=" * 80)

    output_dir = '最终修改'

    # 读取数据
    print("读取数据...")
    oct_data = pd.read_excel(f'{output_dir}/OCT数据_最终清洗.xlsx')
    patient_list = pd.read_excel('名单/患者组名单.xlsx')
    control_list = pd.read_excel('名单/健康对照名单.xlsx')

    # 合并数据
    print("合并数据...")
    patient_oct = oct_data[oct_data['分组'] == '抑郁症'].copy()
    patient_merged = patient_oct.merge(
        patient_list[['姓名', '年龄', '性别', 'PHQ-9', 'GAD-7']],
        left_on='Patient_ID', right_on='姓名', how='left'
    )

    control_oct = oct_data[oct_data['分组'] == '健康对照'].copy()
    control_merged = control_oct.merge(
        control_list[['姓名', '年龄', '性别']],
        left_on='Patient_ID', right_on='姓名', how='left'
    )

    all_data = pd.concat([patient_merged, control_merged], ignore_index=True)

    dep = all_data[all_data['分组'] == '抑郁症']
    ctrl = all_data[all_data['分组'] == '健康对照']

    # 生成基线特征表
    baseline_table = []

    # 样本量
    baseline_table.append({
        '特征': '样本量 (n)',
        '抑郁症组': f"{len(dep)}",
        '对照组': f"{len(ctrl)}",
        'P值': '-'
    })

    # 年龄
    age_dep = dep['年龄'].dropna()
    age_ctrl = ctrl['年龄'].dropna()
    _, p_age = stats.mannwhitneyu(age_dep, age_ctrl)
    baseline_table.append({
        '特征': '年龄 (岁)',
        '抑郁症组': f"{age_dep.mean():.1f} ± {age_dep.std():.1f}",
        '对照组': f"{age_ctrl.mean():.1f} ± {age_ctrl.std():.1f}",
        'P值': f"{p_age:.4f}"
    })

    # 性别
    sex_dep_f = (dep['性别'] == '女').sum()
    sex_dep_m = (dep['性别'] == '男').sum()
    sex_ctrl_f = (ctrl['性别'] == '女').sum()
    sex_ctrl_m = (ctrl['性别'] == '男').sum()
    baseline_table.append({
        '特征': '性别 女/男 (%)',
        '抑郁症组': f"{sex_dep_f}/{sex_dep_m} ({sex_dep_f/len(dep)*100:.1f}%/{sex_dep_m/len(dep)*100:.1f}%)",
        '对照组': f"{sex_ctrl_f}/{sex_ctrl_m} ({sex_ctrl_f/len(ctrl)*100:.1f}%/{sex_ctrl_m/len(ctrl)*100:.1f}%)",
        'P值': '-'
    })

    # PHQ-9
    phq9_dep = dep['PHQ-9'].dropna()
    baseline_table.append({
        '特征': 'PHQ-9评分',
        '抑郁症组': f"{phq9_dep.mean():.2f} ± {phq9_dep.std():.2f} (n={len(phq9_dep)})",
        '对照组': '-',
        'P值': '-'
    })

    # GAD-7
    gad7_dep = dep['GAD-7'].dropna()
    baseline_table.append({
        '特征': 'GAD-7评分',
        '抑郁症组': f"{gad7_dep.mean():.2f} ± {gad7_dep.std():.2f} (n={len(gad7_dep)})",
        '对照组': '-',
        'P值': '-'
    })

    # OCT指标 - 黄斑
    retina_metrics = [
        ('黄斑平均厚度 (μm)', 'Retina_平均厚度'),
        ('黄斑总体积 (mm³)', 'Retina_总体积'),
        ('黄斑中心凹厚度 (μm)', 'Retina_黄斑中心凹')
    ]

    for label, col in retina_metrics:
        if col in all_data.columns:
            dep_vals = dep[col].dropna()
            ctrl_vals = ctrl[col].dropna()
            if len(dep_vals) > 0 and len(ctrl_vals) > 0:
                _, p = stats.mannwhitneyu(dep_vals, ctrl_vals)
                baseline_table.append({
                    '特征': label,
                    '抑郁症组': f"{dep_vals.mean():.2f} ± {dep_vals.std():.2f}",
                    '对照组': f"{ctrl_vals.mean():.2f} ± {ctrl_vals.std():.2f}",
                    'P值': f"{p:.4f}"
                })

    # OCT指标 - 视盘
    disc_metrics = [
        ('RNFL Total (μm)', 'RNFL_Total'),
        ('Disc Area (mm²)', 'Disc Area'),
        ('Cup Area (mm²)', 'Cup Area'),
        ('Rim Area (mm²)', 'Rim Area'),
        ('C/D Area Ratio', 'C/D Area Ratio')
    ]

    for label, col in disc_metrics:
        if col in all_data.columns:
            dep_vals = dep[col].dropna()
            ctrl_vals = ctrl[col].dropna()
            if len(dep_vals) > 0 and len(ctrl_vals) > 0:
                _, p = stats.mannwhitneyu(dep_vals, ctrl_vals)
                baseline_table.append({
                    '特征': label,
                    '抑郁症组': f"{dep_vals.mean():.2f} ± {dep_vals.std():.2f}",
                    '对照组': f"{ctrl_vals.mean():.2f} ± {ctrl_vals.std():.2f}",
                    'P值': f"{p:.4f}"
                })

    # 创建DataFrame
    baseline_df = pd.DataFrame(baseline_table)

    # 显示表格
    print("\n")
    print(baseline_df.to_string(index=False))

    # 保存为Excel
    baseline_df.to_excel(f'{output_dir}/Table1_基线特征表.xlsx', index=False)
    print(f"\n已保存: {output_dir}/Table1_基线特征表.xlsx")

    # 生成Markdown格式
    print("\nMarkdown格式（用于论文）:")
    print("| 特征 | 抑郁症组 (n={}) | 对照组 (n={}) | P值 |".format(len(dep), len(ctrl)))
    print("|------|----------------|--------------|-----|")
    for _, row in baseline_df.iterrows():
        print(f"| {row['特征']} | {row['抑郁症组']} | {row['对照组']} | {row['P值']} |")

    print("\n基线特征表生成完成!")



if __name__ == "__main__":
    main()