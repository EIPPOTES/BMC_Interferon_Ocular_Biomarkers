#!/usr/bin/env python3
"""
修复Figure 6 - 使用正确的箱线图格式
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

output_dir = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版'

def load_data():
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    return mdd

def fix_figure6_boxplot(mdd):
    """Figure 6: 使用正确的箱线图格式"""
    print("修复 Figure 6: 使用正确的箱线图格式...")
    
    mdd_with_phq9 = mdd[mdd['PHQ-9'].notna()].copy()
    
    # 按PHQ-9分层
    def classify_phq9(score):
        if pd.isna(score):
            return None
        elif score < 5:
            return 'Minimal\n(0-4)'
        elif score < 10:
            return 'Mild\n(5-9)'
        elif score < 15:
            return 'Moderate\n(10-14)'
        else:
            return 'Severe\n(≥15)'
    
    mdd_with_phq9['PHQ9_Group'] = mdd_with_phq9['PHQ-9'].apply(classify_phq9)
    
    # 准备数据用于箱线图
    groups = ['Minimal\n(0-4)', 'Mild\n(5-9)', 'Moderate\n(10-14)', 'Severe\n(≥15)']
    data_for_box = []
    group_labels = []
    stats_data = []
    
    print("  各组统计 (箱线图数据):")
    for group in groups:
        subset = mdd_with_phq9[mdd_with_phq9['PHQ9_Group'] == group]
        if len(subset) > 0:
            thickness = subset['Retina_平均厚度'].dropna()
            
            # 计算箱线图统计量
            q1 = thickness.quantile(0.25)
            q3 = thickness.quantile(0.75)
            median = thickness.median()
            iqr = q3 - q1
            lower_whisker = max(thickness.min(), q1 - 1.5 * iqr)
            upper_whisker = min(thickness.max(), q3 + 1.5 * iqr)
            
            n_participants = subset['Patient_ID'].nunique()
            n_eyes = len(thickness)
            
            stats_data.append({
                'group': group,
                'n_participants': n_participants,
                'n_eyes': n_eyes,
                'median': median,
                'q1': q1,
                'q3': q3,
                'iqr': iqr,
                'min': thickness.min(),
                'max': thickness.max(),
                'lower_whisker': lower_whisker,
                'upper_whisker': upper_whisker
            })
            
            print(f"    {group.replace(chr(10), ' ')}: N={n_participants}, n={n_eyes}")
            print(f"      Median={median:.1f}, Q1={q1:.1f}, Q3={q3:.1f}, IQR={iqr:.1f}")
            
            # 为箱线图准备数据
            for val in thickness:
                data_for_box.append(val)
                group_labels.append(group)
    
    # 创建DataFrame用于seaborn
    box_df = pd.DataFrame({
        'Thickness': data_for_box,
        'Group': group_labels
    })
    
    # 绘制
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # 使用seaborn绘制箱线图
    sns.boxplot(data=box_df, x='Group', y='Thickness', 
                palette=['#3498db', '#2ecc71', '#f39c12', '#e74c3c'],
                width=0.6, ax=ax)
    
    # 添加样本量标签
    for i, stat in enumerate(stats_data):
        ax.text(i, stat['upper_whisker'] + 3, 
               f'N={stat["n_participants"]}\nn={stat["n_eyes"]} eyes', 
               ha='center', fontsize=10, fontweight='bold')
    
    ax.set_xlabel('Depression Severity (PHQ-9)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean Macular Thickness (μm)', fontsize=12, fontweight='bold')
    ax.set_title('Figure 6. Distribution of Mean Macular Thickness by Depression Severity\n'
                'Box plots show median (horizontal line), interquartile range (box), '
                'and whiskers (1.5×IQR).',
                fontsize=13, fontweight='bold')
    
    # 添加说明
    ax.text(0.02, 0.98, 'Box = Q1-Q3 (IQR)\nLine = Median\nWhiskers = 1.5×IQR', 
           transform=ax.transAxes, fontsize=9, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure6_Subgroup_Analysis_箱线图版.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 同时创建一个带统计表格的版本
    fig2, (ax2, ax_table) = plt.subplots(2, 1, figsize=(12, 9), 
                                         gridspec_kw={'height_ratios': [3, 1]})
    
    # 上图：箱线图
    sns.boxplot(data=box_df, x='Group', y='Thickness', 
                palette=['#3498db', '#2ecc71', '#f39c12', '#e74c3c'],
                width=0.6, ax=ax2)
    
    for i, stat in enumerate(stats_data):
        ax2.text(i, stat['upper_whisker'] + 3, 
                f'N={stat["n_participants"]}', 
                ha='center', fontsize=10, fontweight='bold')
    
    ax2.set_xlabel('')
    ax2.set_ylabel('Mean Macular Thickness (μm)', fontsize=12, fontweight='bold')
    ax2.set_title('Figure 6. Distribution of Mean Macular Thickness by Depression Severity (MDD Patients)',
                  fontsize=13, fontweight='bold')
    
    # 下图：统计表格
    ax_table.axis('off')
    
    table_data = []
    for stat in stats_data:
        table_data.append([
            stat['group'].replace('\n', ' '),
            f"{stat['n_participants']}",
            f"{stat['n_eyes']}",
            f"{stat['median']:.1f}",
            f"{stat['q1']:.1f}-{stat['q3']:.1f}",
            f"{stat['min']:.1f}-{stat['max']:.1f}"
        ])
    
    table = ax_table.table(cellText=table_data,
                          colLabels=['PHQ-9 Group', 'N', 'n (eyes)', 'Median', 'IQR (Q1-Q3)', 'Range (Min-Max)'],
                          cellLoc='center',
                          loc='center',
                          bbox=[0, 0.3, 1, 0.6])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    for i in range(6):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure6_Subgroup_Analysis_箱线图带表格版.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✅ 已保存两个版本:")
    print(f"     - Figure6_Subgroup_Analysis_箱线图版.png (纯箱线图)")
    print(f"     - Figure6_Subgroup_Analysis_箱线图带表格版.png (带统计表格)")

def main():
    print("="*70)
    print("修复 Figure 6 - 使用正确的箱线图格式")
    print("="*70)
    
    mdd = load_data()
    fix_figure6_boxplot(mdd)
    
    print("\n" + "="*70)
    print("✅ Figure 6 修复完成!")
    print("="*70)
    print("""
改进内容:
  - 使用正确的箱线图格式 (Median, Q1-Q3 IQR, 1.5×IQR whiskers)
  - 不再是均值±SD的柱状图
  - 清晰显示数据分布和异常值
  - 添加统计表格版本便于查看具体数值
    """)

if __name__ == "__main__":
    main()