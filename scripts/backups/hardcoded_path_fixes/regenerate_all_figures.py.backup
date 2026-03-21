#!/usr/bin/env python3
"""
根据原始数据重新生成所有Figures
使用正确的数据：02_OCT数据_完整整合.xlsx
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import roc_curve, auc
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# 设置输出目录
output_dir = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_重新生成_全面版'
import os
os.makedirs(output_dir, exist_ok=True)

def load_data():
    """加载原始数据"""
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    
    # 按人分组
    df_patient = df.drop_duplicates(subset=['Patient_ID'])
    mdd = df_patient[df_patient['分组'] == '抑郁症']
    control = df_patient[df_patient['分组'] == '健康对照']
    
    print(f"数据加载完成:")
    print(f"  总样本: {len(df_patient)}人")
    print(f"  MDD组: {len(mdd)}人")
    print(f"  对照组: {len(control)}人")
    
    return df, df_patient, mdd, control

def generate_figure1(mdd, control):
    """Figure 1: 研究流程图 (简化版)"""
    print("\n生成 Figure 1: Study Flowchart...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.axis('off')
    
    # 流程图数据
    total_screened = 280
    excluded = 29
    final_sample = 251
    mdd_n = len(mdd)
    control_n = len(control)
    
    # 绘制流程图
    ax.text(0.5, 0.9, 'Study Flowchart', fontsize=16, fontweight='bold', ha='center')
    ax.text(0.5, 0.8, f'Initial Assessment\nN = {total_screened}', fontsize=12, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    ax.text(0.75, 0.7, f'Excluded (n={excluded})\n- Ocular disease\n- Poor image quality\n- Other reasons', 
            fontsize=9, ha='left', style='italic')
    ax.text(0.5, 0.6, f'Final Sample\nN = {final_sample} ({mdd_n} MDD + {control_n} Control)', 
            fontsize=12, ha='center', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    ax.text(0.5, 0.45, f'MDD Group: {mdd_n} participants', fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))
    ax.text(0.5, 0.35, f'Control Group: {control_n} participants', fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure1_Study_Flowchart.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure1_Study_Flowchart.png")

def generate_figure2(df, mdd, control):
    """Figure 2: 全面的组间比较"""
    print("\n生成 Figure 2: Group Comparison (全面版)...")
    
    # 选择关键参数 (15个参数)
    key_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('Retina_外环下方', 'Outer Inferior'),
        ('RNFL_Total', 'Total RNFL'),
        ('RNFL_上方', 'Superior RNFL'),
        ('RNFL_颞侧', 'Temporal RNFL'),
        ('GCL+_平均厚度', 'GCL+ Mean'),
        ('GCL++_平均厚度', 'GCL++ Mean'),
        ('Cup Area', 'Cup Area'),
        ('Rim Area', 'Rim Area'),
        ('C/D Area Ratio', 'C/D Ratio'),
        ('Rim Volume', 'Rim Volume'),
        ('Choroid_平均厚度', 'Choroid Mean'),
    ]
    
    # 创建大图
    fig, axes = plt.subplots(5, 3, figsize=(18, 24))
    axes = axes.flatten()
    
    for idx, (col, name) in enumerate(key_params):
        if col in df.columns:
            ax = axes[idx]
            
            mdd_data = mdd[col].dropna()
            ctrl_data = control[col].dropna()
            
            # 计算统计量
            stat, p = stats.mannwhitneyu(mdd_data, ctrl_data, alternative='two-sided')
            
            # 绘制箱线图
            data_to_plot = [ctrl_data, mdd_data]
            bp = ax.boxplot(data_to_plot, labels=['Control', 'MDD'], patch_artist=True)
            
            # 设置颜色
            bp['boxes'][0].set_facecolor('lightblue')
            bp['boxes'][1].set_facecolor('lightcoral')
            
            # 添加统计信息
            ax.set_title(f'{name}\nP={p:.4f}', fontsize=10)
            ax.set_ylabel('Thickness (μm)')
    
    # 删除多余的子图
    for idx in range(len(key_params), len(axes)):
        fig.delaxes(axes[idx])
    
    plt.suptitle('Figure 2. Comprehensive Group Comparison (15 Parameters)', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure2_Group_Comparison_全面版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure2_Group_Comparison_全面版.png (15个参数)")

def generate_figure3(df, mdd, control):
    """Figure 3: ROC曲线 (全面版)"""
    print("\n生成 Figure 3: ROC Curves (全面版)...")
    
    # 准备数据
    df['分组_编码'] = (df['分组'] == '抑郁症').astype(int)
    df_patient = df.groupby('Patient_ID').agg({
        'Retina_外环颞侧': 'mean',
        'Retina_平均厚度': 'mean',
        'Retina_内环颞侧': 'mean',
        'RNFL_Total': 'mean',
        '分组_编码': 'first'
    }).dropna()
    
    y_true = df_patient['分组_编码'].values
    
    # 选择参数
    roc_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('RNFL_Total', 'Total RNFL'),
    ]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = ['blue', 'red', 'green', 'orange']
    
    for (col, name), color in zip(roc_params, colors):
        y_scores = -df_patient[col].values
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        ax.plot(fpr, tpr, color=color, lw=2, 
                label=f'{name} (AUC = {roc_auc:.3f})')
    
    # 对角线
    ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', label='Random')
    
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=12)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax.set_title('Figure 3. ROC Curves for OCT Parameters', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure3_ROC_Curves_全面版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure3_ROC_Curves_全面版.png (4个参数)")

def generate_figure4(mdd):
    """Figure 4: 相关性散点图"""
    print("\n生成 Figure 4: Correlation Analysis...")
    
    mdd_with_phq9 = mdd[mdd['PHQ-9'].notna()]
    
    corr_params = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('RNFL_Total', 'Total RNFL Thickness'),
    ]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, (col, name) in enumerate(corr_params):
        ax = axes[idx]
        
        x = mdd_with_phq9['PHQ-9']
        y = mdd_with_phq9[col]
        
        # 移除缺失值
        mask = x.notna() & y.notna()
        x_clean = x[mask]
        y_clean = y[mask]
        
        if len(x_clean) > 2:
            r, p = stats.spearmanr(x_clean, y_clean)
            
            # 绘制散点图
            ax.scatter(x_clean, y_clean, alpha=0.6, s=50)
            
            # 添加回归线
            z = np.polyfit(x_clean, y_clean, 1)
            p_fit = np.poly1d(z)
            ax.plot(x_clean, p_fit(x_clean), "r--", alpha=0.8)
            
            ax.set_xlabel('PHQ-9 Score', fontsize=11)
            ax.set_ylabel(f'{name} (μm)', fontsize=11)
            ax.set_title(f'{name}\nr = {r:.3f}, P = {p:.3f}', fontsize=11)
            ax.grid(True, alpha=0.3)
    
    plt.suptitle('Figure 4. Correlation between PHQ-9 and OCT Parameters', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure4_Correlation_Scatter_全面版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure4_Correlation_Scatter_全面版.png (3个参数)")

def generate_figure5(df, mdd, control):
    """Figure 5: 效应量森林图"""
    print("\n生成 Figure 5: Forest Plot...")
    
    # 计算效应量
    def cohens_d(x, y):
        nx = len(x)
        ny = len(y)
        dof = nx + ny - 2
        return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1)**2 + (ny-1)*np.std(y, ddof=1)**2) / dof)
    
    effect_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('RNFL_Total', 'Total RNFL'),
        ('Cup Area', 'Cup Area'),
        ('Rim Volume', 'Rim Volume'),
    ]
    
    effects = []
    labels = []
    
    for col, name in effect_params:
        if col in df.columns:
            mdd_data = mdd[col].dropna()
            ctrl_data = control[col].dropna()
            d = cohens_d(mdd_data, ctrl_data)
            effects.append(d)
            labels.append(name)
    
    # 创建森林图
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_pos = np.arange(len(labels))
    colors = ['red' if e < 0 else 'blue' for e in effects]
    
    ax.barh(y_pos, effects, color=colors, alpha=0.7, edgecolor='black')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
    ax.axvline(x=-0.2, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(x=0.2, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # 添加数值标签
    for i, (effect, label) in enumerate(zip(effects, labels)):
        ax.text(effect + 0.02 if effect > 0 else effect - 0.02, i, 
                f'{effect:.2f}', va='center', fontsize=9,
                ha='left' if effect > 0 else 'right')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12)
    ax.set_title('Figure 5. Effect Sizes of Retinal Changes in MDD', fontsize=14, fontweight='bold')
    ax.grid(True, axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure5_Forest_Plot_全面版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure5_Forest_Plot_全面版.png ({len(effects)}个参数)")

def generate_figure6(mdd):
    """Figure 6: 亚组分析"""
    print("\n生成 Figure 6: Subgroup Analysis...")
    
    # 按PHQ-9分层
    def classify_phq9(score):
        if pd.isna(score):
            return 'Unknown'
        elif score < 5:
            return 'Minimal\n(0-4)'
        elif score < 10:
            return 'Mild\n(5-9)'
        elif score < 15:
            return 'Moderate\n(10-14)'
        else:
            return 'Severe\n(≥15)'
    
    mdd_analysis = mdd.copy()
    mdd_analysis['PHQ9_Group'] = mdd_analysis['PHQ-9'].apply(classify_phq9)
    
    # 选择主要参数
    param = 'Retina_平均厚度'
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    groups = ['Minimal\n(0-4)', 'Mild\n(5-9)', 'Moderate\n(10-14)', 'Severe\n(≥15)']
    means = []
    stds = []
    ns = []
    
    for group in groups:
        subset = mdd_analysis[mdd_analysis['PHQ9_Group'] == group]
        if len(subset) > 0:
            values = subset[param].dropna()
            means.append(values.mean())
            stds.append(values.std())
            ns.append(len(values))
        else:
            means.append(0)
            stds.append(0)
            ns.append(0)
    
    x_pos = np.arange(len(groups))
    bars = ax.bar(x_pos, means, yerr=stds, capsize=5, 
                  color=['lightblue', 'lightgreen', 'orange', 'lightcoral'],
                  alpha=0.7, edgecolor='black')
    
    # 添加样本量标签
    for i, (mean, std, n) in enumerate(zip(means, stds, ns)):
        if n > 0:
            ax.text(i, mean + std + 2, f'n={n}', ha='center', fontsize=9)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(groups)
    ax.set_ylabel('Mean Macular Thickness (μm)', fontsize=12)
    ax.set_title('Figure 6. OCT Parameters by Depression Severity', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Figure6_Subgroup_Analysis_全面版.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 已保存: Figure6_Subgroup_Analysis_全面版.png")

def main():
    print("="*70)
    print("重新生成所有Figures (全面版)")
    print("="*70)
    print(f"输出目录: {output_dir}")
    
    # 加载数据
    df, df_patient, mdd, control = load_data()
    
    # 生成所有Figures
    generate_figure1(mdd, control)
    generate_figure2(df, mdd, control)
    generate_figure3(df, mdd, control)
    generate_figure4(mdd)
    generate_figure5(df, mdd, control)
    generate_figure6(mdd)
    
    print("\n" + "="*70)
    print("✅ 所有Figures生成完成!")
    print("="*70)
    
    print(f"\n生成的文件:")
    for f in os.listdir(output_dir):
        if f.endswith('.png'):
            size = os.path.getsize(os.path.join(output_dir, f)) / 1024
            print(f"  ✅ {f} ({size:.1f} KB)")
    
    print(f"\n所有Figures均为300 DPI，符合期刊标准")

if __name__ == "__main__":
    main()