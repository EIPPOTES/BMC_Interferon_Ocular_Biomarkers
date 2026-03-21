#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
制作SCI图表 - Figure 2: Group comparison of key OCT parameters
基于463眼OCT-MDD研究数据
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
from scipy import stats
import matplotlib

# SCI图表设置
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'grid.linewidth': 0.3
})

# 色盲友好配色方案
colors = {
    'mdd': '#D55E00',  # 橙色 - MDD组
    'control': '#0072B2',  # 蓝色 - 对照组
    'significant': '#CC79A7',  # 粉色 - 显著
    'ns': '#999999',  # 灰色 - 不显著
    'retina': '#E69F00',  # 黄色 - 视网膜参数
    'rnfl': '#56B4E9',  # 浅蓝 - RNFL参数
    'gcl': '#009E73',  # 绿色 - GCL参数
    'disc': '#F0E442',  # 黄色 - 视盘参数
    'choroid': '#0072B2'  # 蓝色 - 脉络膜参数
}

def load_group_comparison_data():
    """加载组间比较数据"""
    data_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/03_Tables/Group_Comparison_463eyes_20260315.xlsx"
    df = pd.read_excel(data_path)
    
    print(f"加载组间比较数据: {len(df)} 个指标")
    print(f"显著指标数 (P<0.05): {len(df[df['P_value'] < 0.05])}")
    print(f"高度显著指标数 (P<0.001): {len(df[df['P_value'] < 0.001])}")
    
    # 计算95%置信区间 (基于Cohen's d)
    # 对于Cohen's d, 95% CI ≈ d ± 1.96 * SE, 其中SE = sqrt((n1+n2)/(n1*n2) + d²/(2*(n1+n2)))
    df['Cohens_d_CI_lower'] = df['Cohens_d'] - 1.96 * np.sqrt(
        (df['MDD_n'] + df['Control_n']) / (df['MDD_n'] * df['Control_n']) + 
        df['Cohens_d']**2 / (2 * (df['MDD_n'] + df['Control_n']))
    )
    df['Cohens_d_CI_upper'] = df['Cohens_d'] + 1.96 * np.sqrt(
        (df['MDD_n'] + df['Control_n']) / (df['MDD_n'] * df['Control_n']) + 
        df['Cohens_d']**2 / (2 * (df['MDD_n'] + df['Control_n']))
    )
    
    # 添加类别信息
    def categorize_parameter(param):
        param_lower = param.lower()
        if 'retina' in param_lower:
            return 'Retina'
        elif 'rnfl' in param_lower:
            return 'RNFL'
        elif 'gcl' in param_lower:
            return 'GCL'
        elif any(kw in param_lower for kw in ['disc', 'cup', 'rim', 'area', 'ratio']):
            return 'Optic Disc'
        elif 'choroid' in param_lower:
            return 'Choroid'
        else:
            return 'Other'
    
    df['Category'] = df['Parameter'].apply(categorize_parameter)
    
    # 添加显著性标记
    df['Significance_level'] = 'ns'
    df.loc[df['P_value'] < 0.05, 'Significance_level'] = '*'
    df.loc[df['P_value'] < 0.01, 'Significance_level'] = '**'
    df.loc[df['P_value'] < 0.001, 'Significance_level'] = '***'
    
    # 选择前15个效应量最大的指标
    df_sorted = df.sort_values('Cohens_d').head(15).copy()
    
    # 反转顺序用于绘图（从上到下）
    df_sorted = df_sorted.iloc[::-1].reset_index(drop=True)
    
    return df_sorted

def create_forest_plot(df):
    """创建森林图"""
    # 创建图形 - 单栏宽度 (8.5 cm = 3.35 inches)
    fig, ax = plt.subplots(figsize=(3.35, 4.5))  # 宽度: 单栏, 高度: 根据内容调整
    
    # y轴位置（从上到下）
    y_positions = np.arange(len(df))
    
    # 绘制效应量点
    for i, (idx, row) in enumerate(df.iterrows()):
        # 确定颜色
        if row['Category'] == 'Retina':
            color = colors['retina']
        elif row['Category'] == 'RNFL':
            color = colors['rnfl']
        elif row['Category'] == 'GCL':
            color = colors['gcl']
        elif row['Category'] == 'Optic Disc':
            color = colors['disc']
        else:
            color = colors['choroid']
        
        # 绘制效应量点
        ax.plot(row['Cohens_d'], i, 'o', 
                markersize=8, 
                color=color,
                markeredgecolor='black',
                markeredgewidth=0.5,
                zorder=3)
        
        # 绘制置信区间线
        ax.plot([row['Cohens_d_CI_lower'], row['Cohens_d_CI_upper']], [i, i], 
                color=color, linewidth=1.5, zorder=2)
        
        # 绘制置信区间端帽
        ax.plot([row['Cohens_d_CI_lower'], row['Cohens_d_CI_lower']], [i-0.1, i+0.1], 
                color=color, linewidth=1.5, zorder=2)
        ax.plot([row['Cohens_d_CI_upper'], row['Cohens_d_CI_upper']], [i-0.1, i+0.1], 
                color=color, linewidth=1.5, zorder=2)
        
        # 添加显著性标记
        if row['Significance_level'] != 'ns':
            ax.text(row['Cohens_d'] + 0.02, i, row['Significance_level'], 
                   fontsize=7, fontweight='bold', va='center')
    
    # 添加零效应参考线
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5, zorder=1)
    
    # 设置y轴标签（参数名称）
    # 简化参数名称用于显示
    def simplify_param_name(param):
        # 移除常见前缀
        replacements = [
            ('Retina_', ''),
            ('RNFL_', ''),
            ('GCL+_', 'GCL+ '),
            ('GCL++_', 'GCL++ '),
            ('Choroid_', ''),
            ('_厚度', ''),
            ('_体积', ' volume'),
            ('_面积', ' area'),
            ('_', ' ')
        ]
        
        result = param
        for old, new in replacements:
            result = result.replace(old, new)
        
        # 截断过长的名称
        if len(result) > 25:
            result = result[:22] + '...'
        
        return result
    
    y_labels = [simplify_param_name(p) for p in df['Parameter']]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels)
    
    # 设置x轴
    ax.set_xlabel("Cohen's d (95% CI)", fontsize=9, fontweight='normal')
    ax.set_xlim(-0.7, 0.1)  # 根据数据范围调整
    
    # 添加网格
    ax.grid(True, axis='x', linestyle=':', linewidth=0.3, alpha=0.5)
    
    # 添加类别图例
    categories = df['Category'].unique()
    legend_elements = []
    for cat in categories:
        if cat == 'Retina':
            color = colors['retina']
        elif cat == 'RNFL':
            color = colors['rnfl']
        elif cat == 'GCL':
            color = colors['gcl']
        elif cat == 'Optic Disc':
            color = colors['disc']
        else:
            color = colors['choroid']
        
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                         markerfacecolor=color, 
                                         markeredgecolor='black',
                                         markeredgewidth=0.5,
                                         markersize=8, label=cat))
    
    ax.legend(handles=legend_elements, loc='lower left', frameon=True, 
              framealpha=0.9, edgecolor='black', fontsize=7)
    
    # 添加标题
    ax.set_title('Group Comparison of Key OCT Parameters\n(MDD vs. Control, n=463 eyes)', 
                fontsize=10, fontweight='bold', pad=15)
    
    # 添加样本量信息
    sample_info = f"MDD: {df['MDD_n'].iloc[0]} eyes, Control: {df['Control_n'].iloc[0]} eyes"
    ax.text(0.02, 0.98, sample_info, transform=ax.transAxes,
           fontsize=7, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 添加显著性说明
    sig_text = "*P<0.05, **P<0.01, ***P<0.001"
    ax.text(0.98, 0.02, sig_text, transform=ax.transAxes,
           fontsize=6, horizontalalignment='right',
           verticalalignment='bottom')
    
    # 调整布局
    plt.tight_layout()
    
    return fig

def main():
    """主函数"""
    print("制作Figure 2: Group comparison forest plot...")
    
    # 加载数据
    df = load_group_comparison_data()
    
    # 显示关键信息
    print(f"\n展示前{len(df)}个效应量最大的指标:")
    for idx, row in df.iterrows():
        print(f"  {row['Parameter']}: d={row['Cohens_d']:.3f} ({row['Cohens_d_CI_lower']:.3f} to {row['Cohens_d_CI_upper']:.3f}), P={row['P_value']:.6f}")
    
    # 创建图表
    fig = create_forest_plot(df)
    
    # 保存图表
    output_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/02_Figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存为高分辨率TIFF
    output_path = os.path.join(output_dir, "Figure2_Group_comparison_forest_plot.tiff")
    fig.savefig(output_path, format='tiff', dpi=300)
    print(f"\n图表已保存: {output_path}")
    
    # 保存为PDF格式（矢量图）
    pdf_path = output_path.replace('.tiff', '.pdf')
    fig.savefig(pdf_path, format='pdf')
    print(f"矢量图已保存: {pdf_path}")
    
    # 保存为PNG格式（预览用）
    png_path = output_path.replace('.tiff', '.png')
    fig.savefig(png_path, format='png', dpi=150)
    print(f"预览图已保存: {png_path}")
    
    # 显示图表信息
    print(f"\n图表规格:")
    print(f"  尺寸: {fig.get_size_inches()[0]:.2f} × {fig.get_size_inches()[1]:.2f} inches")
    print(f"  分辨率: 300 DPI")
    print(f"  格式: TIFF, PDF, PNG")
    print(f"  包含: {len(df)} 个OCT指标")
    print(f"  最大效应量: {df['Cohens_d'].min():.3f} ({df.loc[df['Cohens_d'].idxmin(), 'Parameter']})")
    
    plt.close(fig)
    
    # 生成图表说明
    caption = f"""**Figure 2. Forest plot of group comparison for key OCT parameters.** 
Cohen's d with 95% confidence intervals are shown for the 15 OCT parameters with largest effect sizes comparing major depressive disorder (MDD) patients (n={df['MDD_n'].iloc[0]} eyes) with healthy controls (n={df['Control_n'].iloc[0]} eyes). 
Negative values indicate lower values in MDD patients. Significance levels: *P<0.05, **P<0.01, ***P<0.001. 
The largest effect was observed for outer temporal retinal thickness (Cohen's d=-0.50, P<0.001)."""
    
    caption_path = os.path.join(output_dir, "Figure2_caption.txt")
    with open(caption_path, 'w', encoding='utf-8') as f:
        f.write(caption)
    print(f"\n图表说明已保存: {caption_path}")

if __name__ == "__main__":
    main()