#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
制作SCI图表 - Figure 6: Subgroup analysis forest plot
基于463眼OCT-MDD研究数据
展示性别分层分析结果
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
    'font.family': 'DejaVu Sans',  # 使用可用字体
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
    'overall': '#999999',  # 灰色 - 整体
    'male': '#0072B2',  # 蓝色 - 男性
    'female': '#CC79A7',  # 粉色 - 女性
    'significant': '#D55E00',  # 橙色 - 显著
    'ns': '#999999',  # 灰色 - 不显著
    'retina': '#E69F00',  # 黄色 - 视网膜参数
    'reference': '#000000'  # 黑色 - 参考线
}

def load_subgroup_data():
    """加载亚组分析数据"""
    data_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/03_Tables/亚组分析结果_20260315.xlsx"
    df = pd.read_excel(data_path)
    
    print(f"亚组分析数据: {len(df)} 行")
    print(f"亚组类型: {df['亚组'].unique()}")
    print(f"指标数量: {len(df['指标'].unique())}")
    
    # 重命名列以便处理
    df = df.rename(columns={
        '亚组': 'Subgroup',
        '指标': 'Parameter',
        '样本量': 'Sample_size',
        '系数β': 'Beta',
        'P值': 'P_value',
        '95%CI下限': 'CI_lower',
        '95%CI上限': 'CI_upper'
    })
    
    # 检查数据完整性
    print(f"\n数据完整性检查:")
    print(f"  总行数: {len(df)}")
    print(f"  缺失Beta值: {df['Beta'].isna().sum()}")
    print(f"  缺失P值: {df['P_value'].isna().sum()}")
    
    # 过滤数据：只保留整体、男性、女性亚组，并移除缺失值
    sex_subgroups = ['整体', '男性', '女性']
    df_filtered = df[df['Subgroup'].isin(sex_subgroups)].copy()
    df_filtered = df_filtered.dropna(subset=['Beta', 'P_value', 'CI_lower', 'CI_upper'])
    
    print(f"\n过滤后数据:")
    print(f"  保留亚组: {sex_subgroups}")
    print(f"  过滤后行数: {len(df_filtered)}")
    print(f"  每个亚组应有行数: {len(df_filtered['Parameter'].unique())} × 3 = {len(df_filtered['Parameter'].unique()) * 3}")
    
    # 按指标和亚组排序
    # 定义指标顺序（按整体效应量大小）
    overall_data = df_filtered[df_filtered['Subgroup'] == '整体'].copy()
    overall_data = overall_data.sort_values('Beta')
    
    parameter_order = overall_data['Parameter'].tolist()
    print(f"\n指标排序 (按整体效应量):")
    for i, param in enumerate(parameter_order):
        beta = overall_data[overall_data['Parameter'] == param]['Beta'].values[0]
        print(f"  {i+1}. {param}: β={beta:.3f}")
    
    # 创建排序映射
    param_to_order = {param: i for i, param in enumerate(parameter_order)}
    df_filtered['Parameter_order'] = df_filtered['Parameter'].map(param_to_order)
    
    # 按指标和亚组排序
    df_filtered = df_filtered.sort_values(['Parameter_order', 'Subgroup']).reset_index(drop=True)
    
    # 添加显著性标记
    df_filtered['Significance'] = 'ns'
    df_filtered.loc[df_filtered['P_value'] < 0.05, 'Significance'] = '*'
    df_filtered.loc[df_filtered['P_value'] < 0.01, 'Significance'] = '**'
    df_filtered.loc[df_filtered['P_value'] < 0.001, 'Significance'] = '***'
    
    # 添加样本量信息
    print(f"\n样本量分布:")
    for subgroup in df_filtered['Subgroup'].unique():
        subgroup_data = df_filtered[df_filtered['Subgroup'] == subgroup]
        sample_size = subgroup_data['Sample_size'].iloc[0]
        n_params = len(subgroup_data)
        print(f"  {subgroup}: n={sample_size}, {n_params}个指标")
    
    return df_filtered, parameter_order

def simplify_parameter_name(param):
    """简化参数名称用于显示"""
    # 定义映射
    mappings = {
        'Macular_Outer_Temporal_Thickness': 'Outer Temporal Thickness',
        'Macular_Inner_Temporal_Thickness': 'Inner Temporal Thickness',
        'Macular_Outer_Superior_Thickness': 'Outer Superior Thickness',
        'Mean_Macular_Thickness': 'Mean Macular Thickness',
        'Total_Macular_Volume': 'Total Macular Volume'
    }
    
    return mappings.get(param, param)

def get_subgroup_color(subgroup):
    """获取亚组对应的颜色"""
    if subgroup == '整体':
        return colors['overall']
    elif subgroup == '男性':
        return colors['male']
    elif subgroup == '女性':
        return colors['female']
    else:
        return colors['ns']

def create_forest_plot(df, parameter_order):
    """创建森林图"""
    # 创建图形 - 单栏宽度 (3.35 inches)
    fig, ax = plt.subplots(figsize=(3.35, 4.5))
    
    # 准备数据
    # 每个指标有3个亚组（整体、男性、女性）
    n_params = len(parameter_order)
    n_subgroups = 3  # 整体、男性、女性
    
    # 计算y轴位置
    # 每个指标组之间有额外间距
    y_positions = []
    param_labels = []
    
    for i, param in enumerate(parameter_order):
        # 添加指标标签位置
        y_positions.append(i * (n_subgroups + 1) + 1.5)  # 指标标签位置
        param_labels.append(simplify_parameter_name(param))
        
        # 添加三个亚组数据点位置
        for j, subgroup in enumerate(['整体', '男性', '女性']):
            y_pos = i * (n_subgroups + 1) + j * 0.8
            y_positions.append(y_pos)
    
    # 调整y_positions顺序
    y_positions_sorted = []
    subgroup_labels = []
    
    for i, param in enumerate(parameter_order):
        # 首先添加三个亚组数据点
        for j, subgroup in enumerate(['整体', '男性', '女性']):
            y_pos = i * (n_subgroups + 1) + j * 0.8
            y_positions_sorted.append(y_pos)
            subgroup_labels.append(subgroup)
        
        # 然后添加指标标签位置（在亚组上方）
        y_label_pos = i * (n_subgroups + 1) + 2.4
        y_positions_sorted.append(y_label_pos)
        subgroup_labels.append('label')
    
    # 绘制数据
    point_idx = 0
    
    for i, param in enumerate(parameter_order):
        param_data = df[df['Parameter'] == param]
        
        for subgroup in ['整体', '男性', '女性']:
            subgroup_data = param_data[param_data['Subgroup'] == subgroup]
            
            if len(subgroup_data) > 0:
                row = subgroup_data.iloc[0]
                beta = row['Beta']
                ci_lower = row['CI_lower']
                ci_upper = row['CI_upper']
                p_value = row['P_value']
                significance = row['Significance']
                sample_size = row['Sample_size']
                
                # 获取颜色
                color = get_subgroup_color(subgroup)
                
                # 绘制效应量点
                ax.plot(beta, y_positions_sorted[point_idx], 'o',
                       markersize=8 if subgroup == '整体' else 7,
                       color=color,
                       markeredgecolor='black',
                       markeredgewidth=0.5,
                       zorder=3)
                
                # 绘制置信区间线
                ax.plot([ci_lower, ci_upper], 
                       [y_positions_sorted[point_idx], y_positions_sorted[point_idx]],
                       color=color, linewidth=1.5 if subgroup == '整体' else 1.2, 
                       zorder=2)
                
                # 绘制置信区间端帽
                cap_size = 0.08
                ax.plot([ci_lower, ci_lower], 
                       [y_positions_sorted[point_idx] - cap_size, y_positions_sorted[point_idx] + cap_size],
                       color=color, linewidth=1.5 if subgroup == '整体' else 1.2, 
                       zorder=2)
                ax.plot([ci_upper, ci_upper], 
                       [y_positions_sorted[point_idx] - cap_size, y_positions_sorted[point_idx] + cap_size],
                       color=color, linewidth=1.5 if subgroup == '整体' else 1.2, 
                       zorder=2)
                
                # 添加显著性标记
                if significance != 'ns':
                    ax.text(beta + 0.3, y_positions_sorted[point_idx], significance,
                           fontsize=7, fontweight='bold', va='center', ha='left')
                
                # 添加样本量标签（在右侧）
                sample_text = f"n={sample_size}"
                ax.text(8.0, y_positions_sorted[point_idx], sample_text,
                       fontsize=6, va='center', ha='right')
                
                point_idx += 1
        
        # 跳过指标标签位置
        point_idx += 1
    
    # 添加零效应参考线
    ax.axvline(x=0, color=colors['reference'], linestyle='--', 
               linewidth=0.8, alpha=0.5, zorder=1)
    
    # 设置y轴标签
    # 只显示指标名称（在亚组数据上方）
    y_tick_positions = []
    y_tick_labels = []
    
    for i, param in enumerate(parameter_order):
        y_label_pos = i * (n_subgroups + 1) + 2.4
        y_tick_positions.append(y_label_pos)
        y_tick_labels.append(simplify_parameter_name(param))
    
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(y_tick_labels)
    
    # 设置x轴
    ax.set_xlabel('Coefficient β (95% CI)', fontsize=9, fontweight='normal')
    
    # 根据数据范围设置x轴限制
    all_betas = df['Beta'].values
    all_cis = np.concatenate([df['CI_lower'].values, df['CI_upper'].values])
    
    # 移除NaN和Inf值
    all_betas_clean = all_betas[~np.isnan(all_betas) & ~np.isinf(all_betas)]
    all_cis_clean = all_cis[~np.isnan(all_cis) & ~np.isinf(all_cis)]
    
    if len(all_betas_clean) > 0 and len(all_cis_clean) > 0:
        x_min = min(all_cis_clean.min(), -1) - 2
        x_max = max(all_cis_clean.max(), 1) + 2
        ax.set_xlim(x_min, x_max)
    else:
        # 默认范围
        ax.set_xlim(-15, 5)
    
    # 添加网格
    ax.grid(True, axis='x', linestyle=':', linewidth=0.3, alpha=0.5)
    
    # 添加图例
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=colors['overall'], markeredgecolor='black',
                  markeredgewidth=0.5, markersize=8, label='Overall'),
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=colors['male'], markeredgecolor='black',
                  markeredgewidth=0.5, markersize=7, label='Male'),
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=colors['female'], markeredgecolor='black',
                  markeredgewidth=0.5, markersize=7, label='Female')
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', frameon=True,
              framealpha=0.9, edgecolor='black', fontsize=7)
    
    # 添加标题
    ax.set_title('Subgroup Analysis by Sex\n(Linear regression coefficients)', 
                fontsize=10, fontweight='bold', pad=15)
    
    # 添加样本总量信息
    total_samples = df[df['Subgroup'] == '整体']['Sample_size'].sum() / len(parameter_order)
    sample_info = f"Total: n={int(total_samples)} eyes (Male: 126, Female: 337)"
    ax.text(0.02, 0.98, sample_info, transform=ax.transAxes,
           fontsize=7, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 添加显著性说明
    sig_text = "*P<0.05, **P<0.01, ***P<0.001"
    ax.text(0.98, 0.02, sig_text, transform=ax.transAxes,
           fontsize=6, horizontalalignment='right',
           verticalalignment='bottom')
    
    # 添加异质性说明（基于P值差异）
    # 计算性别间差异
    print(f"\n性别差异分析:")
    for param in parameter_order:
        male_data = df[(df['Parameter'] == param) & (df['Subgroup'] == '男性')]
        female_data = df[(df['Parameter'] == param) & (df['Subgroup'] == '女性')]
        
        if len(male_data) > 0 and len(female_data) > 0:
            male_beta = male_data.iloc[0]['Beta']
            female_beta = female_data.iloc[0]['Beta']
            diff = abs(female_beta - male_beta)
            
            # 简单判断：差异大于整体效应的20%视为可能有意义
            overall_data = df[(df['Parameter'] == param) & (df['Subgroup'] == '整体')]
            if len(overall_data) > 0:
                overall_beta = overall_data.iloc[0]['Beta']
                if overall_beta != 0:
                    rel_diff = diff / abs(overall_beta)
                    if rel_diff > 0.2:
                        print(f"  {param}: 性别差异较大 (β差值={diff:.2f}, 相对差异={rel_diff:.1%})")
    
    # 调整布局
    plt.tight_layout()
    
    return fig

def create_simple_forest_plot(df, parameter_order):
    """创建简化版森林图（所有亚组在同一水平线）"""
    fig, ax = plt.subplots(figsize=(3.35, 3.5))
    
    # 每个指标一行，三个亚组并排
    y_positions = np.arange(len(parameter_order))
    
    # 简化指标名称
    param_labels = [simplify_parameter_name(p) for p in parameter_order]
    
    # 设置亚组偏移量
    subgroup_offsets = {'整体': 0, '男性': -0.2, '女性': 0.2}
    
    # 绘制每个指标
    for i, param in enumerate(parameter_order):
        param_data = df[df['Parameter'] == param]
        
        for subgroup in ['整体', '男性', '女性']:
            subgroup_data = param_data[param_data['Subgroup'] == subgroup]
            
            if len(subgroup_data) > 0:
                row = subgroup_data.iloc[0]
                beta = row['Beta']
                ci_lower = row['CI_lower']
                ci_upper = row['CI_upper']
                color = get_subgroup_color(subgroup)
                
                # 计算x位置（考虑偏移）
                x_pos = beta
                y_pos = y_positions[i] + subgroup_offsets[subgroup]
                
                # 绘制点
                ax.plot(x_pos, y_pos, 'o',
                       markersize=7,
                       color=color,
                       markeredgecolor='black',
                       markeredgewidth=0.5,
                       zorder=3,
                       label=subgroup if i == 0 else "")
                
                # 绘制置信区间
                ax.plot([ci_lower, ci_upper], [y_pos, y_pos],
                       color=color, linewidth=1.2, zorder=2)
                
                # 端帽
                cap_size = 0.03
                ax.plot([ci_lower, ci_lower], [y_pos - cap_size, y_pos + cap_size],
                       color=color, linewidth=1.2, zorder=2)
                ax.plot([ci_upper, ci_upper], [y_pos - cap_size, y_pos + cap_size],
                       color=color, linewidth=1.2, zorder=2)
    
    # 设置y轴
    ax.set_yticks(y_positions)
    ax.set_yticklabels(param_labels)
    ax.invert_yaxis()  # 最重要的在顶部
    
    # 设置x轴
    ax.set_xlabel('Coefficient β (95% CI)', fontsize=9)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # 添加图例
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # 去重
    ax.legend(by_label.values(), by_label.keys(), loc='upper right',
              frameon=True, framealpha=0.9, fontsize=7)
    
    # 添加标题
    ax.set_title('Subgroup Analysis by Sex', fontsize=10, fontweight='bold', pad=15)
    
    # 添加样本量信息
    sample_info = "Overall: n=463, Male: n=126, Female: n=337"
    ax.text(0.02, 0.98, sample_info, transform=ax.transAxes,
           fontsize=7, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    return fig

def main():
    """主函数"""
    print("制作Figure 6: Subgroup analysis forest plot...")
    
    # 加载数据
    df, parameter_order = load_subgroup_data()
    
    # 显示关键信息
    print(f"\n关键发现:")
    for param in parameter_order:
        param_data = df[df['Parameter'] == param]
        print(f"\n{simplify_parameter_name(param)}:")
        
        for subgroup in ['整体', '男性', '女性']:
            subgroup_data = param_data[param_data['Subgroup'] == subgroup]
            if len(subgroup_data) > 0:
                row = subgroup_data.iloc[0]
                print(f"  {subgroup}: β={row['Beta']:.3f} ({row['CI_lower']:.3f} to {row['CI_upper']:.3f}), P={row['P_value']:.6f}")
    
    # 创建详细森林图
    fig_detailed = create_forest_plot(df, parameter_order)
    
    # 创建简化版森林图（备选）
    fig_simple = create_simple_forest_plot(df, parameter_order)
    
    # 保存图表
    output_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/02_Figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存详细版
    output_path = os.path.join(output_dir, "Figure6_Subgroup_analysis_forest_plot.tiff")
    fig_detailed.savefig(output_path, format='tiff', dpi=300)
    print(f"\n详细森林图已保存: {output_path}")
    
    pdf_path = output_path.replace('.tiff', '.pdf')
    fig_detailed.savefig(pdf_path, format='pdf')
    print(f"矢量图已保存: {pdf_path}")
    
    png_path = output_path.replace('.tiff', '.png')
    fig_detailed.savefig(png_path, format='png', dpi=150)
    print(f"预览图已保存: {png_path}")
    
    # 保存简化版（备选）
    simple_path = os.path.join(output_dir, "Figure6_Subgroup_analysis_simple.tiff")
    fig_simple.savefig(simple_path, format='tiff', dpi=300)
    print(f"简化森林图已保存: {simple_path}")
    
    # 显示图表信息
    print(f"\n图表规格:")
    print(f"  详细版尺寸: {fig_detailed.get_size_inches()[0]:.2f} × {fig_detailed.get_size_inches()[1]:.2f} inches")
    print(f"  简化版尺寸: {fig_simple.get_size_inches()[0]:.2f} × {fig_simple.get_size_inches()[1]:.2f} inches")
    print(f"  分辨率: 300 DPI")
    print(f"  格式: TIFF, PDF, PNG")
    
    # 总结亚组分析发现
    print(f"\n亚组分析核心发现:")
    
    # 找出性别差异最大的指标
    max_diff = 0
    max_diff_param = None
    
    for param in parameter_order:
        male_data = df[(df['Parameter'] == param) & (df['Subgroup'] == '男性')]
        female_data = df[(df['Parameter'] == param) & (df['Subgroup'] == '女性')]
        
        if len(male_data) > 0 and len(female_data) > 0:
            male_beta = male_data.iloc[0]['Beta']
            female_beta = female_data.iloc[0]['Beta']
            diff = abs(female_beta - male_beta)
            
            if diff > max_diff:
                max_diff = diff
                max_diff_param = param
    
    if max_diff_param:
        print(f"  性别差异最大的指标: {simplify_parameter_name(max_diff_param)} (β差值={max_diff:.3f})")
    
    # 检查整体显著性
    overall_sig_params = df[(df['Subgroup'] == '整体') & (df['P_value'] < 0.05)]['Parameter']
    print(f"  整体显著的指标数: {len(overall_sig_params)}/{len(parameter_order)}")
    
    # 检查男性组显著性
    male_sig_params = df[(df['Subgroup'] == '男性') & (df['P_value'] < 0.05)]['Parameter']
    print(f"  男性组显著的指标数: {len(male_sig_params)}/{len(parameter_order)}")
    
    # 检查女性组显著性
    female_sig_params = df[(df['Subgroup'] == '女性') & (df['P_value'] < 0.05)]['Parameter']
    print(f"  女性组显著的指标数: {len(female_sig_params)}/{len(parameter_order)}")
    
    plt.close('all')
    
    # 生成图表说明
    caption = f"""**Figure 6. Forest plot of subgroup analysis by sex.** 
Linear regression coefficients (β) with 95% confidence intervals are shown for five key macular OCT parameters, stratified by sex. 
The overall analysis included 463 eyes (126 males, 337 females). 
Outer temporal macular thickness showed the strongest association with depression status in the overall sample (β=-6.325, P=0.000139), 
with a more pronounced effect in females (β=-6.946, P=0.000589) than in males (β=-5.274, P=0.074). 
All five parameters showed statistically significant associations in the overall sample (P<0.05), 
with four remaining significant in females but only one (total macular volume) showing a trend in males (P=0.250). 
These findings suggest that retinal thinning associated with major depressive disorder may be more pronounced in females, 
potentially reflecting sex-specific neurobiological mechanisms."""
    
    caption_path = os.path.join(output_dir, "Figure6_caption.txt")
    with open(caption_path, 'w', encoding='utf-8') as f:
        f.write(caption)
    print(f"\n图表说明已保存: {caption_path}")

if __name__ == "__main__":
    main()