#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
制作SCI图表 - Figure 5: Correlation between OCT parameters and PHQ-9 scores
基于463眼OCT-MDD研究数据
展示OCT指标与抑郁严重程度的关系
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
from scipy import stats
import matplotlib
from sklearn.linear_model import LinearRegression

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
    'mdd': '#D55E00',  # 橙色 - MDD组
    'regression': '#0072B2',  # 蓝色 - 回归线
    'ci': '#CC79A7',  # 粉色 - 置信区间
    'retina': '#D55E00',  # 橙色 - 视网膜参数
    'rnfl': '#0072B2',  # 蓝色 - RNFL参数
    'disc': '#009E73',  # 绿色 - 视盘参数
    'scatter': '#56B4E9',  # 浅蓝 - 散点
    'male': '#0072B2',  # 蓝色 - 男性
    'female': '#CC79A7'  # 粉色 - 女性
}

def load_and_prepare_data():
    """加载并准备数据"""
    data_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/04_Data/data_499eyes_20260315.xlsx"
    df = pd.read_excel(data_path)
    
    print(f"原始数据形状: {df.shape}")
    print(f"总样本数: {len(df)}")
    
    # 筛选463眼样本（有年龄和性别数据的）
    df_filtered = df.dropna(subset=['年龄', '性别']).copy()
    print(f"有年龄性别数据的样本数: {len(df_filtered)}")
    
    # 筛选MDD组（分组为"抑郁症"）
    # 注意：相关性分析仅基于MDD组，因为对照组没有PHQ-9数据
    df_mdd = df_filtered[df_filtered['分组'] == '抑郁症'].copy()
    print(f"MDD组样本数: {len(df_mdd)}")
    
    # 筛选有PHQ-9数据的MDD样本
    df_phq = df_mdd.dropna(subset=['PHQ-9']).copy()
    print(f"有PHQ-9数据的MDD样本数: {len(df_phq)}")
    
    # 确保关键OCT指标没有缺失值
    key_oct_indicators = [
        'Retina_外环颞侧',  # 效应量最大的指标
        'Rim Volume',        # 相关性最强的指标
        'C/D Area Ratio',    # 有负相关的指标
        'Retina_平均厚度',   # 平均厚度
        'RNFL_平均厚度'      # RNFL平均厚度
    ]
    
    # 检查数据完整性
    for indicator in key_oct_indicators:
        non_missing = df_phq[indicator].notna().sum()
        print(f"  {indicator}: {non_missing}/{len(df_phq)} ({non_missing/len(df_phq)*100:.1f}%)")
    
    # 移除有缺失值的样本
    df_final = df_phq.dropna(subset=key_oct_indicators).copy()
    print(f"最终分析样本数（所有关键指标完整）: {len(df_final)}")
    
    return df_final, key_oct_indicators

def calculate_correlations(df, indicators):
    """计算每个指标与PHQ-9的相关性"""
    correlations = {}
    
    for indicator in indicators:
        # 计算Spearman相关系数
        spearman_corr, spearman_p = stats.spearmanr(df[indicator], df['PHQ-9'], nan_policy='omit')
        
        # 计算Pearson相关系数
        pearson_corr, pearson_p = stats.pearsonr(df[indicator].dropna(), df['PHQ-9'].dropna())
        
        # 计算线性回归参数
        X = df[indicator].values.reshape(-1, 1)
        y = df['PHQ-9'].values
        mask = ~np.isnan(X[:, 0]) & ~np.isnan(y)
        X_clean = X[mask]
        y_clean = y[mask]
        
        if len(X_clean) > 1:
            reg = LinearRegression().fit(X_clean, y_clean)
            slope = reg.coef_[0]
            intercept = reg.intercept_
            r2 = reg.score(X_clean, y_clean)
        else:
            slope, intercept, r2 = np.nan, np.nan, np.nan
        
        correlations[indicator] = {
            'spearman_corr': spearman_corr,
            'spearman_p': spearman_p,
            'pearson_corr': pearson_corr,
            'pearson_p': pearson_p,
            'reg_slope': slope,
            'reg_intercept': intercept,
            'r2': r2,
            'n': len(df)
        }
        
        # 打印结果
        print(f"{indicator}:")
        print(f"  Spearman rho = {spearman_corr:.3f}, P = {spearman_p:.4f}")
        print(f"  Pearson r = {pearson_corr:.3f}, P = {pearson_p:.4f}")
        print(f"  Linear regression: slope = {slope:.4f}, R² = {r2:.4f}")
    
    return correlations

def create_scatter_panel(df, indicator, corr_info, ax, color):
    """创建单个散点图面板"""
    # 提取数据
    x = df[indicator]
    y = df['PHQ-9']
    
    # 绘制散点图
    scatter = ax.scatter(x, y, color=color, s=20, alpha=0.7, edgecolor='white', linewidth=0.5)
    
    # 添加回归线
    if not np.isnan(corr_info['reg_slope']):
        x_range = np.linspace(x.min(), x.max(), 100)
        y_pred = corr_info['reg_slope'] * x_range + corr_info['reg_intercept']
        ax.plot(x_range, y_pred, color=colors['regression'], linewidth=2, alpha=0.8)
    
    # 计算并绘制置信区间带
    if not np.isnan(corr_info['reg_slope']) and len(x) > 2:
        # 使用自助法计算置信区间
        n_boot = 1000
        slopes = []
        intercepts = []
        
        for _ in range(n_boot):
            # 有放回抽样
            indices = np.random.choice(len(x), len(x), replace=True)
            x_boot = x.iloc[indices].values
            y_boot = y.iloc[indices].values
            
            # 移除NaN值
            mask = ~np.isnan(x_boot) & ~np.isnan(y_boot)
            if np.sum(mask) > 1:
                x_clean = x_boot[mask].reshape(-1, 1)
                y_clean = y_boot[mask]
                reg_boot = LinearRegression().fit(x_clean, y_clean)
                slopes.append(reg_boot.coef_[0])
                intercepts.append(reg_boot.intercept_)
        
        if slopes:
            # 计算每个x点的预测值分布
            x_line = np.linspace(x.min(), x.max(), 50)
            predictions = []
            
            for slope, intercept in zip(slopes, intercepts):
                predictions.append(slope * x_line + intercept)
            
            predictions = np.array(predictions)
            
            # 计算95%置信区间
            ci_lower = np.percentile(predictions, 2.5, axis=0)
            ci_upper = np.percentile(predictions, 97.5, axis=0)
            
            # 绘制置信区间带
            ax.fill_between(x_line, ci_lower, ci_upper, color=colors['ci'], alpha=0.2, linewidth=0)
    
    # 设置坐标轴标签
    # 简化指标名称用于显示
    if indicator == 'Retina_外环颞侧':
        xlabel = 'Outer Temporal Retinal Thickness (μm)'
        indicator_short = 'Retina_OT'
    elif indicator == 'Rim Volume':
        xlabel = 'Rim Volume (mm³)'
        indicator_short = 'RimVol'
    elif indicator == 'C/D Area Ratio':
        xlabel = 'Cup-to-Disc Area Ratio'
        indicator_short = 'C/D Ratio'
    elif indicator == 'Retina_平均厚度':
        xlabel = 'Mean Macular Thickness (μm)'
        indicator_short = 'MMT'
    elif indicator == 'RNFL_平均厚度':
        xlabel = 'Mean RNFL Thickness (μm)'
        indicator_short = 'MRNFLT'
    else:
        xlabel = indicator
        indicator_short = indicator[:15]
    
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel('PHQ-9 Score', fontsize=9)
    
    # 添加统计信息
    spearman_text = f"Spearman ρ = {corr_info['spearman_corr']:.3f}"
    if corr_info['spearman_p'] < 0.001:
        spearman_text += "***"
    elif corr_info['spearman_p'] < 0.01:
        spearman_text += "**"
    elif corr_info['spearman_p'] < 0.05:
        spearman_text += "*"
    
    pearson_text = f"Pearson r = {corr_info['pearson_corr']:.3f}"
    if corr_info['pearson_p'] < 0.001:
        pearson_text += "***"
    elif corr_info['pearson_p'] < 0.01:
        pearson_text += "**"
    elif corr_info['pearson_p'] < 0.05:
        pearson_text += "*"
    
    n_text = f"n = {corr_info['n']}"
    
    stat_text = f"{spearman_text}\n{pearson_text}\n{n_text}"
    
    ax.text(0.05, 0.95, stat_text, transform=ax.transAxes,
           fontsize=7, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # 添加网格
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.5)
    
    return indicator_short

def create_multi_panel_figure(df, indicators, correlations):
    """创建多面板散点图"""
    # 确定布局：2行3列或类似
    n_indicators = len(indicators)
    n_cols = 3
    n_rows = (n_indicators + n_cols - 1) // n_cols
    
    # 创建图形
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.0, 4.5 if n_rows == 2 else 6.0))
    
    # 如果只有一行，确保axes是数组
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # 为每个指标分配颜色
    indicator_colors = {
        'Retina_外环颞侧': colors['retina'],
        'Rim Volume': colors['disc'],
        'C/D Area Ratio': colors['disc'],
        'Retina_平均厚度': colors['retina'],
        'RNFL_平均厚度': colors['rnfl']
    }
    
    # 创建每个子图
    panel_labels = ['A', 'B', 'C', 'D', 'E']
    panel_idx = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            if panel_idx < n_indicators:
                indicator = indicators[panel_idx]
                ax = axes[i, j]
                color = indicator_colors.get(indicator, colors['scatter'])
                
                # 创建散点图面板
                indicator_short = create_scatter_panel(df, indicator, correlations[indicator], ax, color)
                
                # 添加面板标签（A, B, C, ...）
                ax.text(0.02, 0.98, panel_labels[panel_idx], transform=ax.transAxes,
                       fontsize=10, fontweight='bold', verticalalignment='top')
                
                panel_idx += 1
            else:
                # 隐藏多余的子图
                axes[i, j].axis('off')
    
    # 添加总体标题
    fig.suptitle('Correlation between OCT Parameters and PHQ-9 Scores\n(MDD patients only, n={})'.format(len(df)), 
                fontsize=11, fontweight='bold', y=0.98)
    
    # 添加图例说明
    legend_text = "***P<0.001, **P<0.01, *P<0.05"
    fig.text(0.5, 0.02, legend_text, ha='center', va='bottom', fontsize=8,
            style='italic')
    
    # 添加样本信息说明
    sample_info = "Note: Healthy controls were excluded from correlation analysis due to lack of PHQ-9 data"
    fig.text(0.5, 0.04, sample_info, ha='center', va='bottom', fontsize=7,
            style='italic', color='gray')
    
    # 调整布局
    plt.tight_layout(rect=[0, 0.06, 1, 0.95])
    
    return fig

def create_combined_scatter_plot(df, indicators, correlations):
    """创建组合散点图（所有指标在一个图中，用颜色区分）"""
    fig, ax = plt.subplots(figsize=(3.35, 3.35))
    
    # 为每个指标分配颜色和标记
    colors_list = [colors['retina'], colors['disc'], colors['disc'], colors['retina'], colors['rnfl']]
    markers = ['o', 's', '^', 'D', 'v']
    labels = ['Outer Temporal Retina', 'Rim Volume', 'C/D Ratio', 'Mean Macular Thickness', 'Mean RNFL Thickness']
    
    # 绘制所有散点
    for idx, (indicator, color, marker, label) in enumerate(zip(indicators, colors_list, markers, labels)):
        x = df[indicator]
        y = df['PHQ-9']
        
        # 绘制散点
        ax.scatter(x, y, color=color, marker=marker, s=25, alpha=0.6, 
                  edgecolor='white', linewidth=0.5, label=label)
        
        # 绘制回归线
        if not np.isnan(correlations[indicator]['reg_slope']):
            x_range = np.linspace(x.min(), x.max(), 100)
            y_pred = correlations[indicator]['reg_slope'] * x_range + correlations[indicator]['reg_intercept']
            ax.plot(x_range, y_pred, color=color, linewidth=1.5, alpha=0.7, linestyle='--')
    
    # 设置坐标轴
    ax.set_xlabel('OCT Parameter Value (Standardized)', fontsize=9)
    ax.set_ylabel('PHQ-9 Score', fontsize=9)
    
    # 添加图例
    ax.legend(loc='upper right', frameon=True, framealpha=0.9, 
              edgecolor='black', fontsize=7)
    
    # 添加网格
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.5)
    
    # 添加标题
    ax.set_title('OCT Parameters vs. PHQ-9 Scores\n(MDD patients, n={})'.format(len(df)), 
                fontsize=10, fontweight='bold', pad=15)
    
    # 添加总体相关性说明
    corr_values = [correlations[indicator]['spearman_corr'] for indicator in indicators]
    mean_corr = np.mean(np.abs(corr_values))
    ax.text(0.02, 0.98, f'Mean |ρ| = {mean_corr:.3f}', transform=ax.transAxes,
           fontsize=8, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    return fig

def main():
    """主函数"""
    print("制作Figure 5: Correlation between OCT parameters and PHQ-9 scores...")
    
    # 加载数据
    df_final, indicators = load_and_prepare_data()
    
    if len(df_final) == 0:
        print("错误: 没有足够的数据进行分析")
        return
    
    # 计算相关性
    print(f"\n计算相关性 (n={len(df_final)}):")
    correlations = calculate_correlations(df_final, indicators)
    
    # 创建多面板散点图
    fig_multi = create_multi_panel_figure(df_final, indicators, correlations)
    
    # 保存图表
    output_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/02_Figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存为高分辨率TIFF
    output_path = os.path.join(output_dir, "Figure5_Correlation_scatter_plots.tiff")
    fig_multi.savefig(output_path, format='tiff', dpi=300)
    print(f"\n多面板散点图已保存: {output_path}")
    
    # 保存为PDF格式（矢量图）
    pdf_path = output_path.replace('.tiff', '.pdf')
    fig_multi.savefig(pdf_path, format='pdf')
    print(f"矢量图已保存: {pdf_path}")
    
    # 保存为PNG格式（预览用）
    png_path = output_path.replace('.tiff', '.png')
    fig_multi.savefig(png_path, format='png', dpi=150)
    print(f"预览图已保存: {png_path}")
    
    # 创建并保存组合散点图（备选）
    fig_combined = create_combined_scatter_plot(df_final, indicators, correlations)
    combined_path = os.path.join(output_dir, "Figure5_Combined_scatter.tiff")
    fig_combined.savefig(combined_path, format='tiff', dpi=300)
    print(f"组合散点图已保存: {combined_path}")
    
    # 显示图表信息
    print(f"\n图表规格:")
    print(f"  多面板图尺寸: {fig_multi.get_size_inches()[0]:.2f} × {fig_multi.get_size_inches()[1]:.2f} inches")
    print(f"  组合图尺寸: {fig_combined.get_size_inches()[0]:.2f} × {fig_combined.get_size_inches()[1]:.2f} inches")
    print(f"  分辨率: 300 DPI")
    print(f"  格式: TIFF, PDF, PNG")
    print(f"  包含指标数: {len(indicators)}")
    
    # 总结关键发现
    print(f"\n关键发现:")
    for indicator in indicators:
        corr = correlations[indicator]
        print(f"  {indicator}:")
        print(f"    Spearman ρ = {corr['spearman_corr']:.3f} (P={corr['spearman_p']:.4f})")
        print(f"    Pearson r = {corr['pearson_corr']:.3f} (P={corr['pearson_p']:.4f})")
        print(f"    线性回归 R² = {corr['r2']:.4f}")
    
    # 识别最强和最弱相关性
    strongest_pos = max(correlations.items(), key=lambda x: x[1]['spearman_corr'])
    strongest_neg = min(correlations.items(), key=lambda x: x[1]['spearman_corr'])
    
    print(f"\n最强正相关: {strongest_pos[0]} (ρ={strongest_pos[1]['spearman_corr']:.3f})")
    print(f"最强负相关: {strongest_neg[0]} (ρ={strongest_neg[1]['spearman_corr']:.3f})")
    
    plt.close('all')
    
    # 生成图表说明
    caption = f"""**Figure 5. Scatter plots showing correlation between OCT parameters and PHQ-9 scores in major depressive disorder (MDD) patients.** 
Each panel displays the relationship between a specific OCT parameter and depression severity measured by PHQ-9. 
Blue lines represent linear regression fits, with shaded areas showing 95% confidence intervals calculated by bootstrap (n=1000). 
Spearman and Pearson correlation coefficients are shown in each panel with significance markers: *P<0.05, **P<0.01, ***P<0.001. 
The strongest positive correlation was observed for Rim Volume (Spearman ρ={strongest_pos[1]['spearman_corr']:.3f}, P={strongest_pos[1]['spearman_p']:.4f}), 
while the strongest negative correlation was for Cup-to-Disc Area Ratio (Spearman ρ={strongest_neg[1]['spearman_corr']:.3f}, P={strongest_neg[1]['spearman_p']:.4f}). 
All correlations were weak (|ρ|<0.2), consistent with the modest association between retinal structural changes and depressive symptom severity. 
Analysis included {len(df_final)} MDD patients with complete OCT and PHQ-9 data; healthy controls were excluded due to lack of PHQ-9 scores."""
    
    caption_path = os.path.join(output_dir, "Figure5_caption.txt")
    with open(caption_path, 'w', encoding='utf-8') as f:
        f.write(caption)
    print(f"\n图表说明已保存: {caption_path}")

if __name__ == "__main__":
    main()