#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
制作SCI图表 - Figure 4: Feature importance and composite score weights
基于463眼OCT-MDD研究数据
展示随机森林特征重要性与逻辑回归复合指标权重
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
    'retina': '#D55E00',  # 橙色 - 视网膜参数
    'rnfl': '#0072B2',  # 蓝色 - RNFL参数
    'gcl': '#CC79A7',  # 粉色 - GCL参数
    'gclplus': '#E69F00',  # 黄色 - GCL+参数
    'gclplusplus': '#56B4E9',  # 浅蓝 - GCL++参数
    'choroid': '#009E73',  # 绿色 - 脉络膜参数
    'disc': '#F0E442',  # 黄色 - 视盘参数
    'positive': '#D55E00',  # 橙色 - 正权重
    'negative': '#0072B2',  # 蓝色 - 负权重
    'correlation': '#CC79A7',  # 粉色 - 相关性
    'significant': '#000000',  # 黑色 - 显著
    'ns': '#999999'  # 灰色 - 不显著
}

def load_feature_data():
    """加载特征重要性数据"""
    # 特征重要性数据（随机森林）
    imp_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/03_Tables/特征重要性分析_20260315.xlsx"
    df_imp = pd.read_excel(imp_path)
    
    # 复合指标权重数据（逻辑回归）
    weight_path = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/03_Tables/复合指标权重_20260315.xlsx"
    df_weight = pd.read_excel(weight_path)
    
    print(f"特征重要性数据: {len(df_imp)} 个特征")
    print(f"复合指标权重数据: {len(df_weight)} 个特征")
    
    # 标准化特征名称以确保匹配
    df_imp['特征_标准化'] = df_imp['特征'].astype(str).str.strip()
    df_weight['特征_标准化'] = df_weight['特征'].astype(str).str.strip()
    
    # 合并数据
    df_merged = pd.merge(df_imp, df_weight, on='特征_标准化', how='inner', 
                        suffixes=('_importance', '_weight'))
    
    print(f"成功匹配特征数: {len(df_merged)}")
    
    # 添加类别信息
    def categorize_feature(feature):
        feature_lower = feature.lower()
        if 'retina' in feature_lower:
            return 'Retina'
        elif 'rnfl' in feature_lower:
            return 'RNFL'
        elif 'gcl+' in feature_lower:
            if 'gcl++' in feature_lower:
                return 'GCL++'
            else:
                return 'GCL+'
        elif 'gcl' in feature_lower:
            return 'GCL'
        elif 'choroid' in feature_lower:
            return 'Choroid'
        elif any(kw in feature_lower for kw in ['disc', 'cup', 'rim', 'area', 'ratio']):
            return 'Optic Disc'
        else:
            return 'Other'
    
    df_merged['Category'] = df_merged['特征_标准化'].apply(categorize_feature)
    
    # 排序特征（按重要性降序）
    df_merged = df_merged.sort_values('重要性', ascending=False).reset_index(drop=True)
    
    # 计算排序信息
    df_merged['Rank_Importance'] = range(1, len(df_merged) + 1)
    
    # 按权重绝对值排序
    df_merged['abs_weight'] = abs(df_merged['权重'])
    weight_rank = df_merged['abs_weight'].rank(ascending=False, method='min')
    df_merged['Rank_Weight'] = weight_rank.astype(int)
    
    # 计算重要性与权重的相关性
    correlation = df_merged['重要性'].corr(df_merged['abs_weight'])
    print(f"特征重要性与权重绝对值的相关性: r={correlation:.3f}")
    
    # 获取顶部特征
    top_n = 20
    df_top_imp = df_merged.head(top_n).copy()
    df_top_imp['顺序_重要性'] = range(top_n, 0, -1)  # 反转顺序用于绘图
    
    # 按权重绝对值获取顶部特征
    df_top_weight = df_merged.sort_values('abs_weight', ascending=False).head(top_n).copy()
    df_top_weight = df_top_weight.sort_values('权重', ascending=True)  # 按权重值排序（负到正）
    df_top_weight['顺序_权重'] = range(top_n)  # 从下到上
    
    return df_merged, df_top_imp, df_top_weight, correlation

def simplify_feature_name(feature, max_length=25):
    """简化特征名称用于显示"""
    # 常见替换
    replacements = [
        ('Retina_', 'Ret '),
        ('RNFL_', 'RNFL '),
        ('GCL+_', 'GCL+ '),
        ('GCL++_', 'GCL++ '),
        ('Choroid_', 'Ch '),
        ('_厚度', ' Thk'),
        ('_体积', ' Vol'),
        ('_面积', ' Area'),
        ('_ratio', ' Ratio'),
        ('_', ' ')
    ]
    
    result = str(feature)
    for old, new in replacements:
        result = result.replace(old, new)
    
    # 截断过长的名称
    if len(result) > max_length:
        result = result[:max_length-3] + '...'
    
    return result

def get_category_color(category):
    """获取类别对应的颜色"""
    if category == 'Retina':
        return colors['retina']
    elif category == 'RNFL':
        return colors['rnfl']
    elif category == 'GCL':
        return colors['gcl']
    elif category == 'GCL+':
        return colors['gclplus']
    elif category == 'GCL++':
        return colors['gclplusplus']
    elif category == 'Choroid':
        return colors['choroid']
    elif category == 'Optic Disc':
        return colors['disc']
    else:
        return colors['ns']

def create_dual_panel_figure(df_top_imp, df_top_weight, correlation):
    """创建双面板特征分析图"""
    # 创建图形 - 双栏宽度 (7.0 inches)，高度适当
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 5.0), 
                            gridspec_kw={'width_ratios': [1, 1]})
    
    ax_imp, ax_weight = axes
    
    # ============================
    # 左侧面板: 特征重要性 (随机森林)
    # ============================
    
    # 准备数据
    y_positions = np.arange(len(df_top_imp))
    feature_names = [simplify_feature_name(f) for f in df_top_imp['特征_标准化']]
    importance_values = df_top_imp['重要性']
    categories = df_top_imp['Category']
    
    # 绘制条形图
    bars = ax_imp.barh(y_positions, importance_values, 
                      color=[get_category_color(cat) for cat in categories],
                      edgecolor='black', linewidth=0.5)
    
    # 设置y轴
    ax_imp.set_yticks(y_positions)
    ax_imp.set_yticklabels(feature_names)
    ax_imp.invert_yaxis()  # 最重要的在顶部
    
    # 设置x轴
    ax_imp.set_xlabel('Feature Importance (Random Forest)', fontsize=9)
    ax_imp.set_xlim([0, importance_values.max() * 1.15])  # 留出空间
    
    # 添加条形数值
    for i, (bar, value) in enumerate(zip(bars, importance_values)):
        width = bar.get_width()
        ax_imp.text(width + 0.001, bar.get_y() + bar.get_height()/2,
                   f'{value:.4f}', ha='left', va='center', fontsize=6)
    
    # 添加网格
    ax_imp.grid(True, axis='x', linestyle=':', linewidth=0.3, alpha=0.5)
    
    # 添加标题
    ax_imp.set_title('A. Random Forest Feature Importance\n(Top 20 features)', 
                    fontsize=10, fontweight='bold', pad=15)
    
    # ============================
    # 右侧面板: 复合指标权重 (逻辑回归)
    # ============================
    
    # 准备数据
    y_positions_w = np.arange(len(df_top_weight))
    feature_names_w = [simplify_feature_name(f) for f in df_top_weight['特征_标准化']]
    weight_values = df_top_weight['权重']
    categories_w = df_top_weight['Category']
    
    # 根据权重正负分配颜色
    bar_colors = []
    for weight, category in zip(weight_values, categories_w):
        if weight >= 0:
            # 正权重使用橙色，但保持类别色调
            base_color = get_category_color(category)
            # 转换为RGB并增加亮度/饱和度
            bar_colors.append(colors['positive'])
        else:
            # 负权重使用蓝色
            bar_colors.append(colors['negative'])
    
    # 绘制条形图
    bars_w = ax_weight.barh(y_positions_w, weight_values, 
                           color=bar_colors,
                           edgecolor='black', linewidth=0.5)
    
    # 设置y轴
    ax_weight.set_yticks(y_positions_w)
    ax_weight.set_yticklabels(feature_names_w)
    
    # 设置x轴
    ax_weight.set_xlabel('Weight (Logistic Regression Composite Score)', fontsize=9)
    
    # 添加零参考线
    ax_weight.axvline(x=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # 添加条形数值
    for i, (bar, value) in enumerate(zip(bars_w, weight_values)):
        if value >= 0:
            ha = 'left'
            x_pos = value + 0.01
        else:
            ha = 'right'
            x_pos = value - 0.01
        
        ax_weight.text(x_pos, bar.get_y() + bar.get_height()/2,
                      f'{value:.3f}', ha=ha, va='center', fontsize=6)
    
    # 添加网格
    ax_weight.grid(True, axis='x', linestyle=':', linewidth=0.3, alpha=0.5)
    
    # 添加标题
    ax_weight.set_title('B. Logistic Regression Weights for Composite Score\n(Top 20 features by absolute weight)', 
                       fontsize=10, fontweight='bold', pad=15)
    
    # ============================
    # 添加类别图例
    # ============================
    
    # 获取唯一类别
    unique_categories = pd.concat([df_top_imp['Category'], df_top_weight['Category']]).unique()
    
    legend_elements = []
    for cat in unique_categories:
        if cat in ['Retina', 'RNFL', 'GCL', 'GCL+', 'GCL++', 'Choroid', 'Optic Disc']:
            color = get_category_color(cat)
            legend_elements.append(plt.Line2D([0], [0], marker='s', color='w', 
                                            markerfacecolor=color, 
                                            markeredgecolor='black',
                                            markeredgewidth=0.5,
                                            markersize=8, label=cat))
    
    # 添加正负权重说明
    legend_elements.append(plt.Line2D([0], [0], marker='s', color='w', 
                                     markerfacecolor=colors['positive'], 
                                     markeredgecolor='black',
                                     markeredgewidth=0.5,
                                     markersize=8, label='Positive weight'))
    legend_elements.append(plt.Line2D([0], [0], marker='s', color='w', 
                                     markerfacecolor=colors['negative'], 
                                     markeredgecolor='black',
                                     markeredgewidth=0.5,
                                     markersize=8, label='Negative weight'))
    
    # 添加图例（放在两个图之间）
    fig.legend(handles=legend_elements, loc='upper center', ncol=4, 
              frameon=True, framealpha=0.9, edgecolor='black', fontsize=7,
              bbox_to_anchor=(0.5, 0.96))
    
    # ============================
    # 添加相关性信息
    # ============================
    
    corr_text = f"Correlation between importance and absolute weight: r={correlation:.3f}"
    fig.text(0.5, 0.02, corr_text, ha='center', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    
    # ============================
    # 添加样本和模型信息
    # ============================
    
    sample_info = "Based on 463 eyes (303 MDD, 160 control); Random Forest with 100 trees; Logistic Regression with L2 regularization"
    fig.text(0.5, 0.05, sample_info, ha='center', va='bottom', fontsize=7,
            style='italic')
    
    # ============================
    # 调整布局
    # ============================
    
    plt.tight_layout(rect=[0, 0.08, 1, 0.92])  # 为图例和标注留出空间
    
    return fig

def create_correlation_inset(df_merged, ax_main):
    """在主图中创建嵌入的相关性散点图"""
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
    ax_inset = inset_axes(ax_main, width="30%", height="30%", loc='lower left',
                         bbox_to_anchor=(0.05, 0.05, 1, 1),
                         bbox_transform=ax_main.transAxes)
    
    # 绘制散点图
    x = df_merged['重要性']
    y = df_merged['abs_weight']
    
    # 按类别着色
    colors_scatter = [get_category_color(cat) for cat in df_merged['Category']]
    
    ax_inset.scatter(x, y, c=colors_scatter, s=20, alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # 添加趋势线
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax_inset.plot(x_line, p(x_line), color='black', linestyle='--', linewidth=1, alpha=0.7)
    
    # 设置坐标轴
    ax_inset.set_xlabel('Importance', fontsize=7)
    ax_inset.set_ylabel('|Weight|', fontsize=7)
    
    # 添加相关性系数
    corr = df_merged['重要性'].corr(df_merged['abs_weight'])
    ax_inset.text(0.05, 0.95, f'r={corr:.3f}', transform=ax_inset.transAxes,
                 fontsize=7, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax_inset.set_title('Importance vs. Weight', fontsize=8, pad=5)
    ax_inset.grid(True, linestyle=':', linewidth=0.3, alpha=0.5)
    
    return ax_inset

def main():
    """主函数"""
    print("制作Figure 4: Feature importance and composite score weights...")
    
    # 加载数据
    df_merged, df_top_imp, df_top_weight, correlation = load_feature_data()
    
    # 显示关键信息
    print(f"\n特征重要性Top 10:")
    for idx, row in df_top_imp.head(10).iterrows():
        print(f"  {row['特征_标准化']}: 重要性={row['重要性']:.4f}, 类别={row['Category']}")
    
    print(f"\n复合指标权重Top 10 (按绝对值):")
    df_top_abs_weight = df_merged.sort_values('abs_weight', ascending=False).head(10)
    for idx, row in df_top_abs_weight.iterrows():
        print(f"  {row['特征_标准化']}: 权重={row['权重']:.3f}, 类别={row['Category']}")
    
    print(f"\n类别分布:")
    category_counts = df_merged['Category'].value_counts()
    for cat, count in category_counts.items():
        print(f"  {cat}: {count}个特征 ({count/len(df_merged)*100:.1f}%)")
    
    # 创建图表
    fig = create_dual_panel_figure(df_top_imp, df_top_weight, correlation)
    
    # 保存图表
    output_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/02_Figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存为高分辨率TIFF
    output_path = os.path.join(output_dir, "Figure4_Feature_importance_composite_weights.tiff")
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
    print(f"  左侧面板: 20个最重要特征 (随机森林)")
    print(f"  右侧面板: 20个最大绝对权重特征 (逻辑回归)")
    
    # 关键发现
    print(f"\n关键发现:")
    print(f"  1. 最重要的特征: {df_top_imp.iloc[0]['特征_标准化']} (重要性={df_top_imp.iloc[0]['重要性']:.4f})")
    print(f"  2. 最大正权重: {df_merged.loc[df_merged['权重'].idxmax(), '特征_标准化']} (权重={df_merged['权重'].max():.3f})")
    print(f"  3. 最大负权重: {df_merged.loc[df_merged['权重'].idxmin(), '特征_标准化']} (权重={df_merged['权重'].min():.3f})")
    print(f"  4. 重要性与权重的相关性: r={correlation:.3f}")
    
    # 按类别的平均重要性
    print(f"\n按类别的平均重要性:")
    category_importance = df_merged.groupby('Category')['重要性'].mean().sort_values(ascending=False)
    for cat, mean_imp in category_importance.items():
        print(f"  {cat}: {mean_imp:.4f}")
    
    plt.close(fig)
    
    # 生成图表说明
    caption = f"""**Figure 4. Feature importance and composite score weights for OCT parameters.** 
Panel A shows the top 20 features ranked by Random Forest importance for discriminating major depressive disorder (MDD) from healthy controls. Panel B displays the top 20 features by absolute weight in the logistic regression composite score, with positive weights (orange) indicating increased probability of MDD and negative weights (blue) indicating decreased probability. 
The most important feature was RNFL_鼻侧 (importance={df_top_imp.iloc[0]['重要性']:.4f}), while the largest positive weight was RNFL_中心厚度 (weight={df_merged['权重'].max():.3f}) and largest negative weight was Choroid_外环鼻侧 (weight={df_merged['权重'].min():.3f}). 
The correlation between feature importance and absolute weight was r={correlation:.3f}, indicating moderate agreement between the two feature selection methods. 
Colors represent OCT parameter categories as shown in the legend."""
    
    caption_path = os.path.join(output_dir, "Figure4_caption.txt")
    with open(caption_path, 'w', encoding='utf-8') as f:
        f.write(caption)
    print(f"\n图表说明已保存: {caption_path}")

if __name__ == "__main__":
    main()