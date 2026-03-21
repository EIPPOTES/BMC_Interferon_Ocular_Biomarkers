import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("重新生成完整的英文表格")
print("=" * 80)

plt.rcParams['font.family'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

input_dir = '/mnt/c/Users/CUI/Desktop/最终修改'
output_dir = '/mnt/c/Users/CUI/Desktop/论文及图表'

def excel_to_png(df, png_path, title, col_widths=None):
    """将DataFrame转换为PNG图片"""
    try:
        n_rows = len(df) + 1
        n_cols = len(df.columns)
        
        # 根据内容调整列宽
        if col_widths is None:
            col_widths = [max(3, len(str(col)) * 0.15) for col in df.columns]
        
        fig_width = sum(col_widths) * 0.5 + 2
        fig_height = max(4, n_rows * 0.28 + 1)
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.axis('tight')
        ax.axis('off')
        
        cell_text = df.values.tolist()
        col_labels = df.columns.tolist()
        
        # 格式化数值
        formatted_cells = []
        for row in cell_text:
            formatted_row = []
            for cell in row:
                if isinstance(cell, float):
                    if abs(cell) < 0.001 and cell != 0:
                        formatted_row.append(f'{cell:.2e}')
                    elif abs(cell) < 0.01:
                        formatted_row.append(f'{cell:.4f}')
                    elif abs(cell) < 1:
                        formatted_row.append(f'{cell:.3f}')
                    else:
                        formatted_row.append(f'{cell:.2f}')
                else:
                    formatted_row.append(str(cell))
            formatted_cells.append(formatted_row)
        
        table = ax.table(
            cellText=formatted_cells,
            colLabels=col_labels,
            loc='center',
            cellLoc='center'
        )
        
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1, 1.6)
        
        # 设置表头样式
        for i in range(len(col_labels)):
            table[(0, i)].set_facecolor('#4472C4')
            table[(0, i)].set_text_props(color='white', fontweight='bold', fontsize=8)
        
        # 设置交替行颜色
        for i in range(1, len(formatted_cells) + 1):
            for j in range(len(col_labels)):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#E7E6E6')
                else:
                    table[(i, j)].set_facecolor('#FFFFFF')
        
        plt.title(title, fontsize=11, fontweight='bold', pad=12)
        plt.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        
        return True
    except Exception as e:
        print(f"  错误: {e}")
        import traceback
        traceback.print_exc()
        return False

# ==================== Table 1: Baseline Characteristics ====================
print("\n生成 Table 1...")
table1_data = {
    'Characteristic': [
        'Sample Size (n)',
        'Age (years)',
        'Sex, Female/Male (%)',
        'PHQ-9 Score',
        'GAD-7 Score',
        'Mean Macular Thickness (μm)',
        'Total Macular Volume (mm³)',
        'Foveal Thickness (μm)',
        'RNFL Total (μm)',
        'Disc Area (mm²)',
        'Cup Area (mm²)',
        'Rim Area (mm²)',
        'C/D Area Ratio'
    ],
    'MDD Group': [
        '325',
        '38.3 ± 20.2',
        '235/68 (72.3%/20.9%)',
        '8.63 ± 7.47 (n=260)',
        '6.78 ± 6.26 (n=260)',
        '271.45 ± 16.91',
        '7.67 ± 0.48',
        '223.91 ± 28.12',
        '105.89 ± 17.43',
        '2.06 ± 1.05',
        '0.67 ± 0.66',
        '1.39 ± 0.64',
        '0.30 ± 0.19'
    ],
    'Control Group': [
        '174',
        '28.0 ± 14.2',
        '102/60 (58.6%/34.5%)',
        '-',
        '-',
        '278.19 ± 14.89',
        '7.87 ± 0.42',
        '226.14 ± 23.57',
        '109.67 ± 13.28',
        '1.97 ± 0.83',
        '0.54 ± 0.54',
        '1.43 ± 0.58',
        '0.25 ± 0.18'
    ],
    'P-value': [
        '-', '0.0000', '-', '-', '-', '0.0000', '0.0000', '0.2924', '0.0465', '0.6269', '0.0103', '0.2016', '0.0082'
    ]
}
df1 = pd.DataFrame(table1_data)
excel_to_png(df1, f'{output_dir}/Table1_Baseline_Characteristics.png', 
             'Table 1. Baseline Characteristics', [4, 3.5, 3.5, 1.5])
print("  ✓ Table 1 完成")

# ==================== Table 2: Macular Layer Analysis ====================
print("\n生成 Table 2...")
df2 = pd.read_excel(f'{input_dir}/Table2_黄斑所有层_Statsmodels.xlsx')

# 翻译列名
df2.columns = ['Layer', 'Region', 'MDD Group', 'Control Group', 'P-value', "Cohen's d", 'P-value (FDR)', 'Reject H0', 'Significance']

# 翻译区域名称
region_map = {
    '黄斑中心凹': 'Fovea',
    '内环上方': 'Inner Superior',
    '内环颞侧': 'Inner Temporal',
    '内环鼻侧': 'Inner Nasal',
    '内环下方': 'Inner Inferior',
    '外环上方': 'Outer Superior',
    '外环颞侧': 'Outer Temporal',
    '外环鼻侧': 'Outer Nasal',
    '外环下方': 'Outer Inferior'
}
df2['Region'] = df2['Region'].map(region_map)

# 翻译层名称
layer_map = {
    'RNFL': 'RNFL',
    'Retina': 'Retina',
    'GCL+': 'GCL+',
    'GCL++': 'GCL++',
    'Choroid': 'Choroid'
}
df2['Layer'] = df2['Layer'].map(layer_map)

# 翻译显著性
df2['Significance'] = df2['Significance'].replace({'ns': 'ns'})
df2['Reject H0'] = df2['Reject H0'].map({True: 'Yes', False: 'No'})

excel_to_png(df2, f'{output_dir}/Table2_Macular_Layer_Analysis.png',
             'Table 2. Macular Layer Analysis', [1.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.2, 1.2])
print(f"  ✓ Table 2 完成 ({len(df2)} rows)")

# ==================== Table 3: Macular Comprehensive Metrics ====================
print("\n生成 Table 3...")
df3 = pd.read_excel(f'{input_dir}/Table3_黄斑综合指标_Statsmodels.xlsx')
df3.columns = ['Layer', 'Metric', 'MDD Group', 'Control Group', 'P-value', "Cohen's d", 'P-value (FDR)', 'Reject H0', 'Significance']

metric_map = {
    '平均厚度': 'Mean Thickness',
    '中心厚度': 'Central Thickness',
    '总体积': 'Total Volume'
}
df3['Metric'] = df3['Metric'].map(metric_map)
df3['Layer'] = df3['Layer'].map(layer_map)
df3['Reject H0'] = df3['Reject H0'].map({True: 'Yes', False: 'No'})

excel_to_png(df3, f'{output_dir}/Table3_Macular_Comprehensive_Metrics.png',
             'Table 3. Macular Comprehensive Metrics', [1.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.2, 1.2])
print(f"  ✓ Table 3 完成 ({len(df3)} rows)")

# ==================== Table 4: Optic Disc Parameters ====================
print("\n生成 Table 4...")
df4 = pd.read_excel(f'{input_dir}/Table4_视盘所有指标组间比较.xlsx')
df4.columns = ['Category', 'Parameter', 'MDD Group', 'Control Group', 'P-value', "Cohen's d", 'P-value (FDR)', 'Significance']

category_map = {
    'RNFL': 'RNFL',
    '视盘结构': 'Optic Disc'
}
df4['Category'] = df4['Category'].map(category_map)

excel_to_png(df4, f'{output_dir}/Table4_Optic_Disc_Parameters.png',
             'Table 4. Optic Disc Parameters', [1.5, 3, 2.5, 2.5, 1.5, 1.5, 1.5, 1.2])
print(f"  ✓ Table 4 完成 ({len(df4)} rows)")

# ==================== Table 5: Correlation Analysis ====================
print("\n生成 Table 5...")
df5 = pd.read_excel(f'{input_dir}/Table5_全面相关分析.xlsx')
df5.columns = ['Layer/Category', 'Region/Parameter', 'n', 'Spearman r', 'P-value', 'Parameter2', 'Category2']

# 翻译区域和参数
df5['Region/Parameter'] = df5['Region/Parameter'].replace(region_map)
df5['Layer/Category'] = df5['Layer/Category'].replace(layer_map).replace(category_map)

# 选择关键列
df5_display = df5[['Layer/Category', 'Region/Parameter', 'n', 'Spearman r', 'P-value']].head(25)

excel_to_png(df5_display, f'{output_dir}/Table5_Correlation_Analysis.png',
             'Table 5. Correlation Analysis (PHQ-9 vs OCT Parameters, Top 25)', [2, 3, 1, 1.5, 1.5])
print(f"  ✓ Table 5 完成 (showing top 25 of {len(df5)} rows)")

# ==================== Table 6: ROC Analysis ====================
print("\n生成 Table 6...")
df6 = pd.read_excel(f'{input_dir}/Table6_全面ROC分析.xlsx')
df6.columns = ['Parameter', 'AUC', 'Sensitivity', 'Specificity', 'Optimal Cutoff']

# 翻译参数名中的区域
for cn, en in region_map.items():
    df6['Parameter'] = df6['Parameter'].str.replace(cn, en, regex=False)

excel_to_png(df6.head(20), f'{output_dir}/Table6_ROC_Analysis.png',  # 只显示前20行
             'Table 6. ROC Analysis (Top 20 Parameters)', [4, 1.5, 1.5, 1.5, 2])
print(f"  ✓ Table 6 完成 (showing top 20 of {len(df6)} rows)")

# ==================== Table 7: Subgroup Analysis ====================
print("\n生成 Table 7...")
df7 = pd.read_excel(f'{input_dir}/Table7_亚组分析_严重程度.xlsx')
df7.columns = ['Parameter', 'No Depression (n=103)', 'Mild (n=54)', 'Moderate (n=40)', 'Severe (n=63)', 
               'Kruskal P', 'Trend r', 'Trend P', 'No Dep2 (n=102)', 'Mild2 (n=52)', 'Severe2 (n=60)']

# 翻译参数名
for cn, en in region_map.items():
    df7['Parameter'] = df7['Parameter'].str.replace(cn, en, regex=False)

# 只显示关键列
df7_display = df7[['Parameter', 'No Depression (n=103)', 'Mild (n=54)', 'Moderate (n=40)', 'Severe (n=63)', 'Kruskal P']]

excel_to_png(df7_display, f'{output_dir}/Table7_Subgroup_Analysis.png',
             'Table 7. Subgroup Analysis by Depression Severity', [4, 2.2, 2.2, 2.2, 2.2, 1.5])
print(f"  ✓ Table 7 完成 ({len(df7)} rows)")

# ==================== Table 8: Multivariate Regression ====================
print("\n生成 Table 8...")
df8 = pd.read_excel(f'{input_dir}/Table8_多因素回归_Statsmodels.xlsx')
df8.columns = ['Outcome', 'R²', 'Adj R²', 'F-stat', 'F P-value', 'Depression β', 'Depression SE', 'Depression t', 'Depression P', 'Age β', 'Age P', 'Sex β', 'Sex P']

outcome_map = {
    '黄斑平均厚度': 'Mean Macular Thickness',
    '黄斑总体积': 'Total Macular Volume'
}
df8['Outcome'] = df8['Outcome'].map(outcome_map)

excel_to_png(df8, f'{output_dir}/Table8_Multivariate_Regression.png',
             'Table 8. Multivariate Regression Analysis', [4, 1.2, 1.5, 1.5, 1.5, 2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5])
print(f"  ✓ Table 8 完成 ({len(df8)} rows)")

# ==================== Table 9: Inter-eye Consistency ====================
print("\n生成 Table 9...")
df9 = pd.read_excel(f'{input_dir}/Table9_双眼一致性分析.xlsx')
df9.columns = ['Parameter', 'Right Eye', 'Left Eye', 'Pearson r', 'Mean Absolute Difference', 'ICC']

# 翻译参数名
for cn, en in region_map.items():
    df9['Parameter'] = df9['Parameter'].str.replace(cn, en, regex=False)

excel_to_png(df9, f'{output_dir}/Table9_Intereye_Consistency.png',
             'Table 9. Inter-eye Consistency', [4, 2.5, 2.5, 1.5, 3, 1.5])
print(f"  ✓ Table 9 完成 ({len(df9)} rows)")

# ==================== Table 10: ML Model Comparison ====================
print("\n生成 Table 10...")
df10 = pd.read_excel(f'{input_dir}/Table10_机器学习模型比较.xlsx')
df10.columns = ['Model', 'AUC Mean', 'AUC Std', 'Accuracy Mean', 'Accuracy Std']

excel_to_png(df10, f'{output_dir}/Table10_ML_Model_Comparison.png',
             'Table 10. Machine Learning Model Comparison (5-Fold CV)', [4, 2, 2, 2.5, 2.5])
print(f"  ✓ Table 10 完成 ({len(df10)} rows)")

# ==================== Table 11: Feature Importance ====================
print("\n生成 Table 11...")
df11 = pd.read_excel(f'{input_dir}/Table11_特征重要性.xlsx')
df11.columns = ['Feature', 'Importance', 'Rank']

# 翻译特征名中的区域
for cn, en in region_map.items():
    df11['Feature'] = df11['Feature'].str.replace(cn, en, regex=False)

excel_to_png(df11.head(15), f'{output_dir}/Table11_Feature_Importance.png',  # 只显示前15行
             'Table 11. Feature Importance (Random Forest, Top 15)', [5, 2, 1.5])
print(f"  ✓ Table 11 完成 (showing top 15 of {len(df11)} rows)")

print(f"\n{'=' * 80}")
print("所有英文表格生成完成!")
print(f"{'=' * 80}")
