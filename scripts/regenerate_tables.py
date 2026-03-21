import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
def main():
    """主函数，包装原有执行代码"""
    warnings.filterwarnings('ignore')

    print("=" * 80)
    print("重新生成英文表格 - 修复布局问题")
    print("=" * 80)

    plt.rcParams['font.family'] = ['DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    output_dir = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    input_dir = str(/mnt/c/Users/CUI/Desktop/论文及图表)

    def create_table_image(df, output_path, title, fig_width=None, fig_height=None):
        """创建表格图片"""
        n_rows = len(df)
        n_cols = len(df.columns)

        # 自动计算尺寸
        if fig_width is None:
            fig_width = max(10, n_cols * 1.8)
        if fig_height is None:
            fig_height = max(4, n_rows * 0.35 + 1.5)

        fig = plt.figure(figsize=(fig_width, fig_height))
        ax = fig.add_subplot(111)
        ax.axis('off')

        # 准备数据
        cell_text = []
        for _, row in df.iterrows():
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
            cell_text.append(formatted_row)

        # 创建表格
        table = ax.table(
            cellText=cell_text,
            colLabels=df.columns.tolist(),
            loc='upper center',
            cellLoc='center',
            bbox=[0, 0, 1, 1]
        )

        # 设置字体大小
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.8)

        # 表头样式
        for i in range(len(df.columns)):
            table[(0, i)].set_facecolor('#4472C4')
            table[(0, i)].set_text_props(color='white', fontweight='bold', fontsize=9)

        # 交替行颜色
        for i in range(1, n_rows + 1):
            for j in range(len(df.columns)):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#F2F2F2')
                else:
                    table[(i, j)].set_facecolor('#FFFFFF')

        # 添加标题（在图上方，不覆盖表格）
        fig.suptitle(title, fontsize=13, fontweight='bold', y=0.98)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.3)
        plt.close()
        return True

    # ==================== Table 1 ====================
    print("\nTable 1...")
    table1_data = {
        'Characteristic': [
            'Participants (n)', 'Eyes (n)', 'Age (years)', 'Sex, Female/Male (%)',
            'PHQ-9 Score (n=132 patients)', 'GAD-7 Score (n=132 patients)', 'Mean Macular Thickness (μm)',
            'Total Macular Volume (mm³)', 'Foveal Thickness (μm)',
            'RNFL Total (μm)', 'Disc Area (mm²)', 'Cup Area (mm²)',
            'Rim Area (mm²)', 'C/D Area Ratio'
        ],
        'MDD Group': [
            '164', '325', '38.3 ± 20.2', '235/68 (72.3%/20.9%)',
            '8.63 ± 7.47', '6.78 ± 6.26',
            '271.45 ± 16.91', '7.67 ± 0.48', '223.91 ± 28.12',
            '105.89 ± 17.43', '2.06 ± 1.05', '0.67 ± 0.66',
            '1.39 ± 0.64', '0.30 ± 0.19'
        ],
        'Control Group': [
            '87', '174', '28.0 ± 14.2', '102/60 (58.6%/34.5%)',
            '-', '-', '278.19 ± 14.89', '7.87 ± 0.42',
            '226.14 ± 23.57', '109.67 ± 13.28', '1.97 ± 0.83',
            '0.54 ± 0.54', '1.43 ± 0.58', '0.25 ± 0.18'
        ],
        'P-value': ['-', '-', '0.0000', '-', '-', '-', '0.0000', '0.0000', '0.2924', '0.0465', '0.6269', '0.0103', '0.2016', '0.0082']
    }
    df1 = pd.DataFrame(table1_data)
    create_table_image(df1, f'{output_dir}/Table1_Baseline_Characteristics.png', 'Table 1. Baseline Characteristics', 13, 8)
    print("  ✓ Done")

    # ==================== Table 2 ====================
    print("\nTable 2...")
    df2 = pd.read_excel(f'{input_dir}/Table2_黄斑所有层_Statsmodels.xlsx')
    df2.columns = ['Layer', 'Region', 'MDD Group', 'Control Group', 'P-value', "Cohen's d", 'P-value (FDR)', 'Reject H0', 'Significance']

    region_map = {'黄斑中心凹': 'Fovea', '内环上方': 'Inner Superior', '内环颞侧': 'Inner Temporal', '内环鼻侧': 'Inner Nasal', '内环下方': 'Inner Inferior', '外环上方': 'Outer Superior', '外环颞侧': 'Outer Temporal', '外环鼻侧': 'Outer Nasal', '外环下方': 'Outer Inferior'}
    df2['Region'] = df2['Region'].map(region_map)
    df2['Reject H0'] = df2['Reject H0'].map({True: 'Yes', False: 'No'})

    create_table_image(df2, f'{output_dir}/Table2_Macular_Layer_Analysis.png', 'Table 2. Macular Layer Analysis', 14, 18)
    print(f"  ✓ Done ({len(df2)} rows)")

    # ==================== Table 3 ====================
    print("\nTable 3...")
    df3 = pd.read_excel(f'{input_dir}/Table3_黄斑综合指标_Statsmodels.xlsx')
    df3.columns = ['Layer', 'Metric', 'MDD Group', 'Control Group', 'P-value', "Cohen's d", 'P-value (FDR)', 'Reject H0', 'Significance']
    metric_map = {'平均厚度': 'Mean Thickness', '中心厚度': 'Central Thickness', '总体积': 'Total Volume'}
    df3['Metric'] = df3['Metric'].map(metric_map)
    df3['Reject H0'] = df3['Reject H0'].map({True: 'Yes', False: 'No'})

    create_table_image(df3, f'{output_dir}/Table3_Macular_Comprehensive_Metrics.png', 'Table 3. Macular Comprehensive Metrics', 14, 8)
    print(f"  ✓ Done ({len(df3)} rows)")

    # ==================== Table 4 ====================
    print("\nTable 4...")
    df4 = pd.read_excel(f'{input_dir}/Table4_视盘所有指标组间比较.xlsx')
    df4.columns = ['Category', 'Parameter', 'MDD Group', 'Control Group', 'P-value', "Cohen's d", 'P-value (FDR)', 'Significance']
    df4['Category'] = df4['Category'].map({'RNFL': 'RNFL', '视盘结构': 'Optic Disc'})

    # 翻译RNFL参数名
    param_map = {
        'RNFL_Total': 'RNFL Total',
        'RNFL_上方': 'RNFL Superior',
        'RNFL_颞侧': 'RNFL Temporal',
        'RNFL_鼻侧': 'RNFL Nasal',
        'RNFL_下方': 'RNFL Inferior'
    }
    df4['Parameter'] = df4['Parameter'].replace(param_map)

    create_table_image(df4, f'{output_dir}/Table4_Optic_Disc_Parameters.png', 'Table 4. Optic Disc Parameters', 13, 7)
    print(f"  ✓ Done ({len(df4)} rows)")

    # ==================== Table 5 ====================
    print("\nTable 5...")
    df5 = pd.read_excel(f'{input_dir}/Table5_全面相关分析.xlsx')
    df5 = df5[['层', '区域', 'n', 'Spearman_r', 'P值']].head(20)
    df5.columns = ['Layer', 'Region', 'n', 'Spearman r', 'P-value']
    df5['Region'] = df5['Region'].map(region_map)

    create_table_image(df5, f'{output_dir}/Table5_Correlation_Analysis.png', 'Table 5. Correlation Analysis (PHQ-9 vs OCT, Top 20)', 12, 9)
    print(f"  ✓ Done ({len(df5)} rows)")

    # ==================== Table 6 ====================
    print("\nTable 6...")
    df6 = pd.read_excel(f'{input_dir}/Table6_全面ROC分析.xlsx').head(15)
    df6.columns = ['Parameter', 'AUC', 'Sensitivity', 'Specificity', 'Optimal Cutoff']

    # 翻译参数名中的中文区域
    region_map = {'黄斑中心凹': 'Fovea', '内环上方': 'Inner Superior', '内环颞侧': 'Inner Temporal', '内环鼻侧': 'Inner Nasal', '内环下方': 'Inner Inferior', '外环上方': 'Outer Superior', '外环颞侧': 'Outer Temporal', '外环鼻侧': 'Outer Nasal', '外环下方': 'Outer Inferior', '平均厚度': 'Mean Thickness', '中心厚度': 'Central Thickness', '总体积': 'Total Volume', '上方': 'Superior', '下方': 'Inferior', '颞侧': 'Temporal', '鼻侧': 'Nasal'}
    for cn, en in region_map.items():
        df6['Parameter'] = df6['Parameter'].str.replace(cn, en, regex=False)

    create_table_image(df6, f'{output_dir}/Table6_ROC_Analysis.png', 'Table 6. ROC Analysis (Top 15 Parameters)', 12, 8)
    print(f"  ✓ Done ({len(df6)} rows)")

    # ==================== Table 7 ====================
    print("\nTable 7...")
    df7 = pd.read_excel(f'{input_dir}/Table7_亚组分析_严重程度.xlsx')
    df7_display = df7[['指标', '无抑郁(n=103)', '轻度(n=54)', '中度(n=40)', '重度(n=63)', 'Kruskal_P']]
    df7_display.columns = ['Parameter', 'No Depression', 'Mild', 'Moderate', 'Severe', 'P-value']

    # 翻译参数名
    param_map = {'黄斑平均厚度': 'Mean Macular Thickness', '黄斑总体积': 'Total Macular Volume', '外环颞侧': 'Outer Temporal', 'RNFL Total': 'RNFL Total'}
    df7_display['Parameter'] = df7_display['Parameter'].replace(param_map)

    create_table_image(df7_display, f'{output_dir}/Table7_Subgroup_Analysis.png', 'Table 7. Subgroup Analysis by Depression Severity', 13, 5)
    print(f"  ✓ Done ({len(df7_display)} rows)")

    # ==================== Table 8 ====================
    print("\nTable 8...")
    df8 = pd.read_excel(f'{input_dir}/Table8_多因素回归_Statsmodels.xlsx')
    df8.columns = ['Outcome', 'R²', 'Adj R²', 'F-stat', 'F P-value', 'Depression β', 'Depression SE', 'Depression t', 'Depression P', 'Age β', 'Age P', 'Sex β', 'Sex P']
    df8['Outcome'] = df8['Outcome'].map({'黄斑平均厚度': 'Mean Macular Thickness', '黄斑总体积': 'Total Macular Volume'})

    create_table_image(df8, f'{output_dir}/Table8_Multivariate_Regression.png', 'Table 8. Multivariate Regression Analysis', 16, 4)
    print(f"  ✓ Done ({len(df8)} rows)")

    # ==================== Table 9 ====================
    print("\nTable 9...")
    df9 = pd.read_excel(f'{input_dir}/Table9_双眼一致性分析.xlsx')
    df9.columns = ['Parameter', 'Sample Size', 'Pearson r', 'P-value', 'Mean Absolute Difference (μm)', 'ICC']

    # 翻译参数名
    param_map = {'黄斑平均厚度': 'Mean Macular Thickness'}
    df9['Parameter'] = df9['Parameter'].replace(param_map)

    create_table_image(df9, f'{output_dir}/Table9_Intereye_Consistency.png', 'Table 9. Inter-eye Consistency', 14, 4)
    print(f"  ✓ Done ({len(df9)} rows)")

    # ==================== Table 10 ====================
    print("\nTable 10...")
    df10 = pd.read_excel(f'{input_dir}/Table10_机器学习模型比较.xlsx')
    df10.columns = ['Model', 'AUC Mean', 'AUC Std', 'Accuracy Mean', 'Accuracy Std']

    create_table_image(df10, f'{output_dir}/Table10_ML_Model_Comparison.png', 'Table 10. Machine Learning Model Comparison (5-Fold CV)', 12, 5)
    print(f"  ✓ Done ({len(df10)} rows)")

    # ==================== Table 11 ====================
    print("\nTable 11...")
    df11 = pd.read_excel(f'{input_dir}/Table11_特征重要性.xlsx').head(15)
    df11.columns = ['Feature', 'Coefficient', 'Abs Coefficient']

    # 翻译特征名中的中文区域
    region_map = {'黄斑中心凹': 'Fovea', '内环上方': 'Inner Superior', '内环颞侧': 'Inner Temporal', '内环鼻侧': 'Inner Nasal', '内环下方': 'Inner Inferior', '外环上方': 'Outer Superior', '外环颞侧': 'Outer Temporal', '外环鼻侧': 'Outer Nasal', '外环下方': 'Outer Inferior', '平均厚度': 'Mean Thickness', '中心厚度': 'Central Thickness', '总体积': 'Total Volume', '上方': 'Superior', '下方': 'Inferior', '颞侧': 'Temporal', '鼻侧': 'Nasal'}
    for cn, en in region_map.items():
        df11['Feature'] = df11['Feature'].str.replace(cn, en, regex=False)

    create_table_image(df11, f'{output_dir}/Table11_Feature_Importance.png', 'Table 11. Feature Importance (Random Forest, Top 15)', 11, 7)
    print(f"  ✓ Done ({len(df11)} rows)")

    print(f"\n{'=' * 80}")
    print("All tables regenerated successfully!")
    print(f"{'=' * 80}")



if __name__ == "__main__":
    main()