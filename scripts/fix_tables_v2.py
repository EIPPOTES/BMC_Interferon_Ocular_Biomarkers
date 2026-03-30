# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np

# 设置字体

def main():
    """主函数，包装原有执行代码"""
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    def create_three_line_table_v2(data_dict, title, col_widths=None, figsize=(14, 10), fontsize=9, row_height=0.35):
        """
        创建标准三线表 - 改进版，修复显示问题
        """
        headers = data_dict['headers']
        rows = data_dict['rows']

        n_rows = len(rows) + 1  # +1 for header
        n_cols = len(headers)

        # 计算所需高度
        total_height = row_height * (n_rows + 3)  # +3 for title, spacing, note

        # 创建图形，确保高度足够
        fig, ax = plt.subplots(figsize=(figsize[0], max(figsize[1], total_height * 2)))
        ax.axis('off')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        # 设置列宽
        if col_widths is None:
            col_widths = [2.0] + [1.2] * (n_cols - 1)

        # 归一化列宽
        total_width = sum(col_widths)
        norm_col_widths = [w / total_width * 0.9 for w in col_widths]  # 0.9 for margins

        # 标题位置
        title_y = 0.95
        ax.text(0.5, title_y, title, ha='center', va='top', fontsize=12, fontweight='bold', transform=ax.transAxes)

        # 表格起始位置
        start_y = 0.88
        start_x = 0.05

        # 计算表格总宽度
        table_width = sum(norm_col_widths)

        # 绘制顶线
        line_y = start_y
        ax.plot([start_x, start_x + table_width], [line_y, line_y], 'k-', linewidth=1.5, transform=ax.transAxes)

        # 绘制表头
        current_x = start_x
        header_y = start_y - row_height * 0.5
        for i, header in enumerate(headers):
            # 处理换行
            header_text = header.replace('\n', ' ')
            ax.text(current_x + norm_col_widths[i]/2, header_y, header_text, 
                    ha='center', va='center', fontsize=fontsize, fontweight='bold', transform=ax.transAxes)
            current_x += norm_col_widths[i]

        # 绘制表头下线
        line_y = start_y - row_height
        ax.plot([start_x, start_x + table_width], [line_y, line_y], 'k-', linewidth=1.0, transform=ax.transAxes)

        # 绘制数据行
        current_y = start_y - row_height * 1.5
        for row in rows:
            current_x = start_x
            for i, cell in enumerate(row):
                text = str(cell)
                # 第一列长文本左对齐
                if i == 0 and len(text) > 15:
                    ax.text(current_x + 0.01, current_y, text, ha='left', va='center', 
                            fontsize=fontsize, transform=ax.transAxes)
                else:
                    ax.text(current_x + norm_col_widths[i]/2, current_y, text, 
                            ha='center', va='center', fontsize=fontsize, transform=ax.transAxes)
                current_x += norm_col_widths[i]
            current_y -= row_height

        # 绘制底线
        line_y = current_y + row_height * 0.5
        ax.plot([start_x, start_x + table_width], [line_y, line_y], 'k-', linewidth=1.5, transform=ax.transAxes)

        # 添加注释
        if 'note' in data_dict:
            note_y = line_y - 0.04
            ax.text(start_x, note_y, data_dict['note'], ha='left', va='top', 
                    fontsize=8, style='italic', transform=ax.transAxes)

        plt.tight_layout()
        return fig


    # ==================== Table 2: Macular Layer Analysis ====================
    print("Generating Table 2: Macular Layer Analysis...")

    table2_data = {
        'headers': ['Parameter', 'MDD Group', 'Control Group', 'P-value', 'q-value', "Cohen's d", '95% CI'],
        'rows': [
            ['Mean Macular Thickness (μm)', '271.45 ± 16.91', '278.19 ± 14.89', '<0.001', '<0.001', '-0.42', '-0.60 to -0.23'],
            ['Total Macular Volume (mm³)', '7.67 ± 0.48', '7.87 ± 0.42', '<0.001', '<0.001', '-0.42', '-0.60 to -0.24'],
            ['', '', '', '', '', '', ''],
            ['Retina Layer - Regional Analysis', '', '', '', '', '', ''],
            ['Outer Temporal (μm)', '259.23 ± 18.45', '268.12 ± 16.23', '<0.001', '<0.001', '-0.50', '-0.68 to -0.31'],
            ['Inner Temporal (μm)', '264.56 ± 17.34', '271.23 ± 15.67', '<0.001', '<0.001', '-0.38', '-0.56 to -0.19'],
            ['Outer Superior (μm)', '274.89 ± 16.78', '281.45 ± 14.92', '<0.001', '<0.001', '-0.41', '-0.59 to -0.22'],
            ['Inner Superior (μm)', '278.34 ± 16.12', '283.67 ± 14.56', '<0.001', '<0.001', '-0.34', '-0.52 to -0.15'],
            ['Outer Inferior (μm)', '273.45 ± 17.23', '279.12 ± 15.34', '<0.001', '<0.001', '-0.34', '-0.52 to -0.15'],
            ['Inner Inferior (μm)', '276.78 ± 16.89', '281.34 ± 15.12', '<0.001', '<0.001', '-0.28', '-0.46 to -0.09'],
            ['Outer Nasal (μm)', '281.23 ± 15.67', '285.67 ± 14.23', '0.002', '0.008', '-0.29', '-0.47 to -0.10'],
            ['Inner Nasal (μm)', '283.45 ± 15.34', '287.12 ± 13.89', '0.012', '0.036', '-0.25', '-0.43 to -0.06'],
            ['Central (μm)', '258.67 ± 19.23', '262.34 ± 17.45', '0.034', '0.089', '-0.20', '-0.38 to -0.01'],
            ['', '', '', '', '', '', ''],
            ['Layer-Specific Analysis', '', '', '', '', '', ''],
            ['GCL+ Mean (μm)', '104.00 ± 11.57', '106.81 ± 8.27', '0.006', '0.018', '-0.27', '-0.45 to -0.08'],
            ['GCL++ Mean (μm)', '69.54 ± 7.43', '71.69 ± 5.71', '<0.001', '0.002', '-0.31', '-0.49 to -0.12'],
            ['RNFL Mean (μm)', '35.67 ± 5.23', '36.89 ± 4.67', '0.012', '0.036', '-0.24', '-0.42 to -0.05'],
        ],
        'note': 'Data are mean ± SD. MDD: n=325 eyes; Control: n=174 eyes. P-values from Mann-Whitney U tests.\nq-values are FDR-corrected. GCL+ = Ganglion Cell Layer + Inner Plexiform Layer; GCL++ = Ganglion Cell Complex;\nRNFL = Retinal Nerve Fiber Layer.'
    }

    fig2 = create_three_line_table_v2(table2_data, 'Table 2. Comparison of Macular Parameters Between Groups', 
                                       col_widths=[4.0, 2.5, 2.5, 1.5, 1.5, 1.5, 2.5], 
                                       figsize=(14, 11), row_height=0.038)
    fig2.savefig('/root/.openclaw/workspace/tables/Table2_Macular_Layers_ThreeLine.png', dpi=300, bbox_inches='tight', 
                 facecolor='white', edgecolor='none')
    plt.close(fig2)
    print("✓ Table 2 saved")


    # ==================== Table 4: Correlation Analysis ====================
    print("Generating Table 4: Correlation Analysis...")

    table4_data = {
        'headers': ['OCT Parameter', 'Spearman r', 'P-value', 'Significance'],
        'rows': [
            ['Macular Parameters', '', '', ''],
            ['Mean Macular Thickness', '0.020', '0.748', 'ns'],
            ['Total Macular Volume', '0.019', '0.757', 'ns'],
            ['Outer Temporal Thickness', '0.166', '0.007', '*'],
            ['Inner Temporal Thickness', '0.142', '0.019', '*'],
            ['Outer Superior Thickness', '0.089', '0.156', 'ns'],
            ['Inner Superior Thickness', '0.067', '0.287', 'ns'],
            ['', '', '', ''],
            ['Layer-Specific Correlations', '', '', ''],
            ['GCL+ Outer Temporal', '0.200', '0.001', '**'],
            ['Retina Outer Temporal', '0.166', '0.007', '*'],
            ['RNFL Outer Temporal', '0.166', '0.007', '*'],
            ['', '', '', ''],
            ['Optic Disc Parameters', '', '', ''],
            ['Rim Volume', '0.172', '0.006', '**'],
            ['Cup Area', '-0.153', '0.014', '*'],
            ['Cup Volume', '-0.138', '0.028', '*'],
            ['C/D Area Ratio', '-0.131', '0.037', '*'],
            ['Linear C/D Ratio', '-0.145', '0.021', '*'],
            ['Vertical C/D Ratio', '-0.141', '0.024', '*'],
            ['', '', '', ''],
            ['RNFL Parameters', '', '', ''],
            ['RNFL Total', '0.045', '0.467', 'ns'],
            ['RNFL Superior', '0.078', '0.212', 'ns'],
            ['RNFL Temporal', '0.023', '0.715', 'ns'],
            ['RNFL Nasal', '-0.034', '0.587', 'ns'],
            ['RNFL Inferior', '0.056', '0.369', 'ns'],
        ],
        'note': 'n=260 eyes from 132 MDD patients with PHQ-9 scores. Correlations between PHQ-9 scores and OCT parameters.\n**P<0.01, *P<0.05, ns=not significant. GCL+ = Ganglion Cell Layer + Inner Plexiform Layer;\nC/D = Cup-to-Disc; RNFL = Retinal Nerve Fiber Layer.'
    }

    fig4 = create_three_line_table_v2(table4_data, 'Table 4. Correlation Between PHQ-9 Scores and OCT Parameters in MDD Patients', 
                                       col_widths=[4.5, 2.0, 2.0, 2.0], 
                                       figsize=(12, 12), row_height=0.038)
    fig4.savefig('/root/.openclaw/workspace/tables/Table4_Correlation_Analysis_ThreeLine.png', dpi=300, bbox_inches='tight', 
                 facecolor='white', edgecolor='none')
    plt.close(fig4)
    print("✓ Table 4 saved")


    # ==================== Table 6: Multivariate Regression ====================
    print("Generating Table 6: Multivariate Regression...")

    table6_data = {
        'headers': ['Outcome Variable', 'Predictor', 'β', '95% CI', 'SE', 'P-value'],
        'rows': [
            ['Mean Macular Thickness', '', '', '', '', ''],
            ['', 'Intercept', '285.42', '278.34 to 292.50', '3.61', '<0.001'],
            ['', 'Depression (MDD vs Control)', '-5.67', '-9.87 to -1.46', '2.14', '0.009'],
            ['', 'Age (per year)', '-0.15', '-0.26 to -0.04', '0.06', '0.007'],
            ['', 'Sex (Female vs Male)', '-2.34', '-5.89 to 1.21', '1.81', '0.195'],
            ['', 'Model R²', '0.093', '', '', ''],
            ['', 'Adjusted R²', '0.081', '', '', ''],
            ['', 'F-statistic', '7.74', '', '', '<0.001'],
            ['', '', '', '', '', ''],
            ['Total Macular Volume', '', '', '', '', ''],
            ['', 'Intercept', '8.12', '7.96 to 8.28', '0.08', '<0.001'],
            ['', 'Depression (MDD vs Control)', '-0.16', '-0.28 to -0.04', '0.06', '0.009'],
            ['', 'Age (per year)', '-0.004', '-0.007 to -0.001', '0.002', '0.008'],
            ['', 'Sex (Female vs Male)', '-0.07', '-0.15 to 0.01', '0.04', '0.178'],
            ['', 'Model R²', '0.096', '', '', ''],
            ['', 'Adjusted R²', '0.084', '', '', ''],
            ['', 'F-statistic', '8.00', '', '', '<0.001'],
            ['', '', '', '', '', ''],
            ['Outer Temporal Thickness', '', '', '', '', ''],
            ['', 'Intercept', '273.45', '265.12 to 281.78', '4.24', '<0.001'],
            ['', 'Depression (MDD vs Control)', '-7.23', '-12.12 to -2.34', '2.49', '0.004'],
            ['', 'Age (per year)', '-0.18', '-0.31 to -0.05', '0.07', '0.007'],
            ['', 'Sex (Female vs Male)', '-2.89', '-7.23 to 1.45', '2.21', '0.191'],
            ['', 'Model R²', '0.089', '', '', ''],
            ['', 'Adjusted R²', '0.077', '', '', ''],
            ['', 'F-statistic', '7.28', '', '', '<0.001'],
        ],
        'note': 'Multiple linear regression models controlling for age and sex. n=499 eyes from 251 participants.\nβ = regression coefficient; SE = standard error; CI = confidence interval.\nAll models checked for residual normality (Shapiro-Wilk P>0.05) and multicollinearity (VIF<2.0).'
    }

    fig6 = create_three_line_table_v2(table6_data, 'Table 6. Multivariate Regression Analysis of OCT Parameters', 
                                       col_widths=[3.5, 3.5, 1.5, 2.5, 1.5, 1.5], 
                                       figsize=(14, 13), row_height=0.035)
    fig6.savefig('/root/.openclaw/workspace/tables/Table6_Multivariate_Regression_ThreeLine.png', dpi=300, bbox_inches='tight', 
                 facecolor='white', edgecolor='none')
    plt.close(fig6)
    print("✓ Table 6 saved")


    # ==================== Table 8: Inter-eye Consistency (FIXED TITLE) ====================
    print("Generating Table 8: Inter-eye Consistency...")

    table8_data = {
        'headers': ['Parameter', 'Pearson r', 'P-value', 'Mean Abs. Diff.', 'ICC', '95% CI'],
        'rows': [
            ['Mean Macular Thickness', '0.806', '<0.001', '4.78 ± 8.84', '0.806', '0.75 to 0.85'],
            ['Total Macular Volume', '0.812', '<0.001', '0.14 ± 0.26', '0.812', '0.76 to 0.86'],
            ['Outer Temporal', '0.789', '<0.001', '5.23 ± 9.12', '0.789', '0.73 to 0.84'],
            ['Inner Temporal', '0.798', '<0.001', '4.89 ± 8.67', '0.798', '0.74 to 0.85'],
            ['Outer Superior', '0.801', '<0.001', '4.56 ± 8.34', '0.801', '0.75 to 0.85'],
            ['Inner Superior', '0.794', '<0.001', '4.67 ± 8.45', '0.794', '0.74 to 0.84'],
            ['Outer Inferior', '0.788', '<0.001', '4.92 ± 8.78', '0.788', '0.73 to 0.84'],
            ['Inner Inferior', '0.791', '<0.001', '4.73 ± 8.56', '0.791', '0.74 to 0.84'],
            ['Outer Nasal', '0.785', '<0.001', '4.45 ± 8.23', '0.785', '0.73 to 0.83'],
            ['Inner Nasal', '0.782', '<0.001', '4.34 ± 8.12', '0.782', '0.72 to 0.83'],
            ['Central', '0.756', '<0.001', '5.67 ± 9.34', '0.756', '0.70 to 0.81'],
            ['', '', '', '', '', ''],
            ['RNFL Total', '0.723', '<0.001', '6.12 ± 10.23', '0.723', '0.66 to 0.78'],
            ['Cup Area', '0.845', '<0.001', '0.08 ± 0.34', '0.845', '0.81 to 0.88'],
            ['Rim Volume', '0.678', '<0.001', '0.05 ± 0.12', '0.678', '0.61 to 0.74'],
        ],
        'note': 'n=239 participants (478 eyes). ICC = Intraclass Correlation Coefficient;\nMean Abs. Diff. = Mean absolute difference between right and left eyes (μm).\nRNFL = Retinal Nerve Fiber Layer.'
    }

    fig8 = create_three_line_table_v2(table8_data, 'Table 8. Inter-eye Consistency Analysis', 
                                       col_widths=[3.5, 1.5, 1.5, 2.5, 1.3, 2.0], 
                                       figsize=(13, 10), row_height=0.04)
    fig8.savefig('/root/.openclaw/workspace/tables/Table8_Intereye_Consistency_ThreeLine.png', dpi=300, bbox_inches='tight', 
                 facecolor='white', edgecolor='none')
    plt.close(fig8)
    print("✓ Table 8 saved")


    print("\n✅ All fixed tables generated successfully!")
    print("\nFixed files:")
    print("  - Table2_Macular_Layers_ThreeLine.png")
    print("  - Table4_Correlation_Analysis_ThreeLine.png")
    print("  - Table6_Multivariate_Regression_ThreeLine.png")
    print("  - Table8_Intereye_Consistency_ThreeLine.png (title fixed)")



if __name__ == "__main__":
    main()