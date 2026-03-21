import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import numpy as np

# 强制重新初始化字体管理器

def main():
    """主函数，包装原有执行代码"""
    fm.fontManager.__init__()

    # 适配中文字体
    cjk_list = ['CJK', 'Han', 'CN', 'TW', 'JP']
    cjk_fonts = [f.name for f in fm.fontManager.ttflist if any(s.lower() in f.name.lower() for s in cjk_list)]

    plt.rcParams['font.family'] = ['DejaVu Sans'] + cjk_fonts
    plt.rcParams['axes.unicode_minus'] = False

    # Table 7 数据 - 亚组分析
    subgroup_data = {
        'Parameter': [
            'n (eyes)',
            'Mean Macular Thickness (μm)',
            'Total Macular Volume (mm³)',
            'Outer Temporal (μm)',
            'Inner Temporal (μm)',
            'Outer Superior (μm)',
            'Inner Superior (μm)',
            'Outer Inferior (μm)',
            'Inner Inferior (μm)',
            'Outer Nasal (μm)',
            'Inner Nasal (μm)',
            'Central (μm)',
            'RNFL Total (μm)',
            'C/D Area Ratio',
            'Rim Volume (mm³)'
        ],
        'Minimal\n(PHQ-9<5)': [
            '103',
            '272.34 ± 16.12',
            '7.71 ± 0.46',
            '266.70 ± 18.23',
            '269.12 ± 17.45',
            '275.45 ± 16.78',
            '278.89 ± 16.12',
            '274.23 ± 17.45',
            '277.56 ± 16.89',
            '281.12 ± 15.67',
            '283.34 ± 15.34',
            '258.67 ± 19.23',
            '106.34 ± 16.78',
            '0.29 ± 0.18',
            '0.26 ± 0.15'
        ],
        'Mild\n(PHQ-9 5-9)': [
            '54',
            '270.89 ± 17.45',
            '7.68 ± 0.49',
            '268.99 ± 17.12',
            '270.34 ± 16.23',
            '276.89 ± 17.34',
            '279.45 ± 15.67',
            '275.12 ± 16.78',
            '278.23 ± 16.12',
            '282.45 ± 16.23',
            '284.12 ± 14.89',
            '260.12 ± 18.45',
            '105.12 ± 17.34',
            '0.30 ± 0.19',
            '0.25 ± 0.16'
        ],
        'Moderate\n(PHQ-9 10-14)': [
            '40',
            '273.12 ± 15.67',
            '7.74 ± 0.44',
            '274.45 ± 16.78',
            '273.67 ± 16.12',
            '279.12 ± 15.67',
            '281.23 ± 14.89',
            '277.45 ± 15.34',
            '280.12 ± 15.67',
            '284.78 ± 14.56',
            '286.45 ± 14.23',
            '261.45 ± 17.67',
            '107.81 ± 15.89',
            '0.31 ± 0.20',
            '0.24 ± 0.14'
        ],
        'Severe\n(PHQ-9≥15)': [
            '63',
            '271.54 ± 16.89',
            '7.68 ± 0.47',
            '274.12 ± 17.34',
            '272.45 ± 17.89',
            '278.56 ± 16.45',
            '280.67 ± 15.78',
            '276.34 ± 16.89',
            '279.45 ± 16.23',
            '283.67 ± 15.78',
            '285.23 ± 15.12',
            '260.89 ± 18.12',
            '107.23 ± 16.45',
            '0.32 ± 0.21',
            '0.23 ± 0.15'
        ],
        'P-trend': [
            '-',
            '0.502',
            '0.512',
            '0.020',
            '0.089',
            '0.156',
            '0.234',
            '0.189',
            '0.267',
            '0.312',
            '0.289',
            '0.445',
            '0.631',
            '0.212',
            '0.089'
        ]
    }

    df = pd.DataFrame(subgroup_data)

    # 创建图形 - 使用text方式绘制标准三线表
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 20)
    ax.axis('off')

    # 标题
    ax.text(8, 19, 'Table 7. OCT Parameters by Depression Severity (n=260 eyes from 132 patients)', 
            ha='center', va='center', fontsize=13, fontweight='bold')

    # 列宽设置
    col_x = [0.3, 3.0, 5.8, 8.6, 11.4, 14.2]
    col_widths = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5]

    # 表头
    headers = ['Parameter', 'Minimal\n(PHQ-9<5)', 'Mild\n(PHQ-9 5-9)', 'Moderate\n(PHQ-9 10-14)', 'Severe\n(PHQ-9≥15)', 'P-trend']
    for i, (header, x) in enumerate(zip(headers, col_x)):
        ax.text(x + col_widths[i]/2, 17.5, header, ha='center', va='center', 
                fontsize=10, fontweight='bold')

    # 顶部横线
    ax.plot([0.1, 15.9], [18.2, 18.2], 'k-', linewidth=1.5)
    # 表头下方横线
    ax.plot([0.1, 15.9], [16.8, 16.8], 'k-', linewidth=1.0)

    # 数据行 - 减小行间距
    row_height = 0.85
    start_y = 15.9
    for idx, row in df.iterrows():
        y = start_y - idx * row_height
        for i, (col, x) in enumerate(zip(df.columns, col_x)):
            ax.text(x + col_widths[i]/2, y, row[col], ha='center', va='center', fontsize=9)

    # 底部横线
    end_y = start_y - (len(df) - 1) * row_height - 0.4
    ax.plot([0.1, 15.9], [end_y, end_y], 'k-', linewidth=1.5)

    # 添加注释
    ax.text(8, end_y - 0.6, 
            'Data are mean ± SD. P-trend from Jonckheere-Terpstra test.',
            ha='center', va='center', fontsize=9, style='italic')

    plt.tight_layout()
    plt.savefig('/mnt/c/Users/CUI/Desktop/最终版/02_图表_Table/Table7_Subgroup_Analysis_ThreeLine.png', 
                dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
    plt.savefig('/root/.openclaw/workspace/tables/Table7_Subgroup_Analysis_ThreeLine.png', 
                dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
    print("Table 7 (fixed row spacing) generated successfully!")



if __name__ == "__main__":
    main()