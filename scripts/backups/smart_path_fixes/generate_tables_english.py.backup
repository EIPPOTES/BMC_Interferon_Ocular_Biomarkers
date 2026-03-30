import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import mannwhitneyu, chi2_contingency, spearmanr
from sklearn.metrics import roc_curve, auc
import warnings

def main():
    """主函数，包装原有执行代码"""
    warnings.filterwarnings('ignore')

    # 设置样式
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False

    print("=" * 80)
    print("Generating Publication-Ready Tables (All English Labels)")
    print("=" * 80)

    # 数据路径
    data_dir = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据'
    tables_dir = '/root/.openclaw/workspace/tables'

    # 读取数据
    df = pd.read_excel(f'{data_dir}/OCT数据_完整整合.xlsx')
    print(f"Data loaded: {len(df)} rows")

    dep_df = df[df['分组'] == '抑郁症'].copy()
    ctrl_df = df[df['分组'] == '健康对照'].copy()

    # ==================== TABLE 1: Baseline Characteristics ====================
    print("\n【Table 1: Baseline Characteristics】")

    # 创建表格数据
    table1_data = []

    # Sample size
    table1_data.append(['Sample size, n (eyes)', f"{len(dep_df)} ({len(dep_df)})", 
                        f"{len(ctrl_df)} ({len(ctrl_df)})", ''])

    # Age
    age_dep = dep_df['年龄'].dropna()
    age_ctrl = ctrl_df['年龄'].dropna()
    _, p_age = mannwhitneyu(age_dep, age_ctrl)
    table1_data.append(['Age, years, mean ± SD', f"{age_dep.mean():.1f} ± {age_dep.std():.1f}",
                        f"{age_ctrl.mean():.1f} ± {age_ctrl.std():.1f}", f"{p_age:.4f}"])

    # Sex
    sex_dep = dep_df['性别'].value_counts()
    sex_ctrl = ctrl_df['性别'].value_counts()
    sex_table = pd.crosstab(df['分组'], df['性别'])
    _, p_sex, _, _ = chi2_contingency(sex_table)
    male_dep = sex_dep.get('男', 0)
    female_dep = sex_dep.get('女', 0)
    male_ctrl = sex_ctrl.get('男', 0)
    female_ctrl = sex_ctrl.get('女', 0)
    table1_data.append(['Sex, Male/Female, n', f"{male_dep}/{female_dep}",
                        f"{male_ctrl}/{female_ctrl}", f"{p_sex:.4f}"])

    # PHQ-9
    phq9_dep = dep_df['PHQ-9'].dropna()
    table1_data.append(['PHQ-9 score, mean ± SD', f"{phq9_dep.mean():.2f} ± {phq9_dep.std():.2f}",
                        '-', '-'])

    # GAD-7
    gad7_dep = dep_df['GAD-7'].dropna()
    table1_data.append(['GAD-7 score, mean ± SD', f"{gad7_dep.mean():.2f} ± {gad7_dep.std():.2f}",
                        '-', '-'])

    # Mean Macular Thickness
    mt_dep = dep_df['Retina_平均厚度'].dropna()
    mt_ctrl = ctrl_df['Retina_平均厚度'].dropna()
    _, p_mt = mannwhitneyu(mt_dep, mt_ctrl)
    table1_data.append(['Mean Macular Thickness, μm, mean ± SD', f"{mt_dep.mean():.2f} ± {mt_dep.std():.2f}",
                        f"{mt_ctrl.mean():.2f} ± {mt_ctrl.std():.2f}", f"{p_mt:.4f}"])

    # Total Macular Volume
    mv_dep = dep_df['Retina_总体积'].dropna()
    mv_ctrl = ctrl_df['Retina_总体积'].dropna()
    _, p_mv = mannwhitneyu(mv_dep, mv_ctrl)
    table1_data.append(['Total Macular Volume, mm³, mean ± SD', f"{mv_dep.mean():.2f} ± {mv_dep.std():.2f}",
                        f"{mv_ctrl.mean():.2f} ± {mv_ctrl.std():.2f}", f"{p_mv:.4f}"])

    # RNFL Total
    rnfl_dep = dep_df['RNFL_Total'].dropna()
    rnfl_ctrl = ctrl_df['RNFL_Total'].dropna()
    _, p_rnfl = mannwhitneyu(rnfl_dep, rnfl_ctrl)
    table1_data.append(['RNFL Total, μm, mean ± SD', f"{rnfl_dep.mean():.2f} ± {rnfl_dep.std():.2f}",
                        f"{rnfl_ctrl.mean():.2f} ± {rnfl_ctrl.std():.2f}", f"{p_rnfl:.4f}"])

    # Create DataFrame
    table1 = pd.DataFrame(table1_data, columns=['Characteristic', 'MDD Group', 'Control Group', 'P-value'])

    # Save as image
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('tight')
    ax.axis('off')

    # Title
    table_title = 'Table 1. Baseline Characteristics of Study Participants'
    ax.text(0.5, 0.95, table_title, fontsize=14, fontweight='bold', ha='center', transform=ax.transAxes)

    # Create table
    table = ax.table(cellText=table1.values, colLabels=table1.columns, cellLoc='left',
                     loc='center', bbox=[0, 0.1, 1, 0.8])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style header
    for i in range(len(table1.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    plt.savefig(f'{tables_dir}/Table1_Baseline_Characteristics.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Table 1 saved")

    # ==================== TABLE 2: Macular Layer Analysis ====================
    print("\n【Table 2: Macular Layer Analysis】")

    # 读取详细数据
    table2_df = pd.read_excel(f'{data_dir}/Table2_黄斑所有层_Statsmodels.xlsx')

    # 简化表格 - 只保留关键列 (使用原始中文列名)
    table2_display = table2_df[['层', '区域', '抑郁症组', '对照组', 
                                 'Cohen_d', 'P值', 'P值_FDR', '显著性']].head(20)

    # Rename columns
    table2_display.columns = ['Layer', 'Region', 'MDD', 'Control', "Cohen's d", 'P-value', 'P-FDR', 'Sig']

    # Format numbers (only for numeric columns)
    for col in ["Cohen's d"]:
        table2_display[col] = table2_display[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) and isinstance(x, (int, float)) else str(x) if pd.notna(x) else '-')

    for col in ['P-value', 'P-FDR']:
        table2_display[col] = table2_display[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else '-')

    # Save as image
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.axis('tight')
    ax.axis('off')

    ax.text(0.5, 0.98, 'Table 2. Macular Layer Analysis', fontsize=14, fontweight='bold', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.95, 'Comparison of Macular Thickness Between MDD Patients and Controls', 
            fontsize=11, ha='center', transform=ax.transAxes, style='italic')

    table = ax.table(cellText=table2_display.values, colLabels=table2_display.columns, 
                     cellLoc='center', loc='center', bbox=[0, 0.02, 1, 0.9])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    # Style header
    for i in range(len(table2_display.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    # Highlight significant rows
    for i in range(1, len(table2_display) + 1):
        if table2_display.iloc[i-1]['Sig'] == True or str(table2_display.iloc[i-1]['Sig']).lower() == 'true':
            for j in range(len(table2_display.columns)):
                table[(i, j)].set_facecolor('#FFF2CC')

    plt.savefig(f'{tables_dir}/Table2_Macular_Layers.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Table 2 saved")

    # ==================== TABLE 4: Optic Disc Parameters ====================
    print("\n【Table 4: Optic Disc Parameters】")

    table4_df = pd.read_excel(f'{data_dir}/Table4_视盘所有指标组间比较.xlsx')

    # Select and rename columns (使用原始中文列名)
    table4_display = table4_df[['指标', '抑郁症组', '对照组', 'P值', 
                                'Cohen_d', 'P值_FDR', '显著性']].head(15)
    table4_display.columns = ['Parameter', 'MDD (Mean±SD)', 'Control (Mean±SD)', 'P-value', "Cohen's d", 'P-FDR', 'Sig']

    # Format
    for col in ['P-value', 'P-FDR']:
        if col in table4_display.columns:
            table4_display[col] = table4_display[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else '-')

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('tight')
    ax.axis('off')

    ax.text(0.5, 0.98, 'Table 4. Optic Disc Parameters', fontsize=14, fontweight='bold', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.95, 'Comparison of Optic Disc Measurements Between Groups', 
            fontsize=11, ha='center', transform=ax.transAxes, style='italic')

    table = ax.table(cellText=table4_display.values, colLabels=table4_display.columns, 
                     cellLoc='center', loc='center', bbox=[0, 0.05, 1, 0.87])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    for i in range(len(table4_display.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    plt.savefig(f'{tables_dir}/Table4_Optic_Disc.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Table 4 saved")

    # ==================== TABLE 6: ROC Analysis ====================
    print("\n【Table 6: ROC Analysis】")

    table6_df = pd.read_excel(f'{data_dir}/Table6_全面ROC分析.xlsx')

    # Select top parameters by AUC (使用原始中文列名)
    table6_display = table6_df[['指标', 'AUC', '敏感度', '特异度', 
                                '最佳截断值']].head(15)
    table6_display.columns = ['Parameter', 'AUC', 'Sensitivity', 'Specificity', 'Cutoff']

    # Format
    for col in ['AUC', 'Sensitivity', 'Specificity']:
        table6_display[col] = table6_display[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) and isinstance(x, (int, float)) else str(x) if pd.notna(x) else '-')

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('tight')
    ax.axis('off')

    ax.text(0.5, 0.98, 'Table 6. ROC Analysis', fontsize=14, fontweight='bold', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.95, 'Diagnostic Performance of OCT Parameters', 
            fontsize=11, ha='center', transform=ax.transAxes, style='italic')

    table = ax.table(cellText=table6_display.values, colLabels=table6_display.columns, 
                     cellLoc='center', loc='center', bbox=[0, 0.05, 1, 0.87])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    for i in range(len(table6_display.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    plt.savefig(f'{tables_dir}/Table6_ROC_Analysis.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Table 6 saved")

    # ==================== TABLE 7: Subgroup Analysis ====================
    print("\n【Table 7: Subgroup Analysis】")

    table7_df = pd.read_excel(f'{data_dir}/Table7_亚组分析_严重程度.xlsx')

    # 使用原始中文列名
    table7_display = table7_df[['指标', '无抑郁(n=103)', '轻度(n=54)', '中度(n=40)', '重度(n=63)', 'Kruskal_P', '趋势_P']].head(12)
    table7_display.columns = ['Parameter', 'Minimal (n=103)', 'Mild (n=54)', 'Moderate (n=40)', 'Severe (n=63)', 'P-value', 'Trend P']

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('tight')
    ax.axis('off')

    ax.text(0.5, 0.98, 'Table 7. Subgroup Analysis by Depression Severity', fontsize=14, fontweight='bold', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.95, 'OCT Parameters Across PHQ-9 Severity Groups', 
            fontsize=11, ha='center', transform=ax.transAxes, style='italic')

    table = ax.table(cellText=table7_display.values, colLabels=table7_display.columns, 
                     cellLoc='center', loc='center', bbox=[0, 0.05, 1, 0.87])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    for i in range(len(table7_display.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    plt.savefig(f'{tables_dir}/Table7_Subgroup_Analysis.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Table 7 saved")

    # ==================== TABLE 9: Inter-eye Consistency ====================
    print("\n【Table 9: Inter-eye Consistency】")

    table9_df = pd.read_excel(f'{data_dir}/Table9_双眼一致性分析.xlsx')

    # 使用原始中文列名
    table9_display = table9_df[['指标', 'Pearson_r', 'P值', '平均绝对差异_μm', 'ICC估计']].head(12)
    table9_display.columns = ['Parameter', 'Pearson r', 'P-value', 'Mean Diff (μm)', 'ICC']

    # Format
    for col in ['Pearson r', 'P-value', 'Mean Diff (μm)', 'ICC']:
        table9_display[col] = table9_display[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) and isinstance(x, (int, float)) else str(x) if pd.notna(x) else '-')

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('tight')
    ax.axis('off')

    ax.text(0.5, 0.98, 'Table 9. Inter-eye Consistency Analysis', fontsize=14, fontweight='bold', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.95, 'Correlation Between Right and Left Eye Measurements', 
            fontsize=11, ha='center', transform=ax.transAxes, style='italic')

    table = ax.table(cellText=table9_display.values, colLabels=table9_display.columns, 
                     cellLoc='center', loc='center', bbox=[0, 0.1, 1, 0.82])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    for i in range(len(table9_display.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    plt.savefig(f'{tables_dir}/Table9_Intereye_Consistency.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Table 9 saved")

    # ==================== TABLE 11: Feature Importance ====================
    print("\n【Table 11: Feature Importance】")

    table11_df = pd.read_excel(f'{data_dir}/Table11_特征重要性.xlsx')

    # 使用原始中文列名
    table11_display = table11_df[['特征', '绝对系数']].head(15)
    table11_display['Rank'] = range(1, len(table11_display) + 1)
    table11_display = table11_display[['特征', '绝对系数', 'Rank']]
    table11_display.columns = ['Feature', 'Importance Score', 'Rank']

    # Format
    table11_display['Importance Score'] = table11_display['Importance Score'].apply(lambda x: f"{x:.4f}" if pd.notna(x) else '-')

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.axis('tight')
    ax.axis('off')

    ax.text(0.5, 0.98, 'Table 11. Feature Importance', fontsize=14, fontweight='bold', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.95, 'Top Predictive Features from Random Forest Model', 
            fontsize=11, ha='center', transform=ax.transAxes, style='italic')

    table = ax.table(cellText=table11_display.values, colLabels=table11_display.columns, 
                     cellLoc='center', loc='center', bbox=[0, 0.05, 1, 0.87])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.8)

    for i in range(len(table11_display.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    plt.savefig(f'{tables_dir}/Table11_Feature_Importance.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ Table 11 saved")

    print("\n" + "=" * 80)
    print("All 7 Tables generated successfully!")
    print(f"Output directory: {tables_dir}")
    print("=" * 80)



if __name__ == "__main__":
    main()