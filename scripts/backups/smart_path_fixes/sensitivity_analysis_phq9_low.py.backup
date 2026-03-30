import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import warnings

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
def main():
    """主函数，包装原有执行代码"""
    warnings.filterwarnings('ignore')

    print("=" * 100)
    print("Sensitivity Analysis: Excluding PHQ-9 < 5 Patients")
    print("=" * 100)

    # 读取数据
    df = pd.read_excel(str(/mnt/c/Users/CUI/Desktop/论文及图表))

    # 创建二分类标签
    df['Label'] = (df['分组'] == '抑郁症').astype(int)

    print(f"\n总样本: {len(df)} 眼")
    print(f"MDD: {df['Label'].sum()} 眼")
    print(f"对照: {len(df) - df['Label'].sum()} 眼")

    # 筛选有PHQ-9数据的MDD患者
    df_phq = df[df['PHQ-9'].notna()].copy()
    print(f"\n有PHQ-9数据的样本: {len(df_phq)} 眼")

    # 排除PHQ-9 < 5的患者（保留PHQ-9 >= 5的活跃抑郁患者）
    df_active = df_phq[(df_phq['Label'] == 0) | (df_phq['PHQ-9'] >= 5)].copy()
    print(f"排除PHQ-9 < 5后: {len(df_active)} 眼")
    print(f"  MDD (PHQ-9 >= 5): {df_active['Label'].sum()} 眼")
    print(f"  对照: {len(df_active) - df_active['Label'].sum()} 眼")

    # 定义关键结局变量
    key_outcomes = [
        ('Retina_平均厚度', 'Mean_Macular_Thickness'),
        ('Retina_外环颞侧', 'Outer_Temporal_Thickness'),
        ('Retina_总体积', 'Total_Macular_Volume'),
        ('C/D Area Ratio', 'CD_Ratio'),
        ('Rim Volume', 'Rim_Volume')
    ]

    # 存储结果
    results = []

    print("\n" + "=" * 100)
    print("Comparison: Original Analysis vs Excluding PHQ-9 < 5")
    print("=" * 100)

    for col, name in key_outcomes:
        print(f"\n{'='*80}")
        print(f"Outcome: {name}")
        print(f"{'='*80}")

        # 原始分析（所有有PHQ-9数据的患者）
        valid_original = df_phq[['Label', '年龄', '性别', col]].dropna()
        if len(valid_original) >= 30:
            model_orig = smf.ols(f'Q("{col}") ~ Label + 年龄 + C(性别)', data=valid_original).fit()
            beta_orig = model_orig.params['Label']
            p_orig = model_orig.pvalues['Label']
            n_orig = len(valid_original)
        else:
            beta_orig = np.nan
            p_orig = np.nan
            n_orig = 0

        # 排除PHQ-9 < 5后的分析
        valid_active = df_active[['Label', '年龄', '性别', col]].dropna()
        if len(valid_active) >= 30:
            model_active = smf.ols(f'Q("{col}") ~ Label + 年龄 + C(性别)', data=valid_active).fit()
            beta_active = model_active.params['Label']
            p_active = model_active.pvalues['Label']
            n_active = len(valid_active)
        else:
            beta_active = np.nan
            p_active = np.nan
            n_active = 0

        # 计算变化
        if not np.isnan(beta_orig) and not np.isnan(beta_active):
            beta_change = beta_active - beta_orig
            percent_change = (beta_change / abs(beta_orig)) * 100 if beta_orig != 0 else 0
        else:
            beta_change = np.nan
            percent_change = np.nan

        print(f"\nOriginal Analysis (all PHQ-9 data, n={n_orig}):")
        print(f"  β = {beta_orig:.3f}, P = {p_orig:.3f}")

        print(f"\nExcluding PHQ-9 < 5 (active depression only, n={n_active}):")
        print(f"  β = {beta_active:.3f}, P = {p_active:.3f}")

        print(f"\nChange:")
        print(f"  Δβ = {beta_change:.3f} ({percent_change:+.1f}%)")

        # 结论
        if p_orig < 0.05 and p_active < 0.05:
            conclusion = "Significant in both analyses - robust finding"
        elif p_orig < 0.05 and p_active >= 0.05:
            conclusion = "WARNING: Significance lost after exclusion"
        elif p_orig >= 0.05 and p_active < 0.05:
            conclusion = "Significance emerged after exclusion"
        else:
            conclusion = "Non-significant in both analyses"

        print(f"\nConclusion: {conclusion}")

        results.append({
            'Outcome': name,
            'N_Original': n_orig,
            'N_Active': n_active,
            'Beta_Original': round(beta_orig, 3),
            'P_Original': round(p_orig, 3),
            'Beta_Active': round(beta_active, 3),
            'P_Active': round(p_active, 3),
            'Beta_Change': round(beta_change, 3),
            'Percent_Change': round(percent_change, 1),
            'Conclusion': conclusion
        })

    # 创建结果DataFrame
    results_df = pd.DataFrame(results)

    print("\n" + "=" * 100)
    print("Summary Table: Sensitivity Analysis Results")
    print("=" * 100)
    print(results_df.to_string(index=False))

    # 保存结果
    output_path = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    results_df.to_excel(output_path, index=False)
    print(f"\nResults saved to: {output_path}")

    workspace_path = '/root/.openclaw/workspace/revised_paper/Sensitivity_Analysis_Exclude_PHQ9_Low.xlsx'
    results_df.to_excel(workspace_path, index=False)
    print(f"Results saved to: {workspace_path}")

    # 生成论文用表格
    print("\n" + "=" * 100)
    print("Table for Manuscript: Sensitivity Analysis Excluding PHQ-9 < 5")
    print("=" * 100)
    print("\n| Outcome | N (Original) | N (Active) | β (Original) | P (Original) | β (Active) | P (Active) | Δβ | Conclusion |")
    print("|---------|--------------|------------|--------------|--------------|------------|------------|-----|------------|")
    for _, row in results_df.iterrows():
        print(f"| {row['Outcome']} | {row['N_Original']} | {row['N_Active']} | {row['Beta_Original']:.3f} | {row['P_Original']:.3f} | {row['Beta_Active']:.3f} | {row['P_Active']:.3f} | {row['Beta_Change']:+.3f} | {row['Conclusion']} |")

    print("\n" + "=" * 100)
    print("Key Findings:")
    print("=" * 100)

    robust = results_df[results_df['Conclusion'].str.contains('robust')]
    if len(robust) > 0:
        print(f"\nRobust findings (significant after excluding PHQ-9 < 5):")
        for _, row in robust.iterrows():
            print(f"  - {row['Outcome']}: β changed by {row['Percent_Change']:+.1f}%")

    warning = results_df[results_df['Conclusion'].str.contains('WARNING')]
    if len(warning) > 0:
        print(f"\nWARNING: Findings that lost significance:")
        for _, row in warning.iterrows():
            print(f"  - {row['Outcome']}")

    print("\n" + "=" * 100)
    print("Interpretation:")
    print("=" * 100)
    if len(warning) == 0:
        print("Excluding patients with PHQ-9 < 5 (minimal symptoms) did not substantially")
        print("change the main findings. Group differences remained significant for key OCT")
        print("parameters, supporting the robustness of our findings across different depression")
        print("severity levels.")
    else:
        print("Some findings were sensitive to excluding patients with minimal symptoms,")
        print("suggesting that depression severity may modify the association between MDD")
        print("and retinal structural changes.")
    print("=" * 100)



if __name__ == "__main__":
    main()