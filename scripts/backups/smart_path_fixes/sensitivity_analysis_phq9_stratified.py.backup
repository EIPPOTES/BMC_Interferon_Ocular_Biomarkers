import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import warnings

def main():
    """主函数，包装原有执行代码"""
    warnings.filterwarnings('ignore')

    print("=" * 100)
    print("Sensitivity Analysis: MDD Patients - Active vs Remission")
    print("=" * 100)

    # 读取数据
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx')

    # 创建分组标签
    df['Group'] = df['分组']  # 保留原始分组
    df['Label'] = (df['分组'] == '抑郁症').astype(int)

    print(f"\n总样本: {len(df)} 眼")
    print(f"MDD: {df['Label'].sum()} 眼")
    print(f"对照: {len(df) - df['Label'].sum()} 眼")

    # 筛选有PHQ-9数据的MDD患者
    df_mdd_phq = df[(df['Label'] == 1) & (df['PHQ-9'].notna())].copy()
    print(f"\n有PHQ-9数据的MDD患者: {len(df_mdd_phq)} 眼")

    # 按PHQ-9分层
    df_mdd_active = df_mdd_phq[df_mdd_phq['PHQ-9'] >= 5].copy()  # 活跃抑郁
    df_mdd_remission = df_mdd_phq[df_mdd_phq['PHQ-9'] < 5].copy()  # 缓解期

    print(f"  活跃抑郁 (PHQ-9 >= 5): {len(df_mdd_active)} 眼 ({len(df_mdd_active)/len(df_mdd_phq)*100:.1f}%)")
    print(f"  缓解期 (PHQ-9 < 5): {len(df_mdd_remission)} 眼 ({len(df_mdd_remission)/len(df_mdd_phq)*100:.1f}%)")

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
    print("Analysis 1: Active Depression (PHQ-9 >= 5) vs Controls")
    print("=" * 100)

    for col, name in key_outcomes:
        print(f"\n{'='*80}")
        print(f"Outcome: {name}")
        print(f"{'='*80}")

        # 分析1: 活跃抑郁 vs 对照
        df_analysis1 = pd.concat([
            df_mdd_active[['Label', '年龄', '性别', col]],
            df[df['Label'] == 0][['Label', '年龄', '性别', col]]
        ]).dropna()

        if len(df_analysis1) >= 30:
            model1 = smf.ols(f'Q("{col}") ~ Label + 年龄 + C(性别)', data=df_analysis1).fit()
            beta1 = model1.params['Label']
            p1 = model1.pvalues['Label']
            n1 = len(df_analysis1)
            n_mdd1 = df_analysis1['Label'].sum()
            n_ctrl1 = len(df_analysis1) - n_mdd1
        else:
            beta1 = np.nan
            p1 = np.nan
            n1 = n_mdd1 = n_ctrl1 = 0

        # 分析2: 缓解期 vs 对照
        df_analysis2 = pd.concat([
            df_mdd_remission[['Label', '年龄', '性别', col]],
            df[df['Label'] == 0][['Label', '年龄', '性别', col]]
        ]).dropna()

        if len(df_analysis2) >= 30:
            model2 = smf.ols(f'Q("{col}") ~ Label + 年龄 + C(性别)', data=df_analysis2).fit()
            beta2 = model2.params['Label']
            p2 = model2.pvalues['Label']
            n2 = len(df_analysis2)
            n_mdd2 = df_analysis2['Label'].sum()
            n_ctrl2 = len(df_analysis2) - n_mdd2
        else:
            beta2 = np.nan
            p2 = np.nan
            n2 = n_mdd2 = n_ctrl2 = 0

        # 分析3: 所有MDD vs 对照（原始分析）
        df_analysis3 = pd.concat([
            df_mdd_phq[['Label', '年龄', '性别', col]],
            df[df['Label'] == 0][['Label', '年龄', '性别', col]]
        ]).dropna()

        if len(df_analysis3) >= 30:
            model3 = smf.ols(f'Q("{col}") ~ Label + 年龄 + C(性别)', data=df_analysis3).fit()
            beta3 = model3.params['Label']
            p3 = model3.pvalues['Label']
            n3 = len(df_analysis3)
        else:
            beta3 = np.nan
            p3 = np.nan
            n3 = 0

        print(f"\nActive Depression (PHQ-9 >= 5) vs Controls:")
        print(f"  n = {n1} (MDD: {n_mdd1}, Control: {n_ctrl1})")
        print(f"  β = {beta1:.3f}, P = {p1:.3f}")

        print(f"\nRemission (PHQ-9 < 5) vs Controls:")
        print(f"  n = {n2} (MDD: {n_mdd2}, Control: {n_ctrl2})")
        print(f"  β = {beta2:.3f}, P = {p2:.3f}")

        print(f"\nAll MDD vs Controls (original):")
        print(f"  n = {n3}")
        print(f"  β = {beta3:.3f}, P = {p3:.3f}")

        # 结论
        if p1 < 0.05 and p2 < 0.05:
            conclusion = "Significant in both active and remission groups"
        elif p1 < 0.05 and p2 >= 0.05:
            conclusion = "Significant only in active depression"
        elif p1 >= 0.05 and p2 < 0.05:
            conclusion = "Significant only in remission group"
        else:
            conclusion = "Non-significant in both groups"

        print(f"\nConclusion: {conclusion}")

        results.append({
            'Outcome': name,
            'N_Active': n1,
            'N_Remission': n2,
            'N_All': n3,
            'Beta_Active': round(beta1, 3),
            'P_Active': round(p1, 3),
            'Beta_Remission': round(beta2, 3),
            'P_Remission': round(p2, 3),
            'Beta_All': round(beta3, 3),
            'P_All': round(p3, 3),
            'Conclusion': conclusion
        })

    # 创建结果DataFrame
    results_df = pd.DataFrame(results)

    print("\n" + "=" * 100)
    print("Summary Table: Stratified Analysis by Depression Severity")
    print("=" * 100)
    print(results_df.to_string(index=False))

    # 保存结果
    output_path = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/Sensitivity_Analysis_PHQ9_Stratified.xlsx'
    results_df.to_excel(output_path, index=False)
    print(f"\nResults saved to: {output_path}")

    workspace_path = '/root/.openclaw/workspace/revised_paper/Sensitivity_Analysis_PHQ9_Stratified.xlsx'
    results_df.to_excel(workspace_path, index=False)
    print(f"Results saved to: {workspace_path}")

    # 生成论文用表格
    print("\n" + "=" * 100)
    print("Table for Manuscript: Stratified Analysis by PHQ-9 Severity")
    print("=" * 100)
    print("\n| Outcome | Active vs Control | Remission vs Control | All MDD vs Control |")
    print("| Outcome | β (P) | β (P) | β (P) |")
    print("|---------|-------------------|----------------------|-------------------|")
    for _, row in results_df.iterrows():
        print(f"| {row['Outcome']} | {row['Beta_Active']:.3f} ({row['P_Active']:.3f}) | {row['Beta_Remission']:.3f} ({row['P_Remission']:.3f}) | {row['Beta_All']:.3f} ({row['P_All']:.3f}) |")

    print("\n" + "=" * 100)
    print("Key Findings:")
    print("=" * 100)

    active_sig = results_df[results_df['P_Active'] < 0.05]
    remission_sig = results_df[results_df['P_Remission'] < 0.05]

    print(f"\nSignificant in Active Depression (PHQ-9 >= 5): {len(active_sig)}/{len(results_df)} parameters")
    for _, row in active_sig.iterrows():
        print(f"  - {row['Outcome']}: β = {row['Beta_Active']:.3f}")

    print(f"\nSignificant in Remission (PHQ-9 < 5): {len(remission_sig)}/{len(results_df)} parameters")
    for _, row in remission_sig.iterrows():
        print(f"  - {row['Outcome']}: β = {row['Beta_Remission']:.3f}")

    print("\n" + "=" * 100)
    print("Interpretation:")
    print("=" * 100)
    if len(active_sig) > len(remission_sig):
        print("The association between depression and retinal changes appears stronger in")
        print("patients with active depression (PHQ-9 >= 5) compared to those in remission,")
        print("suggesting a dose-response relationship with depression severity.")
    elif len(active_sig) == len(remission_sig):
        print("The association between depression and retinal changes is consistent across")
        print("both active depression and remission groups, suggesting these may be trait")
        print("markers rather than state-dependent changes.")
    else:
        print("Unexpectedly, some associations were stronger in the remission group,")
        print("suggesting potential complex relationships between depression course and")
        print("retinal structural changes.")
    print("=" * 100)



if __name__ == "__main__":
    main()