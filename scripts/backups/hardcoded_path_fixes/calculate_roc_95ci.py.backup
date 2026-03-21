import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
from scipy import stats
import warnings

def main():
    """主函数，包装原有执行代码"""
    warnings.filterwarnings('ignore')

    # 读取数据
    df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/03_OCT数据_最终清洗.xlsx')

    print("=" * 80)
    print("ROC Analysis with 95% Confidence Intervals")
    print("=" * 80)

    # 创建二分类标签 (抑郁症=1, 对照=0)
    df['Label'] = (df['分组'] == '抑郁症').astype(int)

    # 定义关键指标
    key_indicators = [
        ('Retina_平均厚度', 'Mean Macular Thickness'),
        ('Retina_外环颞侧', 'Outer Temporal Thickness'),
        ('Retina_总体积', 'Total Macular Volume'),
        ('RNFL_Total', 'RNFL Total'),
        ('RNFL_上方', 'RNFL Superior'),
        ('C/D Area Ratio', 'Cup-to-Disc Area Ratio'),
        ('Rim Volume', 'Rim Volume')
    ]

    # Bootstrap函数计算AUC的95% CI
    def bootstrap_auc_ci(y_true, y_score, n_bootstraps=1000, ci=0.95, random_state=42):
        """
        使用Bootstrap计算AUC的置信区间
        """
        rng = np.random.RandomState(random_state)
        bootstrapped_aucs = []

        for i in range(n_bootstraps):
            # 有放回抽样
            indices = rng.randint(0, len(y_true), len(y_true))
            if len(np.unique(y_true[indices])) < 2:
                # 如果抽样后只有一个类别，跳过
                continue

            fpr, tpr, _ = roc_curve(y_true[indices], y_score[indices])
            auc_score = auc(fpr, tpr)
            bootstrapped_aucs.append(auc_score)

        # 计算置信区间
        alpha = (1 - ci) / 2
        ci_lower = np.percentile(bootstrapped_aucs, alpha * 100)
        ci_upper = np.percentile(bootstrapped_aucs, (1 - alpha) * 100)

        return ci_lower, ci_upper, bootstrapped_aucs

    # DeLong检验函数（简化版，用于比较两个AUC）
    def delong_auc_variance(y_true, y_score):
        """
        计算AUC的方差（用于DeLong检验）
        """
        n1 = np.sum(y_true == 1)
        n0 = np.sum(y_true == 0)

        # 计算每个样本的rank
        ranks = stats.rankdata(y_score)

        # 计算V10和V01
        ranks_positive = ranks[y_true == 1]
        ranks_negative = ranks[y_true == 0]

        theta = auc(*roc_curve(y_true, y_score)[:2])

        Q1 = theta / (2 - theta)
        Q2 = 2 * theta**2 / (1 + theta)

        var = (theta * (1 - theta) + (n1 - 1) * (Q1 - theta**2) + (n0 - 1) * (Q2 - theta**2)) / (n1 * n0)

        return var

    # 存储结果
    results = []

    print(f"\nSample size: {len(df)} eyes")
    print(f"MDD: {np.sum(df['Label'] == 1)} eyes")
    print(f"Control: {np.sum(df['Label'] == 0)} eyes")
    print("\n" + "=" * 80)
    print(f"{'Indicator':<35} {'AUC':<8} {'95% CI Lower':<12} {'95% CI Upper':<12} {'Interpretation'}")
    print("=" * 80)

    for col, name in key_indicators:
        if col not in df.columns:
            print(f"Warning: {col} not found in data")
            continue

        # 删除缺失值
        valid_data = df[[col, 'Label']].dropna()

        if len(valid_data) < 10:
            print(f"Warning: Insufficient data for {col}")
            continue

        y_true = valid_data['Label'].values
        y_score = valid_data[col].values

        # 计算AUC (注意：如果AUC < 0.5，说明方向反了，需要反转)
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        auc_value = auc(fpr, tpr)

        # 如果AUC < 0.5，反转预测方向
        if auc_value < 0.5:
            auc_value = 1 - auc_value

        # 计算Bootstrap 95% CI (如果AUC<0.5，需要反转)
        if auc(fpr, tpr) < 0.5:
            ci_lower, ci_upper, _ = bootstrap_auc_ci(y_true, -y_score, n_bootstraps=2000)
        else:
            ci_lower, ci_upper, _ = bootstrap_auc_ci(y_true, y_score, n_bootstraps=2000)

        # 解释
        if auc_value < 0.6:
            interpretation = "Poor"
        elif auc_value < 0.7:
            interpretation = "Fair"
        elif auc_value < 0.8:
            interpretation = "Good"
        elif auc_value < 0.9:
            interpretation = "Very Good"
        else:
            interpretation = "Excellent"

        results.append({
            'Indicator': name,
            'Column': col,
            'AUC': auc_value,
            'CI_Lower': ci_lower,
            'CI_Upper': ci_upper,
            'Interpretation': interpretation,
            'N': len(valid_data)
        })

        print(f"{name:<35} {auc_value:.3f}    {ci_lower:.3f}       {ci_upper:.3f}       {interpretation}")

    print("=" * 80)

    # 创建结果DataFrame
    results_df = pd.DataFrame(results)

    # 保存结果
    output_path = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/ROC_Analysis_with_95CI.xlsx'
    results_df.to_excel(output_path, index=False)
    print(f"\nResults saved to: {output_path}")

    # 同时保存到workspace
    workspace_path = '/root/.openclaw/workspace/revised_paper/ROC_Analysis_with_95CI.xlsx'
    results_df.to_excel(workspace_path, index=False)
    print(f"Results saved to: {workspace_path}")

    # 打印详细结果
    print("\n" + "=" * 80)
    print("Detailed Results for Manuscript")
    print("=" * 80)

    for _, row in results_df.iterrows():
        print(f"\n{row['Indicator']}:")
        print(f"  AUC = {row['AUC']:.3f} (95% CI: {row['CI_Lower']:.3f}–{row['CI_Upper']:.3f})")
        print(f"  N = {row['N']} eyes")
        print(f"  Diagnostic accuracy: {row['Interpretation']}")

    print("\n" + "=" * 80)
    print("Note: 95% CI calculated using Bootstrap method (2000 iterations)")
    print("=" * 80)



if __name__ == "__main__":
    main()