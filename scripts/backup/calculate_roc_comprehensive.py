import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 读取数据
df = pd.read_excel('/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/03_OCT数据_最终清洗.xlsx')

print("=" * 100)
print("Comprehensive ROC Analysis with 95% Confidence Intervals")
print("=" * 100)

# 创建二分类标签 (抑郁症=1, 对照=0)
df['Label'] = (df['分组'] == '抑郁症').astype(int)

# 定义关键指标（与论文Table 5一致）
key_indicators = [
    ('Retina_外环颞侧', 'Outer Temporal Thickness'),
    ('Retina_内环颞侧', 'Inner Temporal Thickness'),
    ('Retina_外环上方', 'Outer Superior Thickness'),
    ('Retina_平均厚度', 'Mean Macular Thickness'),
    ('Retina_总体积', 'Total Macular Volume'),
    ('C/D Area Ratio', 'Cup-to-Disc Area Ratio'),
    ('Rim Volume', 'Rim Volume'),
    ('RNFL_Total', 'RNFL Total'),
    ('RNFL_上方', 'RNFL Superior')
]

# Bootstrap函数
def bootstrap_metric_ci(y_true, y_pred, metric_func, n_bootstraps=1000, ci=0.95, random_state=42):
    """
    使用Bootstrap计算指标的置信区间
    """
    rng = np.random.RandomState(random_state)
    bootstrapped_metrics = []
    
    for i in range(n_bootstraps):
        indices = rng.randint(0, len(y_true), len(y_true))
        if len(np.unique(y_true[indices])) < 2:
            continue
        metric = metric_func(y_true[indices], y_pred[indices])
        bootstrapped_metrics.append(metric)
    
    alpha = (1 - ci) / 2
    ci_lower = np.percentile(bootstrapped_metrics, alpha * 100)
    ci_upper = np.percentile(bootstrapped_metrics, (1 - alpha) * 100)
    
    return ci_lower, ci_upper

# 计算敏感度和特异度
def calculate_sensitivity_specificity(y_true, y_pred):
    tp = np.sum((y_true == 1) & (y_pred == 1))
    tn = np.sum((y_true == 0) & (y_pred == 0))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    fn = np.sum((y_true == 1) & (y_pred == 0))
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    return sensitivity, specificity

# 存储结果
results = []

print(f"\nSample size: {len(df)} eyes")
print(f"MDD: {np.sum(df['Label'] == 1)} eyes")
print(f"Control: {np.sum(df['Label'] == 0)} eyes")
print("\n" + "=" * 100)
print(f"{'Indicator':<30} {'AUC':<10} {'95% CI':<20} {'Sens':<10} {'Spec':<10} {'Cutoff'}")
print("=" * 100)

for col, name in key_indicators:
    if col not in df.columns:
        continue
    
    valid_data = df[[col, 'Label']].dropna()
    
    if len(valid_data) < 10:
        continue
    
    y_true = valid_data['Label'].values
    y_score = valid_data[col].values
    
    # 计算ROC
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_value = auc(fpr, tpr)
    
    # 如果AUC < 0.5，反转
    if auc_value < 0.5:
        auc_value = 1 - auc_value
        y_score_for_ci = -y_score
    else:
        y_score_for_ci = y_score
    
    # 计算Bootstrap 95% CI for AUC
    def auc_metric(y_t, y_s):
        fpr, tpr, _ = roc_curve(y_t, y_s)
        return auc(fpr, tpr) if auc(fpr, tpr) >= 0.5 else 1 - auc(fpr, tpr)
    
    auc_ci_lower, auc_ci_upper = bootstrap_metric_ci(y_true, y_score_for_ci, auc_metric, n_bootstraps=2000)
    
    # 找到最佳截断点 (Youden指数)
    youden_index = tpr - fpr
    if auc(fpr, tpr) < 0.5:
        youden_index = (1-fpr) - (1-tpr)
    optimal_idx = np.argmax(youden_index)
    optimal_threshold = thresholds[optimal_idx] if auc(fpr, tpr) >= 0.5 else -thresholds[optimal_idx]
    
    # 计算敏感度和特异度
    if auc(fpr, tpr) < 0.5:
        y_pred = (y_score < optimal_threshold).astype(int)
    else:
        y_pred = (y_score >= optimal_threshold).astype(int)
    
    sensitivity, specificity = calculate_sensitivity_specificity(y_true, y_pred)
    
    results.append({
        'Indicator': name,
        'AUC': round(auc_value, 3),
        'AUC_CI_Lower': round(auc_ci_lower, 3),
        'AUC_CI_Upper': round(auc_ci_upper, 3),
        'Sensitivity': round(sensitivity, 3),
        'Specificity': round(specificity, 3),
        'Optimal_Cutoff': round(optimal_threshold, 2),
        'N': len(valid_data)
    })
    
    print(f"{name:<30} {auc_value:.3f}      [{auc_ci_lower:.3f}–{auc_ci_upper:.3f}]{'' :<6} {sensitivity:.3f}      {specificity:.3f}      {optimal_threshold:.2f}")

print("=" * 100)

# 创建结果DataFrame
results_df = pd.DataFrame(results)

# 保存结果
output_path = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/ROC_Analysis_Comprehensive_with_95CI.xlsx'
results_df.to_excel(output_path, index=False)
print(f"\nResults saved to: {output_path}")

# 同时保存到workspace
workspace_path = '/root/.openclaw/workspace/revised_paper/ROC_Analysis_Comprehensive_with_95CI.xlsx'
results_df.to_excel(workspace_path, index=False)
print(f"Results saved to: {workspace_path}")

# 打印用于论文的格式
print("\n" + "=" * 100)
print("Results Formatted for Manuscript (Table 5 Format)")
print("=" * 100)

for _, row in results_df.iterrows():
    print(f"\n{row['Indicator']}:")
    print(f"  AUC = {row['AUC']:.3f} (95% CI: {row['AUC_CI_Lower']:.3f}–{row['AUC_CI_Upper']:.3f})")
    print(f"  Sensitivity = {row['Sensitivity']:.1%}")
    print(f"  Specificity = {row['Specificity']:.1%}")
    print(f"  Optimal cutoff = {row['Optimal_Cutoff']:.2f}")

# 创建汇总表格
print("\n" + "=" * 100)
print("Summary Table for Manuscript")
print("=" * 100)
print("\n| OCT Parameter | AUC (95% CI) | Sensitivity | Specificity | Optimal Cutoff |")
print("|---------------|--------------|-------------|-------------|----------------|")
for _, row in results_df.iterrows():
    print(f"| {row['Indicator']} | {row['AUC']:.3f} ({row['AUC_CI_Lower']:.3f}–{row['AUC_CI_Upper']:.3f}) | {row['Sensitivity']:.1%} | {row['Specificity']:.1%} | {row['Optimal_Cutoff']:.2f} |")

print("\n" + "=" * 100)
print("Note:")
print("- 95% CI calculated using Bootstrap method (2000 iterations)")
print("- Optimal cutoff determined by Youden index (Sensitivity + Specificity - 1)")
print("- AUC interpretation: <0.6=poor, 0.6-0.7=fair, 0.7-0.8=good, 0.8-0.9=very good, >0.9=excellent")
print("=" * 100)
