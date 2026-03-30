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
def bootstrap_auc_ci(y_true, y_score, n_bootstraps=1000, ci=0.95, random_state=42):
    """
    使用Bootstrap计算AUC的置信区间
    """
    rng = np.random.RandomState(random_state)
    bootstrapped_aucs = []
    
    for i in range(n_bootstraps):
        indices = rng.randint(0, len(y_true), len(y_true))
        if len(np.unique(y_true[indices])) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true[indices], y_score[indices])
        auc_score = auc(fpr, tpr)
        if auc_score < 0.5:
            auc_score = 1 - auc_score
        bootstrapped_aucs.append(auc_score)
    
    alpha = (1 - ci) / 2
    ci_lower = np.percentile(bootstrapped_aucs, alpha * 100)
    ci_upper = np.percentile(bootstrapped_aucs, (1 - alpha) * 100)
    
    return ci_lower, ci_upper

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
    auc_original = auc(fpr, tpr)
    
    # 确定方向（MDD组是否更低）
    mdd_mean = np.mean(y_score[y_true == 1])
    control_mean = np.mean(y_score[y_true == 0])
    mdd_lower = mdd_mean < control_mean
    
    # 计算AUC（确保>0.5）
    if auc_original < 0.5:
        auc_value = 1 - auc_original
        # 反转FPR和TPR用于找最佳截断点
        fpr, tpr = 1 - fpr, 1 - tpr
    else:
        auc_value = auc_original
    
    # 计算Bootstrap 95% CI
    auc_ci_lower, auc_ci_upper = bootstrap_auc_ci(y_true, y_score, n_bootstraps=2000)
    
    # 找到最佳截断点 (Youden指数)
    youden_index = tpr - fpr
    optimal_idx = np.argmax(youden_index)
    optimal_threshold = thresholds[optimal_idx]
    
    # 根据方向确定分类规则
    if mdd_lower:
        # MDD组数值更低，小于截断点为阳性
        y_pred = (y_score <= optimal_threshold).astype(int)
    else:
        # MDD组数值更高，大于截断点为阳性
        y_pred = (y_score >= optimal_threshold).astype(int)
    
    # 计算敏感度和特异度
    tp = np.sum((y_true == 1) & (y_pred == 1))
    tn = np.sum((y_true == 0) & (y_pred == 0))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    fn = np.sum((y_true == 1) & (y_pred == 0))
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    results.append({
        'Indicator': name,
        'AUC': round(auc_value, 3),
        'AUC_CI_Lower': round(auc_ci_lower, 3),
        'AUC_CI_Upper': round(auc_ci_upper, 3),
        'Sensitivity': round(sensitivity, 3),
        'Specificity': round(specificity, 3),
        'Optimal_Cutoff': round(optimal_threshold, 2),
        'N': len(valid_data),
        'MDD_Mean': round(mdd_mean, 2),
        'Control_Mean': round(control_mean, 2)
    })
    
    print(f"{name:<30} {auc_value:.3f}      [{auc_ci_lower:.3f}–{auc_ci_upper:.3f}]{'' :<6} {sensitivity:.3f}      {specificity:.3f}      {optimal_threshold:.2f}")

print("=" * 100)

# 创建结果DataFrame
results_df = pd.DataFrame(results)

# 保存结果
output_path = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/ROC_Analysis_Final_with_95CI.xlsx'
results_df.to_excel(output_path, index=False)
print(f"\nResults saved to: {output_path}")

# 同时保存到workspace
workspace_path = '/root/.openclaw/workspace/revised_paper/ROC_Analysis_Final_with_95CI.xlsx'
results_df.to_excel(workspace_path, index=False)
print(f"Results saved to: {workspace_path}")

# 打印用于论文的格式
print("\n" + "=" * 100)
print("Results Formatted for Manuscript (Table 5 Format)")
print("=" * 100)

for _, row in results_df.iterrows():
    print(f"\n{row['Indicator']}:")
    print(f"  AUC = {row['AUC']:.3f} (95% CI: {row['AUC_CI_Lower']:.3f}–{row['AUC_CI_Upper']:.3f})")
    print(f"  Sensitivity = {row['Sensitivity']:.1%} ({row['Sensitivity']*100:.1f}%)")
    print(f"  Specificity = {row['Specificity']:.1%} ({row['Specificity']*100:.1f}%)")
    print(f"  Optimal cutoff = {row['Optimal_Cutoff']:.2f}")
    print(f"  MDD mean = {row['MDD_Mean']:.2f}, Control mean = {row['Control_Mean']:.2f}")

# 创建汇总表格
print("\n" + "=" * 100)
print("Summary Table for Manuscript (Markdown Format)")
print("=" * 100)
print("\n| OCT Parameter | AUC (95% CI) | Sensitivity | Specificity | Optimal Cutoff |")
print("|---------------|--------------|-------------|-------------|----------------|")
for _, row in results_df.iterrows():
    print(f"| {row['Indicator']} | {row['AUC']:.3f} ({row['AUC_CI_Lower']:.3f}–{row['AUC_CI_Upper']:.3f}) | {row['Sensitivity']:.1%} | {row['Specificity']:.1%} | {row['Optimal_Cutoff']:.2f} |")

print("\n" + "=" * 100)
print("Statistical Notes:")
print("- 95% CI calculated using Bootstrap method (2000 iterations)")
print("- Optimal cutoff determined by Youden index (Sensitivity + Specificity - 1)")
print("- AUC interpretation: <0.6=poor, 0.6-0.7=fair, 0.7-0.8=good, 0.8-0.9=very good, >0.9=excellent")
print("- All reported OCT parameters showed higher values in controls vs MDD (except C/D ratio)")
print("=" * 100)
