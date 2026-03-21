#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
增强Figure 3 - 添加更多ROC参数
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

# 色盲友好配色
cb_palette = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

output_dir = '/mnt/c/Users/CUI/Desktop/投稿版/Figures_修订版'

def load_data():
    """加载数据"""
    data_file = '/mnt/c/Users/CUI/Desktop/最终版/04_原始数据/02_OCT数据_完整整合.xlsx'
    df = pd.read_excel(data_file)
    df['分组_编码'] = (df['分组'] == '抑郁症').astype(int)
    df_patient = df.groupby('Patient_ID').agg({
        'Retina_外环颞侧': 'mean',
        'Retina_平均厚度': 'mean',
        'Retina_内环颞侧': 'mean',
        'Retina_外环上方': 'mean',
        'RNFL_Total': 'mean',
        '分组_编码': 'first'
    }).dropna()
    return df_patient

def enhance_figure3():
    """增强Figure 3 - 添加更多ROC参数"""
    print("增强 Figure 3: 添加更多ROC参数...")
    
    df_patient = load_data()
    y_true = df_patient['分组_编码'].values
    
    # 选择6个ROC参数（从2个增加到6个）
    roc_params = [
        ('Retina_外环颞侧', 'Outer Temporal'),
        ('Retina_平均厚度', 'Mean Macular'),
        ('Retina_内环颞侧', 'Inner Temporal'),
        ('Retina_外环上方', 'Outer Superior'),
        ('RNFL_Total', 'Total RNFL'),
    ]
    
    # 创建大图
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], hspace=0.3, wspace=0.3)
    
    # 上图: ROC曲线
    ax1 = fig.add_subplot(gs[0, :])
    
    results_table = []
    
    for (col, name), color in zip(roc_params, cb_palette):
        y_scores = -df_patient[col].values
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        # Bootstrap CI
        n_bootstraps = 1000
        rng = np.random.RandomState(42)
        bootstrapped_scores = []
        
        for i in range(n_bootstraps):
            indices = rng.randint(0, len(y_true), len(y_true))
            if len(np.unique(y_true[indices])) < 2:
                continue
            fpr_boot, tpr_boot, _ = roc_curve(y_true[indices], y_scores[indices])
            bootstrapped_scores.append(auc(fpr_boot, tpr_boot))
        
        ci_lower = np.percentile(bootstrapped_scores, 2.5)
        ci_upper = np.percentile(bootstrapped_scores, 97.5)
        
        # 找到最佳阈值
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        sensitivity = tpr[optimal_idx]
        specificity = 1 - fpr[optimal_idx]
        
        results_table.append({
            'Parameter': name,
            'AUC': f'{roc_auc:.3f}',
            '95% CI': f'({ci_lower:.3f}-{ci_upper:.3f})',
            'Sensitivity': f'{sensitivity:.3f}',
            'Specificity': f'{specificity:.3f}',
            'Threshold': f'{optimal_threshold:.1f}'
        })
        
        # 绘制ROC曲线
        ax1.plot(fpr, tpr, color=color, lw=2.5, 
                label=f'{name} (AUC = {roc_auc:.3f})')
    
    # 对角线
    ax1.plot([0, 1], [0, 1], color='gray', lw=1.5, linestyle='--', 
            label='Random (AUC = 0.500)')
    
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=12)
    ax1.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax1.set_title('A. ROC Curves for OCT Parameters (N=5 Parameters)', 
                 fontsize=13, fontweight='bold', loc='left')
    ax1.legend(loc='lower right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 左下图: AUC对比条形图
    ax2 = fig.add_subplot(gs[1, 0])
    
    auc_values = [float(r['AUC']) for r in results_table]
    param_names = [r['Parameter'] for r in results_table]
    
    bars = ax2.barh(param_names, auc_values, color=cb_palette[:len(auc_values)], 
                   alpha=0.7, edgecolor='black')
    ax2.axvline(x=0.5, color='gray', linestyle='--', label='Random (0.5)')
    ax2.axvline(x=0.7, color='red', linestyle='--', label='Good (0.7)')
    ax2.set_xlim(0.4, 0.8)
    ax2.set_xlabel('AUC', fontsize=11)
    ax2.set_title('B. AUC Comparison', fontsize=12, fontweight='bold', loc='left')
    
    # 添加数值标签
    for i, (bar, auc_val) in enumerate(zip(bars, auc_values)):
        ax2.text(auc_val + 0.01, i, f'{auc_val:.3f}', 
                va='center', fontsize=9, fontweight='bold')
    
    ax2.legend(fontsize=9)
    
    # 右下图: 临床性能表格
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    table_data = [[r['Parameter'], r['AUC'], r['Sensitivity'], r['Specificity']] 
                  for r in results_table]
    
    table = ax3.table(cellText=table_data,
                     colLabels=['Parameter', 'AUC', 'Sensitivity', 'Specificity'],
                     cellLoc='center',
                     loc='center',
                     bbox=[0, 0.2, 1, 0.7])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)
    
    for i in range(4):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax3.set_title('C. Clinical Performance at Optimal Threshold', 
                 fontsize=12, fontweight='bold', pad=20)
    
    # 总标题
    fig.suptitle('Figure 3. Diagnostic Performance of OCT Parameters for MDD Detection\n'
                'Comprehensive ROC Analysis with 5 Key Parameters',
                fontsize=14, fontweight='bold', y=0.98)
    
    # 图注
    fig.text(0.5, 0.02, 
             'Note: AUC = area under the curve; CI = confidence interval from 1000 bootstrap samples. '
             'Optimal threshold determined by Youden index. '
             'Red dashed line (AUC=0.7) indicates good diagnostic performance.',
             ha='center', fontsize=9, style='italic',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    plt.savefig(f'{output_dir}/Figure3_ROC_Curves_增强版.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✅ 已保存: Figure3_ROC_Curves_增强版.png")
    print(f"     包含5个参数:")
    for r in results_table:
        print(f"       - {r['Parameter']}: AUC={r['AUC']}")

def main():
    print("="*70)
    print("增强 Figure 3 - 添加更多ROC数据")
    print("="*70)
    
    enhance_figure3()
    
    print("\n" + "="*70)
    print("✅ Figure 3 增强完成!")
    print("="*70)
    print("""
增强内容:
  - 从2个参数增加到5个参数
  - 新增: Inner Temporal, Outer Superior, Total RNFL
  - 添加AUC对比条形图
  - 添加临床性能表格
  - 三面板布局，信息更丰富
    """)

if __name__ == "__main__":
    main()