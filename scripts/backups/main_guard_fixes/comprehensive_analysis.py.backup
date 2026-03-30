import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_curve, auc
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("全面细化统计分析")
print("=" * 80)

output_dir = '最终修改'

# 读取数据
df = pd.read_excel(f'{output_dir}/OCT数据_完整整合.xlsx')

# 分离两组
dep = df[df['分组'] == '抑郁症'].copy()
ctrl = df[df['分组'] == '健康对照'].copy()

print(f"数据概况: {len(df)} 行")
print(f"抑郁症组: {len(dep)} 行")
print(f"对照组: {len(ctrl)} 行")

# FDR校正函数
def fdr_correction(p_values, alpha=0.05):
    p_values = np.array(p_values)
    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]
    
    adjusted = np.zeros(n)
    for i in range(n-1, -1, -1):
        if i == n-1:
            adjusted[sorted_idx[i]] = sorted_p[i]
        else:
            adjusted[sorted_idx[i]] = min(sorted_p[i] * n / (i+1), adjusted[sorted_idx[i+1]])
    
    return adjusted

# 定义所有要分析的指标
# 1. 黄斑5层所有指标
layers = ['RNFL', 'Retina', 'GCL+', 'GCL++', 'Choroid']
regions = ['黄斑中心凹', '内环上方', '内环颞侧', '内环鼻侧', '内环下方', 
           '外环上方', '外环颞侧', '外环鼻侧', '外环下方']
metrics = ['平均厚度', '中心厚度', '总体积']

# 2. 视盘RNFL指标
rnfl_metrics = ['RNFL_Total', 'RNFL_上方', 'RNFL_颞侧', 'RNFL_鼻侧', 'RNFL_下方']

# 3. 视盘结构指标
disc_metrics = ['Disc Area', 'Cup Area', 'Rim Area', 'Cup Volume', 'Rim Volume', 
                'C/D Area Ratio', 'Linear C/D Ratio', 'Vertical C/D Ratio']

# ==================== 1. 黄斑所有层所有区域组间比较 ====================
print("\n" + "=" * 80)
print("【Table 2. 黄斑所有层所有区域组间比较】")
print("=" * 80)

macula_results = []

for layer in layers:
    for region in regions:
        col = f'{layer}_{region}'
        if col in df.columns:
            dep_vals = dep[col].dropna()
            ctrl_vals = ctrl[col].dropna()
            
            if len(dep_vals) > 10 and len(ctrl_vals) > 10:
                u_stat, p_val = stats.mannwhitneyu(dep_vals, ctrl_vals, alternative='two-sided')
                pooled_std = np.sqrt(((len(dep_vals)-1)*dep_vals.std()**2 + (len(ctrl_vals)-1)*ctrl_vals.std()**2) / (len(dep_vals)+len(ctrl_vals)-2))
                cohens_d = (dep_vals.mean() - ctrl_vals.mean()) / pooled_std
                
                macula_results.append({
                    '层': layer,
                    '区域': region,
                    '抑郁症组': f"{dep_vals.mean():.2f}±{dep_vals.std():.2f}",
                    '对照组': f"{ctrl_vals.mean():.2f}±{ctrl_vals.std():.2f}",
                    'P值': p_val,
                    'Cohen_d': cohens_d
                })

# FDR校正
if macula_results:
    p_values = [r['P值'] for r in macula_results]
    p_adj = fdr_correction(p_values)
    
    for i, r in enumerate(macula_results):
        r['P值_FDR'] = p_adj[i]
        if p_adj[i] < 0.001:
            r['显著性'] = '***'
        elif p_adj[i] < 0.01:
            r['显著性'] = '**'
        elif p_adj[i] < 0.05:
            r['显著性'] = '*'
        else:
            r['显著性'] = 'ns'

macula_df = pd.DataFrame(macula_results)
print(f"\n黄斑分析: {len(macula_df)} 个指标")
print(f"显著指标: {(macula_df['显著性'] != 'ns').sum()} 个")
print("\n前10个最显著的指标:")
print(macula_df.nsmallest(10, 'P值_FDR')[['层', '区域', 'P值_FDR', 'Cohen_d', '显著性']].to_string(index=False))

macula_df.to_excel(f'{output_dir}/Table2_黄斑所有层组间比较.xlsx', index=False)
print(f"\n✓ 已保存: {output_dir}/Table2_黄斑所有层组间比较.xlsx")

# ==================== 2. 黄斑平均厚度、总体积、中心厚度 ====================
print("\n" + "=" * 80)
print("【Table 3. 黄斑各层综合指标组间比较】")
print("=" * 80)

comprehensive_results = []

for layer in layers:
    for metric in metrics:
        col = f'{layer}_{metric}'
        if col in df.columns:
            dep_vals = dep[col].dropna()
            ctrl_vals = ctrl[col].dropna()
            
            if len(dep_vals) > 10 and len(ctrl_vals) > 10:
                u_stat, p_val = stats.mannwhitneyu(dep_vals, ctrl_vals, alternative='two-sided')
                pooled_std = np.sqrt(((len(dep_vals)-1)*dep_vals.std()**2 + (len(ctrl_vals)-1)*ctrl_vals.std()**2) / (len(dep_vals)+len(ctrl_vals)-2))
                cohens_d = (dep_vals.mean() - ctrl_vals.mean()) / pooled_std
                
                comprehensive_results.append({
                    '层': layer,
                    '指标': metric,
                    '抑郁症组': f"{dep_vals.mean():.2f}±{dep_vals.std():.2f}",
                    '对照组': f"{ctrl_vals.mean():.2f}±{ctrl_vals.std():.2f}",
                    'P值': p_val,
                    'Cohen_d': cohens_d
                })

if comprehensive_results:
    p_values = [r['P值'] for r in comprehensive_results]
    p_adj = fdr_correction(p_values)
    
    for i, r in enumerate(comprehensive_results):
        r['P值_FDR'] = p_adj[i]
        if p_adj[i] < 0.001:
            r['显著性'] = '***'
        elif p_adj[i] < 0.01:
            r['显著性'] = '**'
        elif p_adj[i] < 0.05:
            r['显著性'] = '*'
        else:
            r['显著性'] = 'ns'

comp_df = pd.DataFrame(comprehensive_results)
print(f"\n综合指标分析: {len(comp_df)} 个指标")
print(comp_df.to_string(index=False))

comp_df.to_excel(f'{output_dir}/Table3_黄斑综合指标组间比较.xlsx', index=False)
print(f"\n✓ 已保存: {output_dir}/Table3_黄斑综合指标组间比较.xlsx")

# ==================== 3. 视盘所有指标 ====================
print("\n" + "=" * 80)
print("【Table 4. 视盘所有指标组间比较】")
print("=" * 80)

disc_results = []

# RNFL指标
for col in rnfl_metrics:
    if col in df.columns:
        dep_vals = dep[col].dropna()
        ctrl_vals = ctrl[col].dropna()
        
        if len(dep_vals) > 10 and len(ctrl_vals) > 10:
            u_stat, p_val = stats.mannwhitneyu(dep_vals, ctrl_vals, alternative='two-sided')
            pooled_std = np.sqrt(((len(dep_vals)-1)*dep_vals.std()**2 + (len(ctrl_vals)-1)*ctrl_vals.std()**2) / (len(dep_vals)+len(ctrl_vals)-2))
            cohens_d = (dep_vals.mean() - ctrl_vals.mean()) / pooled_std
            
            disc_results.append({
                '类别': 'RNFL',
                '指标': col,
                '抑郁症组': f"{dep_vals.mean():.2f}±{dep_vals.std():.2f}",
                '对照组': f"{ctrl_vals.mean():.2f}±{ctrl_vals.std():.2f}",
                'P值': p_val,
                'Cohen_d': cohens_d
            })

# 视盘结构指标
for col in disc_metrics:
    if col in df.columns:
        dep_vals = dep[col].dropna()
        ctrl_vals = ctrl[col].dropna()
        
        if len(dep_vals) > 10 and len(ctrl_vals) > 10:
            u_stat, p_val = stats.mannwhitneyu(dep_vals, ctrl_vals, alternative='two-sided')
            pooled_std = np.sqrt(((len(dep_vals)-1)*dep_vals.std()**2 + (len(ctrl_vals)-1)*ctrl_vals.std()**2) / (len(dep_vals)+len(ctrl_vals)-2))
            cohens_d = (dep_vals.mean() - ctrl_vals.mean()) / pooled_std
            
            disc_results.append({
                '类别': '视盘结构',
                '指标': col,
                '抑郁症组': f"{dep_vals.mean():.2f}±{dep_vals.std():.2f}",
                '对照组': f"{ctrl_vals.mean():.2f}±{ctrl_vals.std():.2f}",
                'P值': p_val,
                'Cohen_d': cohens_d
            })

if disc_results:
    p_values = [r['P值'] for r in disc_results]
    p_adj = fdr_correction(p_values)
    
    for i, r in enumerate(disc_results):
        r['P值_FDR'] = p_adj[i]
        if p_adj[i] < 0.001:
            r['显著性'] = '***'
        elif p_adj[i] < 0.01:
            r['显著性'] = '**'
        elif p_adj[i] < 0.05:
            r['显著性'] = '*'
        else:
            r['显著性'] = 'ns'

disc_df = pd.DataFrame(disc_results)
print(f"\n视盘分析: {len(disc_df)} 个指标")
print(disc_df.to_string(index=False))

disc_df.to_excel(f'{output_dir}/Table4_视盘所有指标组间比较.xlsx', index=False)
print(f"\n✓ 已保存: {output_dir}/Table4_视盘所有指标组间比较.xlsx")

# ==================== 4. 全面相关分析 ====================
print("\n" + "=" * 80)
print("【Table 5. PHQ-9与所有OCT参数相关分析】")
print("=" * 80)

dep_with_phq9 = dep[dep['PHQ-9'].notna()].copy()
print(f"有PHQ-9数据的患者: {len(dep_with_phq9)} 人")

corr_results = []

# 黄斑所有层所有区域
for layer in layers:
    for region in regions:
        col = f'{layer}_{region}'
        if col in dep.columns:
            valid_data = dep_with_phq9[[col, 'PHQ-9']].dropna()
            if len(valid_data) > 20:
                rho, p_val = stats.spearmanr(valid_data[col], valid_data['PHQ-9'])
                corr_results.append({
                    '层': layer,
                    '区域': region,
                    'n': len(valid_data),
                    'Spearman_r': rho,
                    'P值': p_val
                })

# 黄斑综合指标
for layer in layers:
    for metric in metrics:
        col = f'{layer}_{metric}'
        if col in dep.columns:
            valid_data = dep_with_phq9[[col, 'PHQ-9']].dropna()
            if len(valid_data) > 20:
                rho, p_val = stats.spearmanr(valid_data[col], valid_data['PHQ-9'])
                corr_results.append({
                    '层': layer,
                    '指标': metric,
                    'n': len(valid_data),
                    'Spearman_r': rho,
                    'P值': p_val
                })

# 视盘指标
for col in rnfl_metrics + disc_metrics:
    if col in dep.columns:
        valid_data = dep_with_phq9[[col, 'PHQ-9']].dropna()
        if len(valid_data) > 20:
            rho, p_val = stats.spearmanr(valid_data[col], valid_data['PHQ-9'])
            corr_results.append({
                '类别': '视盘',
                '指标': col,
                'n': len(valid_data),
                'Spearman_r': rho,
                'P值': p_val
            })

corr_df = pd.DataFrame(corr_results)
print(f"\n相关分析: {len(corr_df)} 个指标")
print(f"显著相关: {(corr_df['P值'] < 0.05).sum()} 个")

if len(corr_df) > 0:
    print("\n最显著的相关（P<0.05）:")
    sig_corr = corr_df[corr_df['P值'] < 0.05]
    if len(sig_corr) > 0:
        print(sig_corr.to_string(index=False))
    else:
        print("无显著相关")

corr_df.to_excel(f'{output_dir}/Table5_全面相关分析.xlsx', index=False)
print(f"\n✓ 已保存: {output_dir}/Table5_全面相关分析.xlsx")

# ==================== 5. 全面ROC分析 ====================
print("\n" + "=" * 80)
print("【Table 6. 全面ROC分析】")
print("=" * 80)

# 准备数据（双眼取平均）
df['分组_编码'] = (df['分组'] == '抑郁症').astype(int)

# 所有要分析的列
all_cols = []
for layer in layers:
    for region in regions:
        all_cols.append(f'{layer}_{region}')
    for metric in metrics:
        all_cols.append(f'{layer}_{metric}')
all_cols.extend(rnfl_metrics)
all_cols.extend(disc_metrics)

# 按Patient_ID分组
df_patient = df.groupby(['Patient_ID', '分组_编码']).agg({col: 'mean' for col in all_cols if col in df.columns}).reset_index()

print(f"患者数: {len(df_patient)}")
print(f"抑郁症: {df_patient['分组_编码'].sum()} 人")
print(f"对照: {len(df_patient) - df_patient['分组_编码'].sum()} 人")

roc_results = []

for col in all_cols:
    if col in df_patient.columns:
        valid_data = df_patient[[col, '分组_编码']].dropna()
        
        if len(valid_data) > 30:
            y_true = valid_data['分组_编码'].values
            y_scores = valid_data[col].values
            
            # 计算ROC
            fpr, tpr, thresholds = roc_curve(y_true, -y_scores)
            roc_auc = auc(fpr, tpr)
            
            # 最佳截断点
            youden = tpr - fpr
            best_idx = np.argmax(youden)
            best_threshold = -thresholds[best_idx]
            sensitivity = tpr[best_idx]
            specificity = 1 - fpr[best_idx]
            
            roc_results.append({
                '指标': col,
                'AUC': roc_auc,
                '敏感度': sensitivity,
                '特异度': specificity,
                '最佳截断值': best_threshold
            })

roc_df = pd.DataFrame(roc_results)
roc_df = roc_df.sort_values('AUC', ascending=False)

print(f"\nROC分析: {len(roc_df)} 个指标")
print(f"AUC>0.6: {(roc_df['AUC'] > 0.6).sum()} 个")
print(f"AUC>0.7: {(roc_df['AUC'] > 0.7).sum()} 个")

print("\nTop 10 AUC指标:")
print(roc_df.head(10).to_string(index=False))

roc_df.to_excel(f'{output_dir}/Table6_全面ROC分析.xlsx', index=False)
print(f"\n✓ 已保存: {output_dir}/Table6_全面ROC分析.xlsx")

# ==================== 6. 生成总结报告 ====================
print("\n" + "=" * 80)
print("【全面分析总结】")
print("=" * 80)

print(f"""
1. 黄斑分析:
   - 总指标数: {len(macula_df)}
   - 显著指标: {(macula_df['显著性'] != 'ns').sum()} ({(macula_df['显著性'] != 'ns').sum()/len(macula_df)*100:.1f}%)
   - 最显著: {macula_df.loc[macula_df['P值_FDR'].idxmin(), '层']}_{macula_df.loc[macula_df['P值_FDR'].idxmin(), '区域']} (P={macula_df['P值_FDR'].min():.2e})

2. 黄斑综合指标:
   - 总指标数: {len(comp_df)}
   - 显著指标: {(comp_df['显著性'] != 'ns').sum()}

3. 视盘分析:
   - 总指标数: {len(disc_df)}
   - 显著指标: {(disc_df['显著性'] != 'ns').sum()}

4. 相关分析:
   - 总指标数: {len(corr_df)}
   - 显著相关: {(corr_df['P值'] < 0.05).sum()}

5. ROC分析:
   - 总指标数: {len(roc_df)}
   - AUC>0.6: {(roc_df['AUC'] > 0.6).sum()}
   - AUC>0.7: {(roc_df['AUC'] > 0.7).sum()}
   - 最佳AUC: {roc_df['AUC'].max():.3f} ({roc_df.loc[roc_df['AUC'].idxmax(), '指标']})
""")

print("\n" + "=" * 80)
print("全面细化分析完成!")
print("=" * 80)
