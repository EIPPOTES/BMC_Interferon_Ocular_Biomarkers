#!/usr/bin/env python3

def main():
    """主函数，包装原有执行代码"""
    """
    统一分析脚本：使用同一数据集重新运行所有分析
    基于data.xlsx (499眼版本)，确保所有分析使用相同样本
    """

    import pandas as pd
    import numpy as np
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from scipy import stats
    from sklearn.metrics import roc_curve, auc
    import warnings
    warnings.filterwarnings('ignore')

    print("=" * 80)
    print("统一分析：使用同一数据集重新运行所有分析")
    print("数据源: data.xlsx (499眼版本)")
    print("=" * 80)

    # 读取数据
    df = pd.read_excel('data.xlsx')
    print(f"原始数据: {len(df)} 眼")
    print(f"抑郁症: {df[df['分组'] == '抑郁症'].shape[0]} 眼")
    print(f"健康对照: {df[df['分组'] == '健康对照'].shape[0]} 眼")

    # 1. 数据准备和清洗
    # 创建抑郁状态编码
    df['抑郁状态'] = (df['分组'] == '抑郁症').astype(int)

    # 创建PHQ-9严重程度分组
    def categorize_phq9(score):
        if pd.isna(score):
            return np.nan
        elif score < 5:
            return '无症状'
        elif score < 10:
            return '轻度'
        elif score < 15:
            return '中度'
        else:
            return '重度'

    df['PHQ9_严重程度'] = df['PHQ-9'].apply(categorize_phq9)

    # 创建性别编码（女=0，男=1）
    df['性别_编码'] = df['性别'].map({'女': 0, '男': 1})

    # 选择分析样本：有年龄和性别数据的样本
    df_analysis = df.dropna(subset=['年龄', '性别'])
    print(f"\n分析样本（有年龄和性别数据）: {len(df_analysis)} 眼")
    print(f"抑郁症: {df_analysis[df_analysis['分组'] == '抑郁症'].shape[0]} 眼")
    print(f"健康对照: {df_analysis[df_analysis['分组'] == '健康对照'].shape[0]} 眼")
    print(f"年龄范围: {df_analysis['年龄'].min():.1f} - {df_analysis['年龄'].max():.1f} 岁")

    # 定义关键OCT指标（基于之前分析选择最显著的）
    key_indicators = [
        'Retina_外环颞侧',
        'Retina_外环上方',
        'Retina_平均厚度',
        'Retina_总体积',
        'Retina_内环颞侧',
        'RNFL_平均厚度',
        'RNFL_外环上方',
        'GCL+_平均厚度',
        'Cup_to_Disc_Ratio',
        'Rim_Volume'
    ]

    # 检查指标是否存在
    available_indicators = [col for col in key_indicators if col in df_analysis.columns]
    print(f"\n可用OCT指标: {len(available_indicators)}/{len(key_indicators)}")
    print("主要分析指标:", available_indicators[:5])

    # 2. 多变量回归分析（线性模型，控制年龄、性别）
    print("\n" + "=" * 80)
    print("2. 多变量回归分析（控制年龄、性别）")
    print("=" * 80)

    linear_results = []

    for indicator in available_indicators:
        # 删除当前指标的缺失值
        df_temp = df_analysis.dropna(subset=[indicator, '年龄', '性别_编码', '抑郁状态'])

        if len(df_temp) < 50:  # 样本太小
            continue

        # 标准化变量（年龄标准化以便解释）
        df_temp['年龄_标准化'] = (df_temp['年龄'] - df_temp['年龄'].mean()) / df_temp['年龄'].std()

        # 准备回归变量
        X = pd.DataFrame({
            '抑郁状态': df_temp['抑郁状态'],
            '年龄': df_temp['年龄_标准化'],
            '性别': df_temp['性别_编码'],
            '截距': 1
        })
        y = df_temp[indicator]

        # 线性回归
        model = sm.OLS(y, X).fit()

        # 提取结果
        beta = model.params['抑郁状态']
        p_value = model.pvalues['抑郁状态']
        ci_low, ci_high = model.conf_int().loc['抑郁状态']

        # 标准化系数
        y_std = y.std()
        beta_std = beta / y_std

        linear_results.append({
            '指标': indicator,
            '系数β': round(beta, 3),
            '标准化β': round(beta_std, 3),
            'P值': p_value,
            '95%CI下限': round(ci_low, 3),
            '95%CI上限': round(ci_high, 3),
            '样本量': len(df_temp),
            '模型R²': round(model.rsquared, 3)
        })

    # 创建线性回归结果DataFrame
    linear_df = pd.DataFrame(linear_results)
    linear_df = linear_df.sort_values('P值')
    print(f"\n分析完成: {len(linear_df)} 个指标")
    print(linear_df[['指标', '系数β', '标准化β', 'P值', '样本量']].head(10).to_string(index=False))

    # 3. 混合效应模型（考虑双眼相关性）
    print("\n" + "=" * 80)
    print("3. 混合效应模型（考虑双眼相关性）")
    print("=" * 80)

    mixed_results = []

    # 选择几个关键指标进行混合效应模型（计算量较大）
    key_mixed_indicators = available_indicators[:5]  # 前5个指标

    for indicator in key_mixed_indicators:
        # 删除缺失值
        df_temp = df_analysis.dropna(subset=[indicator, '年龄', '性别_编码', '抑郁状态', 'Patient_ID'])

        if len(df_temp) < 50:
            continue

        try:
            # 标准化年龄
            df_temp['年龄_标准化'] = (df_temp['年龄'] - df_temp['年龄'].mean()) / df_temp['年龄'].std()

            # 混合效应模型
            formula = f"{indicator} ~ 抑郁状态 + 年龄_标准化 + 性别_编码"
            model = smf.mixedlm(formula, df_temp, groups=df_temp["Patient_ID"]).fit()

            # 提取结果
            beta = model.params['抑郁状态']
            p_value = model.pvalues['抑郁状态']

            # 计算置信区间（近似）
            se = model.bse['抑郁状态']
            ci_low = beta - 1.96 * se
            ci_high = beta + 1.96 * se

            mixed_results.append({
                '指标': indicator,
                '系数β': round(beta, 3),
                'P值': p_value,
                '95%CI下限': round(ci_low, 3),
                '95%CI上限': round(ci_high, 3),
                '样本量': len(df_temp),
                '患者数': df_temp['Patient_ID'].nunique()
            })

        except Exception as e:
            print(f"  指标 {indicator} 混合效应模型失败: {str(e)[:50]}")

    mixed_df = pd.DataFrame(mixed_results)
    print(f"\n混合效应模型完成: {len(mixed_df)} 个指标")
    if len(mixed_df) > 0:
        print(mixed_df[['指标', '系数β', 'P值', '样本量', '患者数']].to_string(index=False))

    # 4. 亚组分析（性别、年龄、PHQ-9）
    print("\n" + "=" * 80)
    print("4. 亚组分析")
    print("=" * 80)

    # 4.1 性别亚组
    print("\n4.1 性别亚组分析")
    gender_subgroups = []

    # 年龄中位数用于年龄分层
    age_median = df_analysis['年龄'].median()
    print(f"年龄中位数: {age_median:.1f} 岁")

    for gender in ['女', '男']:
        gender_df = df_analysis[df_analysis['性别'] == gender]

        if len(gender_df) < 30:
            continue

        for indicator in available_indicators[:5]:  # 分析前5个指标
            # 删除缺失值
            df_temp = gender_df.dropna(subset=[indicator, '年龄', '抑郁状态'])

            if len(df_temp) < 20:
                continue

            # 标准化年龄
            df_temp['年龄_标准化'] = (df_temp['年龄'] - df_temp['年龄'].mean()) / df_temp['年龄'].std()

            # 回归分析
            X = pd.DataFrame({
                '抑郁状态': df_temp['抑郁状态'],
                '年龄': df_temp['年龄_标准化'],
                '截距': 1
            })
            y = df_temp[indicator]

            try:
                model = sm.OLS(y, X).fit()
                beta = model.params['抑郁状态']
                p_value = model.pvalues['抑郁状态']

                gender_subgroups.append({
                    '亚组': gender,
                    '指标': indicator,
                    '系数β': round(beta, 3),
                    'P值': p_value,
                    '样本量': len(df_temp),
                    '标准化β': round(beta / y.std(), 3)
                })
            except:
                pass

    gender_subgroup_df = pd.DataFrame(gender_subgroups)
    print(f"性别亚组分析完成: {len(gender_subgroup_df)} 个结果")

    # 4.2 年龄亚组（基于中位数）
    print("\n4.2 年龄亚组分析")
    age_subgroups = []

    for age_group, condition in [('年轻组 (<中位数)', f'年龄 < {age_median}'), 
                                 ('年长组 (≥中位数)', f'年龄 >= {age_median}')]:
        age_df = df_analysis.query(condition)

        if len(age_df) < 30:
            continue

        for indicator in available_indicators[:5]:
            df_temp = age_df.dropna(subset=[indicator, '性别_编码', '抑郁状态'])

            if len(df_temp) < 20:
                continue

            # 标准化变量
            df_temp['性别_编码'] = df_temp['性别_编码']

            # 回归分析
            X = pd.DataFrame({
                '抑郁状态': df_temp['抑郁状态'],
                '性别': df_temp['性别_编码'],
                '截距': 1
            })
            y = df_temp[indicator]

            try:
                model = sm.OLS(y, X).fit()
                beta = model.params['抑郁状态']
                p_value = model.pvalues['抑郁状态']

                age_subgroups.append({
                    '亚组': age_group,
                    '指标': indicator,
                    '系数β': round(beta, 3),
                    'P值': p_value,
                    '样本量': len(df_temp)
                })
            except:
                pass

    age_subgroup_df = pd.DataFrame(age_subgroups)
    print(f"年龄亚组分析完成: {len(age_subgroup_df)} 个结果")

    # 5. ROC分析
    print("\n" + "=" * 80)
    print("5. ROC分析（诊断性能评估）")
    print("=" * 80)

    # 按患者取平均（避免双眼数据的相关性）
    patient_df = df_analysis.groupby(['Patient_ID', '抑郁状态']).agg({
        '年龄': 'mean',
        '性别': lambda x: x.mode()[0] if not x.mode().empty else np.nan
    }).reset_index()

    # 为每个指标添加平均值
    for indicator in available_indicators[:7]:  # 分析前7个指标
        indicator_means = df_analysis.groupby('Patient_ID')[indicator].mean()
        patient_df = patient_df.merge(indicator_means.rename(indicator), left_on='Patient_ID', right_index=True, how='left')

    # 删除缺失值
    patient_df_clean = patient_df.dropna(subset=available_indicators[:7] + ['抑郁状态'])

    print(f"患者样本: {len(patient_df_clean)} 人")
    print(f"抑郁症患者: {patient_df_clean['抑郁状态'].sum()} 人")

    roc_results = []

    for indicator in available_indicators[:7]:
        if indicator not in patient_df_clean.columns:
            continue

        # 删除该指标的缺失值
        roc_data = patient_df_clean.dropna(subset=[indicator, '抑郁状态'])

        if len(roc_data) < 50:
            continue

        y_true = roc_data['抑郁状态']
        y_score = roc_data[indicator]

        # 计算ROC曲线
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)

        # 寻找最佳截断点（Youden's J statistic）
        youden_j = tpr - fpr
        optimal_idx = np.argmax(youden_j)
        optimal_threshold = thresholds[optimal_idx]

        # 预测
        y_pred = (y_score >= optimal_threshold).astype(int)

        # 计算敏感度和特异度
        from sklearn.metrics import confusion_matrix
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        roc_results.append({
            '指标': indicator,
            'AUC': round(roc_auc, 3),
            '敏感度': round(sensitivity, 3),
            '特异度': round(specificity, 3),
            '最佳截断值': round(optimal_threshold, 2),
            '样本量': len(roc_data)
        })

    roc_df = pd.DataFrame(roc_results)
    roc_df = roc_df.sort_values('AUC', ascending=False)
    print(f"\nROC分析完成: {len(roc_df)} 个指标")
    print(roc_df[['指标', 'AUC', '敏感度', '特异度', '样本量']].to_string(index=False))

    # 6. 相关性分析（OCT vs PHQ-9）
    print("\n" + "=" * 80)
    print("6. 相关性分析（OCT指标 vs PHQ-9评分）")
    print("=" * 80)

    # 选择有PHQ-9数据的抑郁患者
    phq9_df = df_analysis[(df_analysis['分组'] == '抑郁症') & (df_analysis['PHQ-9'].notna())]

    print(f"有PHQ-9数据的抑郁患者: {len(phq9_df)} 眼")
    if len(phq9_df) > 0:
        print(f"PHQ-9评分范围: {phq9_df['PHQ-9'].min():.1f} - {phq9_df['PHQ-9'].max():.1f}")
        print(f"PHQ-9均值: {phq9_df['PHQ-9'].mean():.2f} ± {phq9_df['PHQ-9'].std():.2f}")

    correlation_results = []

    for indicator in available_indicators[:7]:
        corr_data = phq9_df.dropna(subset=[indicator, 'PHQ-9'])

        if len(corr_data) < 30:
            continue

        # Spearman相关性
        rho, p_value = stats.spearmanr(corr_data[indicator], corr_data['PHQ-9'])

        # 效应大小分类
        effect_size = "可忽略"
        if abs(rho) >= 0.3:
            effect_size = "中"
        elif abs(rho) >= 0.1:
            effect_size = "小"

        correlation_results.append({
            '指标': indicator,
            'Spearman_rho': round(rho, 3),
            'Spearman_P值': p_value,
            '效应大小': effect_size,
            '样本量': len(corr_data)
        })

    corr_df = pd.DataFrame(correlation_results)
    corr_df = corr_df.sort_values('Spearman_P值')
    print(f"\n相关性分析完成: {len(corr_df)} 个指标")
    print(corr_df[['指标', 'Spearman_rho', 'Spearman_P值', '效应大小', '样本量']].to_string(index=False))

    # 7. 保存所有结果
    print("\n" + "=" * 80)
    print("7. 保存分析结果")
    print("=" * 80)

    # 创建输出目录
    import os
    output_dir = 'unified_analysis_results_20260315'
    os.makedirs(output_dir, exist_ok=True)

    # 保存结果到Excel
    with pd.ExcelWriter(f'{output_dir}/统一分析结果汇总.xlsx') as writer:
        linear_df.to_excel(writer, sheet_name='多变量回归_线性模型', index=False)
        if len(mixed_df) > 0:
            mixed_df.to_excel(writer, sheet_name='多变量回归_混合效应模型', index=False)

        # 合并亚组分析结果
        all_subgroups = pd.concat([
            gender_subgroup_df.assign(分析类型='性别亚组'),
            age_subgroup_df.assign(分析类型='年龄亚组')
        ], ignore_index=True)
        all_subgroups.to_excel(writer, sheet_name='亚组分析', index=False)

        roc_df.to_excel(writer, sheet_name='ROC分析', index=False)
        corr_df.to_excel(writer, sheet_name='相关性分析', index=False)

    print(f"所有结果已保存到: {output_dir}/统一分析结果汇总.xlsx")

    # 8. 生成总结报告
    print("\n" + "=" * 80)
    print("8. 分析总结报告")
    print("=" * 80)

    print(f"\n📊 统一分析完成 - 使用样本: {len(df_analysis)} 眼")
    print(f"   抑郁症: {df_analysis[df_analysis['分组'] == '抑郁症'].shape[0]} 眼")
    print(f"   健康对照: {df_analysis[df_analysis['分组'] == '健康对照'].shape[0]} 眼")

    print(f"\n🔍 主要发现:")
    print(f"   1. 最显著的OCT指标: {linear_df.iloc[0]['指标']} (β={linear_df.iloc[0]['系数β']}, P={linear_df.iloc[0]['P值']:.2e})")
    print(f"   2. 最佳诊断指标: {roc_df.iloc[0]['指标']} (AUC={roc_df.iloc[0]['AUC']})")
    print(f"   3. 与PHQ-9最相关的指标: {corr_df.iloc[0]['指标']} (ρ={corr_df.iloc[0]['Spearman_rho']})")

    print(f"\n📈 关键统计量:")
    print(f"   - 多变量回归分析: {len(linear_df)} 个指标")
    print(f"   - 混合效应模型: {len(mixed_df)} 个指标")
    print(f"   - ROC分析: {len(roc_df)} 个指标 (最高AUC: {roc_df.iloc[0]['AUC']})")
    print(f"   - 相关性分析: {len(corr_df)} 个指标")

    print(f"\n💾 生成文件:")
    print(f"   - {output_dir}/统一分析结果汇总.xlsx (所有结果)")
    print(f"   - 各工作表: 多变量回归、亚组分析、ROC分析、相关性分析")

    print(f"\n🎯 后续建议:")
    print(f"   1. 使用统一结果更新论文数据")
    print(f"   2. 检查与之前版本的差异")
    print(f"   3. 确保方法学一致性")

    print("\n" + "=" * 80)
    print("统一分析完成!")
    print("=" * 80)


if __name__ == "__main__":
    main()