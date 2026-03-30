#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *

def main():
    """主函数，包装原有执行代码"""
    """
    基于463眼版本更新ROC表格数据
    重新计算主要OCT指标的ROC性能
    """

    import pandas as pd
    import numpy as np
    from sklearn.metrics import roc_curve, auc, roc_auc_score
    from sklearn.utils import resample
    import os
    from datetime import datetime

    print("=" * 80)
    print("基于463眼版本更新ROC表格数据")
    print(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    # 路径设置
    data_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改"
    data_path = os.path.join(data_dir, "04_Data", "data_499eyes_20260315.xlsx")
    tables_dir = os.path.join(data_dir, "03_Tables")
    output_dir = os.path.join(data_dir, "01_Manuscript")

    # 1. 读取数据
    print(f"\n📊 读取数据...")
    df = pd.read_excel(data_path)

    # 筛选463眼有完整年龄性别数据的样本
    df_463 = df.dropna(subset=['年龄', '性别']).copy()
    print(f"原始数据: {df.shape[0]} 眼")
    print(f"463眼样本: {df_463.shape[0]} 眼")

    # 创建分组变量：抑郁症=1，健康对照=0
    df_463['group_numeric'] = df_463['分组'].apply(lambda x: 1 if x == '抑郁症' else 0)

    # 2. 识别关键OCT指标
    print(f"\n🔍 识别关键OCT指标...")

    # 基于论文中提到的关键指标
    key_indicators = [
        'Retina_外环颞侧',    # Outer temporal thickness
        'Retina_平均厚度',     # Mean macular thickness  
        'Retina_内环颞侧',    # Inner temporal thickness
        'Retina_总体积',      # Total macular volume
        'Retina_外环下方',    # Outer inferior thickness
        'RNFL_平均厚度',      # Peripapillary RNFL (average)
        'GCL+_平均厚度',      # GCL+ thickness
        'GCL++_平均厚度',     # GCL++ thickness
        'Choroid_平均厚度',   # Mean choroidal thickness
        'Rim Volume',         # Rim volume
        'C/D Area Ratio'      # Cup-to-disc ratio
    ]

    # 检查哪些指标在数据中
    available_indicators = []
    for indicator in key_indicators:
        if indicator in df_463.columns:
            available_indicators.append(indicator)
        else:
            # 尝试查找相似的列名
            matching_cols = [col for col in df_463.columns if indicator.lower() in col.lower()]
            if matching_cols:
                available_indicators.append(matching_cols[0])
                print(f"  使用替代列名: {indicator} -> {matching_cols[0]}")
            else:
                print(f"  警告: 未找到指标 {indicator}")

    print(f"找到 {len(available_indicators)}/{len(key_indicators)} 个关键指标")

    # 3. 计算每个指标的ROC性能
    print(f"\n📈 计算ROC性能...")

    roc_results = []

    for indicator in available_indicators:
        # 筛选有该指标和分组数据的样本
        temp_df = df_463.dropna(subset=[indicator, 'group_numeric'])

        if len(temp_df) < 50:  # 样本量太小
            continue

        # 获取数据和标签
        X = temp_df[indicator].values
        y = temp_df['group_numeric'].values

        # 计算AUC
        try:
            auc_score = roc_auc_score(y, X)

            # 计算敏感度、特异度、最佳截断值
            fpr, tpr, thresholds = roc_curve(y, X)

            # 计算Youden指数
            youden_index = tpr - fpr
            optimal_idx = np.argmax(youden_index)

            optimal_threshold = thresholds[optimal_idx]
            sensitivity = tpr[optimal_idx]
            specificity = 1 - fpr[optimal_idx]

            # 计算95%置信区间（bootstrap方法）
            n_bootstraps = 1000
            bootstrap_scores = []

            for i in range(n_bootstraps):
                # bootstrap抽样
                indices = resample(range(len(X)), replace=True, n_samples=len(X))
                if len(np.unique(y[indices])) < 2:
                    continue  # 需要两个类别

                try:
                    score = roc_auc_score(y[indices], X[indices])
                    bootstrap_scores.append(score)
                except:
                    continue

            if bootstrap_scores:
                sorted_scores = np.array(bootstrap_scores)
                sorted_scores.sort()
                ci_lower = np.percentile(sorted_scores, 2.5)
                ci_upper = np.percentile(sorted_scores, 97.5)
            else:
                ci_lower = ci_upper = np.nan

            # 格式化指标名称
            indicator_name = indicator
            if indicator == 'Retina_外环颞侧':
                display_name = 'Outer temporal thickness'
            elif indicator == 'Retina_平均厚度':
                display_name = 'Mean macular thickness'
            elif indicator == 'Retina_内环颞侧':
                display_name = 'Inner temporal thickness'
            elif indicator == 'Retina_总体积':
                display_name = 'Total macular volume'
            elif indicator == 'Retina_外环下方':
                display_name = 'Outer inferior thickness'
            elif indicator == 'RNFL_平均厚度':
                display_name = 'Peripapillary RNFL (average)'
            elif indicator == 'GCL+_平均厚度':
                display_name = 'GCL+ thickness'
            elif indicator == 'GCL++_平均厚度':
                display_name = 'GCL++ thickness'
            elif indicator == 'Choroid_平均厚度':
                display_name = 'Mean choroidal thickness'
            elif indicator == 'Rim Volume':
                display_name = 'Rim volume'
            elif indicator == 'C/D Area Ratio':
                display_name = 'Cup-to-disc ratio'
            else:
                display_name = indicator

            # 添加结果
            roc_results.append({
                'Parameter': display_name,
                'Original_name': indicator,
                'Sample_size': len(temp_df),
                'AUC': auc_score,
                'AUC_CI_lower': ci_lower,
                'AUC_CI_upper': ci_upper,
                'Sensitivity': sensitivity,
                'Specificity': specificity,
                'Youden_Index': youden_index[optimal_idx],
                'Optimal_Cutoff': optimal_threshold
            })

            print(f"  {indicator}: AUC={auc_score:.3f}, n={len(temp_df)}, 敏感度={sensitivity:.3f}, 特异度={specificity:.3f}")

        except Exception as e:
            print(f"  错误计算 {indicator}: {e}")

    # 按AUC降序排序
    roc_results.sort(key=lambda x: x['AUC'], reverse=True)

    # 4. 生成Table 5数据
    print(f"\n📋 生成Table 5数据...")

    # 创建DataFrame
    roc_df = pd.DataFrame(roc_results)

    # 选择前6个最佳指标（与论文一致）
    top_roc_df = roc_df.head(6).copy()

    # 格式化显示
    table5_rows = []

    for idx, row in top_roc_df.iterrows():
        # 格式化AUC和置信区间
        auc_str = f"{row['AUC']:.3f}"
        ci_str = f"{row['AUC_CI_lower']:.3f}–{row['AUC_CI_upper']:.3f}"

        # 确定单位
        if 'thickness' in row['Parameter'].lower():
            unit = 'μm'
            cutoff = f"{row['Optimal_Cutoff']:.2f} {unit}"
        elif 'volume' in row['Parameter'].lower():
            unit = 'mm³'
            cutoff = f"{row['Optimal_Cutoff']:.2f} {unit}"
        elif 'ratio' in row['Parameter'].lower():
            unit = ''
            cutoff = f"{row['Optimal_Cutoff']:.3f}"
        else:
            unit = ''
            cutoff = f"{row['Optimal_Cutoff']:.2f}"

        table5_rows.append({
            'Parameter': row['Parameter'],
            'AUC_CI': f"{auc_str} ({ci_str})",
            'Sensitivity_percent': f"{row['Sensitivity']*100:.1f}",
            'Specificity_percent': f"{row['Specificity']*100:.1f}",
            'Youden_Index': f"{row['Youden_Index']:.3f}",
            'Optimal_Cutoff': cutoff,
            'Sample_size': row['Sample_size']
        })

    # 5. 生成Markdown格式的Table 5
    print(f"\n📝 生成Markdown格式的Table 5...")

    md_table5 = "## Table 5. Diagnostic performance of OCT parameters (based on 463 eyes with complete age and sex data)\n\n"
    md_table5 += "| Parameter | AUC (95% CI) | Sensitivity (%) | Specificity (%) | Youden Index | Optimal Cutoff |\n"
    md_table5 += "|-----------|--------------|-----------------|-----------------|--------------|----------------|\n"

    for row in table5_rows:
        md_table5 += f"| {row['Parameter']} | {row['AUC_CI']} | {row['Sensitivity_percent']} | {row['Specificity_percent']} | {row['Youden_Index']} | {row['Optimal_Cutoff']} |\n"

    # 添加脚注
    md_table5 += f"\n*Note: Analysis based on {roc_df['Sample_size'].min()}-{roc_df['Sample_size'].max()} eyes with complete OCT and demographic data from the 463-eye sample. AUC = area under the receiver operating characteristic curve; CI = confidence interval.*\n"

    # 6. 保存结果
    print(f"\n💾 保存结果...")

    # 保存完整ROC结果
    roc_output_file = os.path.join(tables_dir, f"ROC_Analysis_463eyes_{datetime.now().strftime('%Y%m%d')}.xlsx")
    roc_df.to_excel(roc_output_file, index=False)
    print(f"✅ 完整ROC结果: {roc_output_file}")

    # 保存Table 5数据
    table5_output_file = os.path.join(output_dir, f"Table5_ROC_Analysis_463eyes_{datetime.now().strftime('%Y%m%d')}.xlsx")
    top_roc_df.to_excel(table5_output_file, index=False)
    print(f"✅ Table 5数据: {table5_output_file}")

    # 保存Markdown格式
    table5_md_file = os.path.join(output_dir, f"Table5_ROC_Markdown_{datetime.now().strftime('%Y%m%d')}.md")
    with open(table5_md_file, 'w', encoding='utf-8') as f:
        f.write(md_table5)
    print(f"✅ Table 5 Markdown: {table5_md_file}")

    # 7. 显示关键结果
    print(f"\n🔑 关键结果摘要:")
    print("=" * 60)
    print(f"分析样本量: {roc_df['Sample_size'].min()}-{roc_df['Sample_size'].max()} 眼 (来自463眼样本)")
    print(f"最佳指标: {table5_rows[0]['Parameter']} (AUC: {table5_rows[0]['AUC_CI']})")
    print(f"平均AUC: {roc_df['AUC'].mean():.3f}")

    # 对比旧结果
    old_best_auc = 0.646  # 论文中旧的outer temporal thickness AUC
    new_best_auc = top_roc_df.iloc[0]['AUC']
    difference = new_best_auc - old_best_auc

    print(f"\n📊 与旧结果对比:")
    print(f"  旧最佳AUC: {old_best_auc:.3f}")
    print(f"  新最佳AUC: {new_best_auc:.3f}")
    print(f"  差异: {difference:.3f} ({difference/old_best_auc*100:.1f}%)")

    # 8. 生成更新说明
    summary = f"""# ROC表格数据更新说明
    ## 基于463眼版本重新计算
    ### 更新日期: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

    ---

    ## 📊 更新概览

    ### 样本基础
    - **总样本**: 463眼 (有完整年龄性别数据)
    - **ROC分析样本**: {roc_df['Sample_size'].min()}-{roc_df['Sample_size'].max()}眼 (各指标完整数据)
    - **对比旧版本**: 原基于251眼，现基于463眼样本

    ### 关键结果变化
    | 指标 | 旧AUC (251眼) | 新AUC (463眼) | 变化 |
    |------|---------------|---------------|------|
    | **最佳指标** | Outer temporal thickness (0.646) | {table5_rows[0]['Parameter']} ({table5_rows[0]['AUC_CI'].split('(')[0].strip()}) | {difference:.3f} |
    | **平均AUC** | ~0.620 | {roc_df['AUC'].mean():.3f} | {roc_df['AUC'].mean()-0.620:.3f} |
    | **样本量** | 251眼 | {int(roc_df['Sample_size'].mean()):.0f}眼 | +{int(roc_df['Sample_size'].mean()-251):.0f}眼 |

    ### 新Table 5数据
    {md_table5}

    ### 方法学改进
    1. **样本一致性**: 所有分析基于相同463眼样本
    2. **数据完整性**: 排除缺失年龄/性别的病例  
    3. **统计稳健性**: Bootstrap方法计算95%置信区间
    4. **可比性提高**: 与传统ROC分析保持一致，便于与机器学习结果对比

    ### 对论文的影响
    1. **需要更新**: Table 5中的ROC数据
    2. **需要更新**: 3.6 Diagnostic Performance章节中的具体数值
    3. **需要确认**: 与3.7机器学习结果的一致性
    4. **建议更新**: 图片引用部分(Table 5引用)

    ### 生成文件
    1. `ROC_Analysis_463eyes_{datetime.now().strftime('%Y%m%d')}.xlsx` - 完整ROC结果
    2. `Table5_ROC_Analysis_463eyes_{datetime.now().strftime('%Y%m%d')}.xlsx` - Table 5数据
    3. `Table5_ROC_Markdown_{datetime.now().strftime('%Y%m%d')}.md` - 可直接插入论文的Table 5

    ---

    ## 🎯 论文更新建议

    ### 立即更新内容
    1. **替换Table 5**: 使用新生成的Table 5数据
    2. **更新3.6章节**: 更新AUC、敏感度、特异度等具体数值
    3. **更新样本量说明**: 明确说明基于463眼样本

    ### 审稿应对准备
    - **解释样本量变化**: 从251眼到463眼的方法学改进
    - **保持一致性**: 确保传统ROC分析与机器学习结果部分逻辑一致
    - **透明度**: 报告数据筛选过程和样本量变化

    ---
    *更新完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
    *数据来源: data_499eyes_20260315.xlsx → 463眼有完整年龄性别数据样本*
    """

    summary_file = os.path.join(output_dir, f"ROC_Table_Update_Summary_{datetime.now().strftime('%Y%m%d')}.md")
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write(summary)

    print(f"✅ 更新说明: {summary_file}")

    print(f"\n" + "=" * 80)
    print("🎉 ROC表格数据更新完成!")
    print("基于463眼版本的ROC分析数据已生成")
    print("需要更新论文中的Table 5和3.6章节")
    print("=" * 80)


if __name__ == "__main__":
    main()