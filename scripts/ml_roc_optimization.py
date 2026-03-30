#!/usr/bin/env python3

# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *

def main():
    """主函数，包装原有执行代码"""
    """
    机器学习模型优化ROC分析 - 使用交叉验证提高诊断性能
    目标：通过机器学习模型和交叉验证优化抑郁诊断的ROC性能
    """

    import pandas as pd
    import numpy as np
    import os
    from datetime import datetime
    import warnings
    warnings.filterwarnings('ignore')

    print("=" * 80)
    print("机器学习模型优化ROC分析")
    print(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    # 路径设置
    base_dir = "/mnt/c/Users/CUI/Desktop/投稿、数据修改"
    data_path = os.path.join(base_dir, "04_Data", "data_499eyes_20260315.xlsx")
    tables_dir = os.path.join(base_dir, "03_Tables")
    output_dir = os.path.join(base_dir, "分析报告", "02_分析报告文档")

    # 确保目录存在
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # 1. 数据准备
    print(f"\n📊 数据准备...")
    df = pd.read_excel(data_path)
    print(f"原始数据: {df.shape[0]} 行 × {df.shape[1]} 列")

    # 筛选463眼有完整年龄性别数据的样本
    df_463 = df.dropna(subset=['年龄', '性别']).copy()
    print(f"463眼样本: {df_463.shape[0]} 眼")
    print(f"分组分布: 抑郁症 {df_463[df_463['分组'] == '抑郁症'].shape[0]}眼, "
          f"健康对照 {df_463[df_463['分组'] == '健康对照'].shape[0]}眼")

    # 识别OCT指标列
    print(f"\n🔍 识别OCT指标...")
    oct_columns = []
    for col in df_463.columns:
        if any(keyword in col for keyword in ['Retina', 'RNFL', 'GCL', 'Choroid', 'Disc', 'Cup', 'Rim']):
            oct_columns.append(col)

    print(f"找到 {len(oct_columns)} 个OCT指标")

    # 检查每个OCT指标的缺失情况
    missing_info = []
    for col in oct_columns[:15]:  # 只检查前15个
        missing = df_463[col].isna().sum()
        pct = missing / df_463.shape[0] * 100
        if missing > 0:
            missing_info.append((col, missing, pct))

    if missing_info:
        print(f"部分OCT指标缺失情况:")
        for col, missing, pct in missing_info[:10]:
            print(f"  {col}: 缺失 {missing}眼 ({pct:.1f}%)")

    # 选择缺失最少的指标进行机器学习分析
    # 计算每个指标的缺失比例
    missing_rates = {}
    for col in oct_columns:
        missing_rate = df_463[col].isna().sum() / df_463.shape[0]
        missing_rates[col] = missing_rate

    # 选择缺失率<5%的指标
    selected_features = [col for col, rate in missing_rates.items() if rate < 0.05]
    print(f"\n✅ 选择 {len(selected_features)} 个缺失率<5%的OCT指标进行机器学习分析")

    # 创建完整数据集（所有选定特征都非缺失）
    df_complete = df_463.dropna(subset=selected_features + ['分组'])
    print(f"完整数据集: {df_complete.shape[0]} 眼")
    print(f"  抑郁症组: {df_complete[df_complete['分组'] == '抑郁症'].shape[0]} 眼")
    print(f"  健康对照组: {df_complete[df_complete['分组'] == '健康对照'].shape[0]} 眼")

    # 准备特征和目标变量
    X = df_complete[selected_features]
    y = (df_complete['分组'] == '抑郁症').astype(int)  # 1=抑郁, 0=对照

    print(f"\n📈 特征矩阵: {X.shape[0]} 样本 × {X.shape[1]} 特征")
    print(f"目标变量: 抑郁症 {y.sum()}眼, 健康对照 {len(y)-y.sum()}眼")

    # 2. 机器学习模型导入
    print(f"\n🤖 导入机器学习库...")
    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.svm import SVC
        from sklearn.neighbors import KNeighborsClassifier
        from sklearn.tree import DecisionTreeClassifier
        from sklearn.naive_bayes import GaussianNB
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_val_score
        from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, precision_score, recall_score, f1_score
        from sklearn.pipeline import Pipeline
        import xgboost as xgb
        xgb_available = True
        print("✅ XGBoost 可用")
    except ImportError as e:
        print(f"⚠️ 部分库不可用: {e}")
        xgb_available = False

    # 3. 定义机器学习模型
    models = {
        '逻辑回归 (L2正则化)': LogisticRegression(penalty='l2', C=1.0, max_iter=1000, random_state=42),
        '随机森林': RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42),
        '支持向量机 (线性)': SVC(kernel='linear', probability=True, random_state=42),
        'K近邻 (K=5)': KNeighborsClassifier(n_neighbors=5),
        '决策树': DecisionTreeClassifier(max_depth=5, random_state=42),
        '朴素贝叶斯': GaussianNB(),
    }

    if xgb_available:
        models['XGBoost'] = xgb.XGBClassifier(n_estimators=100, max_depth=3, learning_rate=0.1, random_state=42)

    # 4. 5折交叉验证
    print(f"\n🔄 开始5折分层交叉验证...")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    results = []
    roc_results = {}

    for model_name, model in models.items():
        print(f"  训练 {model_name}...")

        # 创建包含标准化的pipeline
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', model)
        ])

        try:
            # 交叉验证预测概率
            y_pred_proba = cross_val_predict(pipeline, X, y, cv=cv, method='predict_proba')[:, 1]

            # 计算AUC
            auc = roc_auc_score(y, y_pred_proba)

            # 计算其他指标
            y_pred = (y_pred_proba >= 0.5).astype(int)
            accuracy = accuracy_score(y, y_pred)
            precision = precision_score(y, y_pred)
            recall = recall_score(y, y_pred)
            f1 = f1_score(y, y_pred)

            # ROC曲线数据
            fpr, tpr, thresholds = roc_curve(y, y_pred_proba)
            youden_idx = np.argmax(tpr - fpr)
            optimal_threshold = thresholds[youden_idx]
            sensitivity = tpr[youden_idx]
            specificity = 1 - fpr[youden_idx]

            results.append({
                '模型': model_name,
                'AUC': auc,
                '敏感度': sensitivity,
                '特异度': specificity,
                '准确率': accuracy,
                '精确率': precision,
                '召回率': recall,
                'F1分数': f1,
                '最佳阈值': optimal_threshold
            })

            roc_results[model_name] = {
                'fpr': fpr,
                'tpr': tpr,
                'thresholds': thresholds,
                'auc': auc
            }

            print(f"    AUC: {auc:.3f}, 敏感度: {sensitivity:.3f}, 特异度: {specificity:.3f}")

        except Exception as e:
            print(f"   ⚠️ {model_name} 训练失败: {e}")

    # 5. 结果分析
    print(f"\n📊 机器学习模型性能对比")
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('AUC', ascending=False)
    print(results_df.to_string(index=False))

    # 6. 与单一指标对比
    print(f"\n🔍 与单一OCT指标对比...")

    # 计算每个单一指标的AUC
    single_auc_results = []
    for feature in selected_features[:10]:  # 只测试前10个指标
        try:
            # 使用简单的逻辑回归（仅该特征）
            X_single = df_complete[[feature]]
            model_single = LogisticRegression(penalty='l2', C=1.0, max_iter=1000, random_state=42)

            # 交叉验证
            y_pred_single = cross_val_predict(
                Pipeline([('scaler', StandardScaler()), ('clf', model_single)]),
                X_single, y, cv=cv, method='predict_proba'
            )[:, 1]

            auc_single = roc_auc_score(y, y_pred_single)

            # ROC曲线
            fpr_s, tpr_s, _ = roc_curve(y, y_pred_single)
            youden_idx_s = np.argmax(tpr_s - fpr_s)
            sensitivity_s = tpr_s[youden_idx_s]
            specificity_s = 1 - fpr_s[youden_idx_s]

            single_auc_results.append({
                '指标': feature,
                'AUC': auc_single,
                '敏感度': sensitivity_s,
                '特异度': specificity_s
            })
        except:
            continue

    single_df = pd.DataFrame(single_auc_results)
    single_df = single_df.sort_values('AUC', ascending=False)

    print(f"单一指标最佳AUC: {single_df['AUC'].max():.3f} ({single_df.iloc[0]['指标']})")
    print(f"机器学习最佳AUC: {results_df['AUC'].max():.3f} ({results_df.iloc[0]['模型']})")
    print(f"AUC提升: {results_df['AUC'].max() - single_df['AUC'].max():.3f}")

    # 7. 特征重要性分析（随机森林）
    print(f"\n🎯 特征重要性分析...")
    rf_model = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
    rf_pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', rf_model)
    ])
    rf_pipeline.fit(X, y)

    feature_importance = pd.DataFrame({
        '特征': selected_features,
        '重要性': rf_pipeline.named_steps['classifier'].feature_importances_
    })
    feature_importance = feature_importance.sort_values('重要性', ascending=False)

    print(f"最重要的10个特征:")
    print(feature_importance.head(10).to_string(index=False))

    # 8. 复合指标创建
    print(f"\n🔧 创建复合指标...")
    # 使用逻辑回归权重创建复合指标
    lr_model = LogisticRegression(penalty='l2', C=1.0, max_iter=1000, random_state=42)
    lr_pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', lr_model)
    ])
    lr_pipeline.fit(X, y)

    # 获取标准化后的权重
    scaler = lr_pipeline.named_steps['scaler']
    X_scaled = scaler.transform(X)
    lr_final = lr_pipeline.named_steps['classifier']

    # 复合指标公式: weighted sum of standardized features
    composite_weights = lr_final.coef_[0]
    composite_score = X_scaled.dot(composite_weights)

    # 计算复合指标的ROC
    from sklearn.metrics import roc_auc_score, roc_curve
    auc_composite = roc_auc_score(y, composite_score)
    fpr_c, tpr_c, thresholds_c = roc_curve(y, composite_score)
    youden_idx_c = np.argmax(tpr_c - fpr_c)
    optimal_threshold_c = thresholds_c[youden_idx_c]
    sensitivity_c = tpr_c[youden_idx_c]
    specificity_c = 1 - fpr_c[youden_idx_c]

    print(f"复合指标性能:")
    print(f"  AUC: {auc_composite:.3f}")
    print(f"  敏感度: {sensitivity_c:.3f}")
    print(f"  特异度: {specificity_c:.3f}")
    print(f"  最佳阈值: {optimal_threshold_c:.3f}")

    # 9. 保存结果
    print(f"\n💾 保存结果...")
    date_str = datetime.now().strftime("%Y%m%d")

    # 保存模型性能对比
    model_perf_file = os.path.join(tables_dir, f"机器学习模型性能对比_{date_str}.xlsx")
    results_df.to_excel(model_perf_file, index=False)
    print(f"✅ 模型性能对比: {model_perf_file}")

    # 保存单一指标对比
    single_auc_file = os.path.join(tables_dir, f"单一指标AUC对比_{date_str}.xlsx")
    single_df.to_excel(single_auc_file, index=False)
    print(f"✅ 单一指标对比: {single_auc_file}")

    # 保存特征重要性
    feature_imp_file = os.path.join(tables_dir, f"特征重要性分析_{date_str}.xlsx")
    feature_importance.to_excel(feature_imp_file, index=False)
    print(f"✅ 特征重要性: {feature_imp_file}")

    # 保存复合指标权重
    composite_weights_df = pd.DataFrame({
        '特征': selected_features,
        '权重': composite_weights,
        '绝对值权重': np.abs(composite_weights)
    })
    composite_weights_df = composite_weights_df.sort_values('绝对值权重', ascending=False)
    composite_file = os.path.join(tables_dir, f"复合指标权重_{date_str}.xlsx")
    composite_weights_df.to_excel(composite_file, index=False)
    print(f"✅ 复合指标权重: {composite_file}")

    # 10. 生成详细报告
    print(f"\n📋 生成分析报告...")

    report_text = f"""# 机器学习模型优化ROC分析报告
    ## 使用交叉验证提高抑郁诊断性能
    ### 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

    ---

    ## 📊 数据概况

    ### 样本信息
    - **原始数据**: {df.shape[0]} 眼
    - **463眼样本**: {df_463.shape[0]} 眼（有完整年龄性别数据）
    - **完整数据集**: {df_complete.shape[0]} 眼（所有选定OCT指标非缺失）
    - **特征数量**: {len(selected_features)} 个OCT指标（缺失率<5%）

    ### 分组分布
    - **抑郁症组**: {df_complete[df_complete['分组'] == '抑郁症'].shape[0]} 眼
    - **健康对照组**: {df_complete[df_complete['分组'] == '健康对照'].shape[0]} 眼

    ---

    ## 🤖 机器学习模型性能

    ### 模型对比（5折交叉验证）
    {results_df.to_string(index=False)}

    ### 关键发现
    1. **最佳模型**: {results_df.iloc[0]['模型']} (AUC={results_df.iloc[0]['AUC']:.3f})
    2. **AUC范围**: {results_df['AUC'].min():.3f} - {results_df['AUC'].max():.3f}
    3. **敏感度-特异度平衡**: 
       - 最佳敏感度: {results_df['敏感度'].max():.3f} ({results_df.loc[results_df['敏感度'].idxmax(), '模型']})
       - 最佳特异度: {results_df['特异度'].max():.3f} ({results_df.loc[results_df['特异度'].idxmax(), '模型']})

    ---

    ## 🔍 与单一指标对比

    ### 单一指标最佳性能
    {single_df.head(5).to_string(index=False)}

    ### 性能提升分析
    - **单一指标最佳AUC**: {single_df['AUC'].max():.3f} ({single_df.iloc[0]['指标']})
    - **机器学习最佳AUC**: {results_df['AUC'].max():.3f} ({results_df.iloc[0]['模型']})
    - **AUC绝对提升**: {results_df['AUC'].max() - single_df['AUC'].max():.3f}
    - **AUC相对提升**: {(results_df['AUC'].max() - single_df['AUC'].max()) / single_df['AUC'].max() * 100:.1f}%

    ---

    ## 🎯 特征重要性分析

    ### 最重要的10个特征
    {feature_importance.head(10).to_string(index=False)}

    ### 特征类别分布
    ```python
    # 特征分类统计
    视网膜厚度指标: {len([f for f in selected_features if 'Retina' in f])} 个
    RNFL指标: {len([f for f in selected_features if 'RNFL' in f])} 个
    视盘参数: {len([f for f in selected_features if any(kw in f for kw in ['Disc', 'Cup', 'Rim'])])} 个
    ```

    ---

    ## 🔧 复合指标分析

    ### 逻辑回归加权复合指标
    - **AUC**: {auc_composite:.3f}
    - **敏感度**: {sensitivity_c:.3f}
    - **特异度**: {specificity_c:.3f}
    - **最佳阈值**: {optimal_threshold_c:.3f}

    ### 复合指标公式
    ```
    复合分数 = Σ(权重ᵢ × 标准化特征ᵢ)
    ```

    ### 最重要的5个权重
    {composite_weights_df.head(5).to_string(index=False)}

    ---

    ## 📈 临床意义分析

    ### 诊断性能评估
    1. **当前最佳单一指标**: Retina_外环颞侧 (AUC≈0.650)
    2. **机器学习优化后**: AUC提升至 {results_df['AUC'].max():.3f}
    3. **临床可接受性**: 
       - AUC>0.70通常被认为具有"可接受的"诊断性能
       - 当前最佳: {results_df['AUC'].max():.3f} ({'达到' if results_df['AUC'].max() > 0.70 else '未达到'}可接受标准)

    ### 敏感度-特异度权衡
    - **高敏感度策略**: 适合筛查，减少漏诊
    - **高特异度策略**: 适合确诊辅助，减少误诊
    - **平衡策略**: Youden指数最大化（当前采用）

    ---

    ## ⚠️ 方法学注意事项

    ### 交叉验证优势
    1. **避免过拟合**: 5折分层交叉验证确保结果稳健
    2. **偏差-方差权衡**: 平衡模型复杂度和泛化能力
    3. **可重复性**: 固定随机种子确保结果可重复

    ### 模型选择考虑
    1. **线性模型**: 逻辑回归易于解释，适合临床转化
    2. **非线性模型**: 随机森林/XGBoost可能捕获复杂模式但可解释性差
    3. **过拟合风险**: 特征数({len(selected_features)}) vs 样本数({df_complete.shape[0]})

    ### 数据限制
    1. **样本量**: n={df_complete.shape[0]}，对于复杂模型可能偏小
    2. **类别平衡**: 抑郁症:健康对照 ≈ {df_complete[df_complete['分组'] == '抑郁症'].shape[0]}:{df_complete[df_complete['分组'] == '健康对照'].shape[0]}
    3. **缺失数据**: 使用完整病例分析，可能引入选择偏倚

    ---

    ## 🎯 优化建议

    ### 立即实施
    1. **采用复合指标**: 逻辑回归加权复合指标 (AUC={auc_composite:.3f})
    2. **阈值优化**: 根据临床需求调整敏感度/特异度平衡
    3. **特征简化**: 使用最重要的5-10个特征以减少过拟合风险

    ### 进一步优化
    1. **特征工程**: 创建交互项、比值特征等
    2. **集成方法**: 模型堆叠或投票集成
    3. **超参数调优**: 网格搜索优化模型参数
    4. **外部验证**: 在独立数据集上验证性能

    ### 临床转化
    1. **简化模型**: 选择最重要的3-5个特征创建临床实用评分
    2. **阈值设定**: 根据临床场景设定不同阈值
    3. **验证研究**: 前瞻性研究验证诊断性能

    ---

    ## 📁 输出文件清单

    ### Excel文件
    1. `机器学习模型性能对比_{date_str}.xlsx` - 所有模型性能对比
    2. `单一指标AUC对比_{date_str}.xlsx` - 单一指标性能
    3. `特征重要性分析_{date_str}.xlsx` - 随机森林特征重要性
    4. `复合指标权重_{date_str}.xlsx` - 逻辑回归复合指标权重

    ### 关键结果摘要
    - **最佳模型**: {results_df.iloc[0]['模型']}
    - **最佳AUC**: {results_df.iloc[0]['AUC']:.3f}
    - **AUC提升**: {results_df['AUC'].max() - single_df['AUC'].max():.3f}
    - **推荐方案**: 逻辑回归复合指标 (AUC={auc_composite:.3f})

    ---

    ## 📞 使用建议

    ### 论文撰写
    - **方法部分**: 描述使用的机器学习方法和交叉验证策略
    - **结果部分**: 报告最佳模型性能和与单一指标的对比
    - **讨论部分**: 讨论机器学习在OCT诊断中的潜在价值

    ### 审稿准备
    - **优势强调**: 交叉验证确保结果稳健性
    - **透明报告**: 详细报告所有模型性能
    - **局限性**: 承认样本量限制和过拟合风险

    ### 临床应用考虑
    - **可解释性**: 逻辑回归模型权重易于临床解释
    - **实用性**: 复合指标可转化为临床评分系统
    - **验证需求**: 强调需要外部验证研究

    ---
    *分析完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
    *交叉验证策略: 5折分层交叉验证*
    *样本规模: {df_complete.shape[0]}眼 (所有选定OCT指标完整)*
    """

    # 保存报告
    report_file = os.path.join(output_dir, f"机器学习优化ROC分析报告_{date_str}.md")
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report_text)

    print(f"✅ 详细分析报告: {report_file}")

    print(f"\n" + "=" * 80)
    print("🎉 机器学习模型优化完成!")
    print(f"最佳模型: {results_df.iloc[0]['模型']} (AUC={results_df.iloc[0]['AUC']:.3f})")
    print(f"AUC提升: {results_df['AUC'].max() - single_df['AUC'].max():.3f}")
    print("=" * 80)

    # 11. 生成简化的临床评分系统
    print(f"\n🔬 生成简化临床评分系统...")

    # 选择最重要的5个特征
    top_features = feature_importance.head(5)['特征'].tolist()
    print(f"最重要的5个特征: {top_features}")

    # 使用这5个特征重新训练简化模型
    X_simple = df_complete[top_features]
    lr_simple = LogisticRegression(penalty='l2', C=1.0, max_iter=1000, random_state=42)
    lr_simple_pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', lr_simple)
    ])

    # 交叉验证评估简化模型
    y_pred_simple = cross_val_predict(lr_simple_pipeline, X_simple, y, cv=cv, method='predict_proba')[:, 1]
    auc_simple = roc_auc_score(y, y_pred_simple)
    print(f"简化模型(5个特征) AUC: {auc_simple:.3f}")

    # 训练最终模型获取权重
    lr_simple_pipeline.fit(X_simple, y)
    scaler_simple = lr_simple_pipeline.named_steps['scaler']
    lr_final_simple = lr_simple_pipeline.named_steps['classifier']

    # 生成临床评分公式
    print(f"\n📝 临床评分公式:")
    print("抑郁风险评分 = ")
    for i, (feature, weight) in enumerate(zip(top_features, lr_final_simple.coef_[0])):
        sign = "+" if weight >= 0 else "-"
        print(f"    {sign} {abs(weight):.3f} × (标准化{feature})")

    print(f"\n使用说明:")
    print("1. 对每个特征进行标准化: (原始值 - 均值)/标准差")
    print("2. 计算加权和得到抑郁风险评分")
    print("3. 与最佳阈值比较进行诊断")

    # 保存简化模型信息
    simple_model_file = os.path.join(tables_dir, f"简化临床评分系统_{date_str}.xlsx")
    simple_model_df = pd.DataFrame({
        '特征': top_features,
        '原始均值': [df_complete[feat].mean() for feat in top_features],
        '原始标准差': [df_complete[feat].std() for feat in top_features],
        '标准化权重': lr_final_simple.coef_[0],
        '重要性排名': range(1, 6)
    })
    simple_model_df.to_excel(simple_model_file, index=False)
    print(f"✅ 简化临床评分系统: {simple_model_file}")


if __name__ == "__main__":
    main()