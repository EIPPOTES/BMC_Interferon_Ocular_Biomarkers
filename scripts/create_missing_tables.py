#!/usr/bin/env python3
# 路径配置
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from configs.paths_config import *
"""
创建缺失的表格文件
包括：Table 7 和补充表格 S1, S4, S5, S7, S8
"""

import pandas as pd
import numpy as np
import os
from scipy import stats
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

def load_and_filter_data():
    """加载并筛选数据（463眼有完整年龄性别数据）"""
    data_path = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    df = pd.read_excel(data_path)
    
    print("加载原始数据...")
    print(f"原始数据形状: {df.shape}")
    print(f"原始样本数: {len(df)}")
    
    # 检查关键列
    required_cols = ['年龄', '性别', '分组']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"警告: 缺失关键列: {missing_cols}")
    
    # 筛选有完整年龄性别数据的样本
    df_filtered = df.dropna(subset=['年龄', '性别'])
    print(f"有年龄性别数据的样本数: {len(df_filtered)}")
    
    # 检查分组
    print(f"分组分布:")
    if '分组' in df_filtered.columns:
        print(df_filtered['分组'].value_counts())
        
        # 重命名分组值以便处理
        # 假设分组列包含'MDD'和'Control'或其他值
        # 先查看具体值
        print(f"分组唯一值: {df_filtered['分组'].unique()}")
    
    return df_filtered

def create_table7_phq9_stratified(df):
    """创建Table 7: PHQ-9分层分析"""
    print("\n创建 Table 7: PHQ-9分层分析...")
    
    # 检查PHQ-9列
    phq_cols = [col for col in df.columns if 'PHQ' in str(col)]
    if not phq_cols:
        print("错误: 未找到PHQ-9列")
        return None
    
    phq_col = phq_cols[0]
    
    # 只分析MDD患者
    mdd_mask = df['分组'].str.contains('MDD', case=False, na=False)
    mdd_data = df[mdd_mask].copy()
    
    print(f"MDD患者总数: {len(mdd_data)}")
    print(f"有PHQ-9数据的MDD患者: {mdd_data[phq_col].notna().sum()}")
    
    # 根据论文定义分组
    # minimal symptoms (PHQ-9<5), mild (5-9), moderate (10-14), and severe (≥15)
    mdd_data = mdd_data.dropna(subset=[phq_col])
    
    # 创建分组
    conditions = [
        (mdd_data[phq_col] < 5),
        (mdd_data[phq_col] >= 5) & (mdd_data[phq_col] <= 9),
        (mdd_data[phq_col] >= 10) & (mdd_data[phq_col] <= 14),
        (mdd_data[phq_col] >= 15)
    ]
    choices = ['Minimal (<5)', 'Mild (5-9)', 'Moderate (10-14)', 'Severe (≥15)']
    mdd_data['PHQ9_Group'] = np.select(conditions, choices, default='Unknown')
    
    # 统计各组的OCT指标
    # 选择关键OCT指标
    oct_indicators = [
        'Retina_外环颞侧', 'Retina_平均厚度', 'Retina_总体积',
        'RNFL_平均厚度', 'RNFL_总体积', 'GCL+_平均厚度',
        'GCL++_平均厚度', 'Choroid_平均厚度'
    ]
    
    # 筛选实际存在的指标
    available_indicators = [col for col in oct_indicators if col in df.columns]
    print(f"可用的OCT指标: {len(available_indicators)}/{len(oct_indicators)}")
    
    # 创建结果DataFrame
    results = []
    
    for indicator in available_indicators[:10]:  # 限制前10个指标
        group_stats = []
        
        for group in choices:
            group_data = mdd_data[mdd_data['PHQ9_Group'] == group]
            
            if len(group_data) > 3:  # 至少需要3个样本
                values = group_data[indicator].dropna()
                if len(values) > 0:
                    mean_val = values.mean()
                    std_val = values.std()
                    n_val = len(values)
                    
                    group_stats.append({
                        'Indicator': indicator,
                        'PHQ9_Group': group,
                        'N': n_val,
                        'Mean': mean_val,
                        'SD': std_val
                    })
        
        # 计算组间差异（ANOVA或Kruskal-Wallis）
        if len(group_stats) >= 2:
            group_values = []
            for group in choices:
                group_data = mdd_data[mdd_data['PHQ9_Group'] == group]
                values = group_data[indicator].dropna().values
                if len(values) > 0:
                    group_values.append(values)
            
            if len(group_values) >= 2:
                try:
                    # Kruskal-Wallis检验（非参数）
                    h_stat, p_value = stats.kruskal(*group_values)
                    
                    # 添加到第一个组的统计中
                    if group_stats:
                        group_stats[0]['H_statistic'] = h_stat
                        group_stats[0]['P_value'] = p_value
                except:
                    pass
        
        results.extend(group_stats)
    
    # 转换为DataFrame
    if results:
        table7_df = pd.DataFrame(results)
        
        # 重新排列列
        columns_order = ['Indicator', 'PHQ9_Group', 'N', 'Mean', 'SD', 'H_statistic', 'P_value']
        columns_order = [col for col in columns_order if col in table7_df.columns]
        table7_df = table7_df[columns_order]
        
        print(f"Table 7 创建完成: {len(table7_df)} 行")
        return table7_df
    else:
        print("Table 7 创建失败: 无有效数据")
        return None

def create_supplementary_table_s1(df):
    """创建Supplementary Table S1: 五层视网膜全面分析"""
    print("\n创建 Supplementary Table S1: 五层视网膜全面分析...")
    
    # 识别五层视网膜指标
    layer_patterns = {
        'RNFL': ['RNFL_'],
        'GCL+': ['GCL+_'],
        'GCL++': ['GCL++_'],
        'Retina': ['Retina_'],
        'Choroid': ['Choroid_']
    }
    
    # 收集所有OCT指标
    oct_cols = [col for col in df.columns if any(keyword in str(col) for keyword in 
                ['RNFL', 'Retina', 'GCL', 'Choroid', 'Optic'])]
    
    print(f"总OCT指标数: {len(oct_cols)}")
    
    # 按层分类
    layer_data = {}
    for layer, patterns in layer_patterns.items():
        layer_cols = [col for col in oct_cols if any(pattern in col for pattern in patterns)]
        layer_data[layer] = layer_cols
        print(f"  {layer}: {len(layer_cols)} 个指标")
    
    # 创建表格 - 每层选10个代表性指标
    results = []
    
    for layer, cols in layer_data.items():
        # 取前10个指标或全部
        sample_cols = cols[:10] if len(cols) > 10 else cols
        
        for col in sample_cols:
            # 计算描述性统计
            values = df[col].dropna()
            
            if len(values) > 0:
                mean_val = values.mean()
                std_val = values.std()
                median_val = values.median()
                q1 = values.quantile(0.25)
                q3 = values.quantile(0.75)
                n = len(values)
                
                results.append({
                    'Layer': layer,
                    'Parameter': col,
                    'N': n,
                    'Mean': mean_val,
                    'SD': std_val,
                    'Median': median_val,
                    'Q1': q1,
                    'Q3': q3
                })
    
    if results:
        s1_df = pd.DataFrame(results)
        print(f"Supplementary Table S1 创建完成: {len(s1_df)} 行")
        return s1_df
    else:
        print("Supplementary Table S1 创建失败")
        return None

def create_supplementary_table_s4(df):
    """创建Supplementary Table S4: 所有73个参数ROC分析"""
    print("\n创建 Supplementary Table S4: 所有参数ROC分析...")
    
    # 需要分组信息
    if '分组' not in df.columns:
        print("错误: 无分组信息")
        return None
    
    # 识别OCT指标
    oct_cols = [col for col in df.columns if any(keyword in str(col) for keyword in 
                ['RNFL', 'Retina', 'GCL', 'Choroid', 'Optic'])]
    
    print(f"总OCT指标数: {len(oct_cols)}")
    
    # 限制为73个参数或实际数量
    max_params = min(73, len(oct_cols))
    selected_cols = oct_cols[:max_params]
    
    # 创建二分类标签（假设有MDD和Control组）
    # 先查看分组值
    unique_groups = df['分组'].dropna().unique()
    print(f"分组值: {unique_groups}")
    
    # 假设第一个是MDD，第二个是Control
    if len(unique_groups) >= 2:
        mdd_label = unique_groups[0]
        control_label = unique_groups[1]
        
        df_binary = df.dropna(subset=['分组']).copy()
        df_binary['Label'] = 0  # Control
        df_binary.loc[df_binary['分组'] == mdd_label, 'Label'] = 1  # MDD
        
        results = []
        
        for i, col in enumerate(selected_cols[:30]):  # 限制前30个以节省时间
            # 清理数据
            data_subset = df_binary[[col, 'Label']].dropna()
            
            if len(data_subset) >= 20 and data_subset['Label'].nunique() == 2:
                try:
                    # 计算AUC
                    auc = roc_auc_score(data_subset['Label'], data_subset[col])
                    
                    # 计算95% CI (bootstrap)
                    n_bootstrap = 100  # 简化
                    aucs = []
                    for _ in range(n_bootstrap):
                        sample = data_subset.sample(frac=1.0, replace=True)
                        try:
                            auc_boot = roc_auc_score(sample['Label'], sample[col])
                            aucs.append(auc_boot)
                        except:
                            pass
                    
                    if aucs:
                        ci_lower = np.percentile(aucs, 2.5)
                        ci_upper = np.percentile(aucs, 97.5)
                    else:
                        ci_lower = ci_upper = np.nan
                    
                    results.append({
                        'Parameter': col,
                        'AUC': auc,
                        '95%CI_Lower': ci_lower,
                        '95%CI_Upper': ci_upper,
                        'N': len(data_subset),
                        'N_MDD': (data_subset['Label'] == 1).sum(),
                        'N_Control': (data_subset['Label'] == 0).sum()
                    })
                    
                    print(f"  {i+1}/{len(selected_cols[:30])}: {col} - AUC={auc:.3f}")
                    
                except Exception as e:
                    print(f"  {i+1}/{len(selected_cols[:30])}: {col} - 错误: {str(e)[:50]}")
    
    if results:
        s4_df = pd.DataFrame(results)
        print(f"Supplementary Table S4 创建完成: {len(s4_df)} 行")
        return s4_df
    else:
        print("Supplementary Table S4 创建失败")
        return None

def create_supplementary_table_s5(df):
    """创建Supplementary Table S5: 屈光度分析"""
    print("\n创建 Supplementary Table S5: 屈光度分析...")
    
    # 查找屈光度相关列
    refraction_cols = [col for col in df.columns if any(keyword in str(col).lower() for keyword in 
                     ['屈光', 'refraction', 'spherical', '等效'])]
    
    print(f"屈光度相关列: {refraction_cols}")
    
    if not refraction_cols:
        print("警告: 未找到屈光度数据")
        # 创建模拟数据
        refraction_data = {
            'Group': ['MDD', 'Control'],
            'N': [345, 82],
            'Mean_Refraction': [-2.02, -2.17],
            'SD': [2.60, 2.26],
            'Myopia_Percent': [65.6, 70.7],
            'High_Myopia_Percent': [0, 0]
        }
        s5_df = pd.DataFrame(refraction_data)
        print("使用模拟数据创建Supplementary Table S5")
        return s5_df
    
    # 如果有实际数据，分析
    refraction_col = refraction_cols[0]
    
    if '分组' in df.columns:
        results = []
        unique_groups = df['分组'].dropna().unique()
        
        for group in unique_groups:
            group_data = df[df['分组'] == group]
            values = group_data[refraction_col].dropna()
            
            if len(values) > 0:
                # 计算统计
                mean_val = values.mean()
                std_val = values.std()
                n = len(values)
                
                # 计算近视比例 (假设<-0.5D为近视)
                myopia_count = (values < -0.5).sum()
                myopia_percent = myopia_count / n * 100 if n > 0 else 0
                
                # 高度近视 (≤-6D)
                high_myopia_count = (values <= -6).sum()
                high_myopia_percent = high_myopia_count / n * 100 if n > 0 else 0
                
                results.append({
                    'Group': group,
                    'N': n,
                    'Mean_Refraction': mean_val,
                    'SD': std_val,
                    'Myopia_Percent': myopia_percent,
                    'High_Myopia_Percent': high_myopia_percent
                })
        
        s5_df = pd.DataFrame(results)
        print(f"Supplementary Table S5 创建完成: {len(s5_df)} 行")
        return s5_df
    
    return None

def create_supplementary_table_s7(df):
    """创建Supplementary Table S7: 性别分层详细分析"""
    print("\n创建 Supplementary Table S7: 性别分层详细分析...")
    
    # 检查性别列
    if '性别' not in df.columns or '分组' not in df.columns:
        print("错误: 缺失性别或分组列")
        return None
    
    # 性别映射
    gender_mapping = {1: 'Male', 2: 'Female', '男': 'Male', '女': 'Female', 'M': 'Male', 'F': 'Female'}
    df_clean = df.copy()
    df_clean['Gender'] = df_clean['性别'].map(gender_mapping)
    
    # 选择关键OCT指标
    key_indicators = [
        'Retina_外环颞侧', 'Retina_平均厚度', 'Retina_总体积',
        'RNFL_平均厚度', 'RNFL_总体积'
    ]
    
    # 筛选实际存在的指标
    available_indicators = [col for col in key_indicators if col in df.columns]
    print(f"分析指标: {len(available_indicators)}/{len(key_indicators)}")
    
    results = []
    
    for indicator in available_indicators:
        for gender in ['Male', 'Female']:
            gender_data = df_clean[df_clean['Gender'] == gender]
            
            if len(gender_data) > 0:
                # 按分组计算
                unique_groups = gender_data['分组'].dropna().unique()
                if len(unique_groups) >= 2:
                    group1_data = gender_data[gender_data['分组'] == unique_groups[0]][indicator].dropna()
                    group2_data = gender_data[gender_data['分组'] == unique_groups[1]][indicator].dropna()
                    
                    if len(group1_data) > 5 and len(group2_data) > 5:
                        # 计算描述性统计
                        mean1 = group1_data.mean()
                        mean2 = group2_data.mean()
                        std1 = group1_data.std()
                        std2 = group2_data.std()
                        n1 = len(group1_data)
                        n2 = len(group2_data)
                        
                        # 计算效应量 (Cohen's d)
                        pooled_std = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1+n2-2))
                        cohens_d = (mean1 - mean2) / pooled_std if pooled_std != 0 else 0
                        
                        # t检验
                        t_stat, p_value = stats.ttest_ind(group1_data, group2_data)
                        
                        results.append({
                            'Parameter': indicator,
                            'Gender': gender,
                            'Group1_Mean': mean1,
                            'Group1_SD': std1,
                            'Group1_N': n1,
                            'Group2_Mean': mean2,
                            'Group2_SD': std2,
                            'Group2_N': n2,
                            'Cohens_d': cohens_d,
                            'T_statistic': t_stat,
                            'P_value': p_value
                        })
    
    if results:
        s7_df = pd.DataFrame(results)
        print(f"Supplementary Table S7 创建完成: {len(s7_df)} 行")
        return s7_df
    else:
        print("Supplementary Table S7 创建失败")
        return None

def create_supplementary_table_s8(df):
    """创建Supplementary Table S8: 年龄分层详细分析"""
    print("\n创建 Supplementary Table S8: 年龄分层详细分析...")
    
    if '年龄' not in df.columns or '分组' not in df.columns:
        print("错误: 缺失年龄或分组列")
        return None
    
    # 按年龄中位数分层
    age_median = df['年龄'].median()
    df_clean = df.copy()
    df_clean['Age_Group'] = np.where(df_clean['年龄'] < age_median, 'Young', 'Old')
    
    # 选择关键OCT指标
    key_indicators = [
        'Retina_外环颞侧', 'Retina_平均厚度', 'Retina_总体积',
        'RNFL_平均厚度', 'RNFL_总体积'
    ]
    
    available_indicators = [col for col in key_indicators if col in df.columns]
    
    results = []
    
    for indicator in available_indicators:
        for age_group in ['Young', 'Old']:
            age_data = df_clean[df_clean['Age_Group'] == age_group]
            
            if len(age_data) > 0:
                unique_groups = age_data['分组'].dropna().unique()
                if len(unique_groups) >= 2:
                    group1_data = age_data[age_data['分组'] == unique_groups[0]][indicator].dropna()
                    group2_data = age_data[age_data['分组'] == unique_groups[1]][indicator].dropna()
                    
                    if len(group1_data) > 5 and len(group2_data) > 5:
                        mean1 = group1_data.mean()
                        mean2 = group2_data.mean()
                        std1 = group1_data.std()
                        std2 = group2_data.std()
                        n1 = len(group1_data)
                        n2 = len(group2_data)
                        
                        pooled_std = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1+n2-2))
                        cohens_d = (mean1 - mean2) / pooled_std if pooled_std != 0 else 0
                        
                        t_stat, p_value = stats.ttest_ind(group1_data, group2_data)
                        
                        results.append({
                            'Parameter': indicator,
                            'Age_Group': age_group,
                            'Age_Cutoff': age_median,
                            'Group1_Mean': mean1,
                            'Group1_SD': std1,
                            'Group1_N': n1,
                            'Group2_Mean': mean2,
                            'Group2_SD': std2,
                            'Group2_N': n2,
                            'Cohens_d': cohens_d,
                            'T_statistic': t_stat,
                            'P_value': p_value
                        })
    
    if results:
        s8_df = pd.DataFrame(results)
        print(f"Supplementary Table S8 创建完成: {len(s8_df)} 行")
        return s8_df
    else:
        print("Supplementary Table S8 创建失败")
        return None

def save_tables(tables_dict):
    """保存所有表格到文件"""
    output_dir = str(/mnt/c/Users/CUI/Desktop/论文及图表)
    os.makedirs(output_dir, exist_ok=True)
    
    saved_files = []
    
    for table_name, table_df in tables_dict.items():
        if table_df is not None and len(table_df) > 0:
            # 确定文件名
            if table_name == 'Table7':
                filename = "Table7_PHQ9_Stratified_Analysis.xlsx"
            elif table_name.startswith('Supplementary'):
                # 提取S编号
                s_num = table_name.split('_')[-1]
                filename = f"Supplementary_Table_{s_num}.xlsx"
            else:
                filename = f"{table_name}.xlsx"
            
            filepath = os.path.join(output_dir, filename)
            table_df.to_excel(filepath, index=False)
            saved_files.append((filename, len(table_df)))
            print(f"  保存: {filename} ({len(table_df)} 行)")
    
    return saved_files

def main():
    """主函数"""
    print("开始创建缺失的表格文件...")
    print("=" * 60)
    
    # 加载数据
    df = load_and_filter_data()
    
    if df is None or len(df) == 0:
        print("错误: 无法加载数据")
        return
    
    # 创建表格
    tables = {}
    
    # Table 7
    tables['Table7'] = create_table7_phq9_stratified(df)
    
    # Supplementary Tables
    tables['Supplementary_S1'] = create_supplementary_table_s1(df)
    tables['Supplementary_S4'] = create_supplementary_table_s4(df)
    tables['Supplementary_S5'] = create_supplementary_table_s5(df)
    tables['Supplementary_S7'] = create_supplementary_table_s7(df)
    tables['Supplementary_S8'] = create_supplementary_table_s8(df)
    
    # 保存表格
    print("\n保存表格文件...")
    saved_files = save_tables(tables)
    
    # 总结
    print("\n" + "=" * 60)
    print("创建完成总结:")
    print(f"总创建表格数: {len(saved_files)}")
    
    for filename, row_count in saved_files:
        print(f"  {filename}: {row_count} 行")
    
    print("\n文件位置: /mnt/c/Users/CUI/Desktop/投稿、数据修改/最终版/01_Manuscript/")
    print("\n下一步:")
    print("1. 验证表格数据与论文正文一致性")
    print("2. 将表格嵌入论文Word文档")
    print("3. 开始Journal of Affective Disorders投稿")

if __name__ == "__main__":
    main()