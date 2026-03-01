"""
阴性结果分析使用示例
展示如何在小样本眼科研究中使用阴性结果分析模块
"""

import pandas as pd
import numpy as np
import sys
import os

# 添加src目录到路径
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

try:
    from src.negative_result_analysis import NegativeResultAnalyzer
    print("成功导入阴性结果分析模块")
except ImportError as e:
    print(f"导入失败: {e}")
    print("请确保在BMC_minimal目录下运行此示例")
    sys.exit(1)

def create_example_ophthalmology_data():
    """创建示例眼科数据（干眼症小样本研究）"""
    
    np.random.seed(2026)
    
    # 小样本设计：4个对照，4个干眼症患者
    n_samples = 8
    n_genes = 1500
    
    # 基础表达水平
    expression = np.random.lognormal(mean=5, sigma=1, size=(n_samples, n_genes))
    
    # 添加少量差异表达（模拟真实小样本研究）
    # 只有2%的基因有差异（小样本常见情况）
    n_de_genes = int(n_genes * 0.02)  # 30个基因
    de_indices = np.random.choice(n_genes, n_de_genes, replace=False)
    
    # 干眼症样本（后4个）
    case_indices = range(4, 8)
    for idx in de_indices:
        # 小到中等效应
        fold_change = np.random.uniform(1.3, 2.0)
        expression[case_indices, idx] *= fold_change
    
    # 基因名
    gene_names = [f"GENE_{i:06d}" for i in range(n_genes)]
    
    # 样本名
    sample_names = [f"Control_{i}" for i in range(4)] + [f"DryEye_{i}" for i in range(4)]
    
    # 创建DataFrame
    expr_df = pd.DataFrame(expression, index=sample_names, columns=gene_names)
    
    # 样本信息
    sample_info = pd.DataFrame({
        'sample_id': sample_names,
        'group': ['control']*4 + ['dry_eye']*4,
        'age': np.random.randint(40, 70, n_samples),
        'sex': np.random.choice(['M', 'F'], n_samples),
        'disease_duration_years': [0]*4 + np.random.randint(1, 10, 4).tolist()
    })
    
    return expr_df, sample_info

def main():
    """主函数：演示阴性结果分析"""
    
    print("="*80)
    print("阴性结果分析模块使用示例")
    print("小样本眼科研究：干眼症转录组分析")
    print("="*80)
    
    # 1. 准备数据
    print("\n1. 准备示例数据...")
    expression_df, sample_df = create_example_ophthalmology_data()
    
    print(f"   表达矩阵: {expression_df.shape}")
    print(f"   样本信息: {sample_df.shape}")
    print(f"   分组: {sample_df['group'].value_counts().to_dict()}")
    
    # 2. 运行阴性结果分析
    print("\n2. 运行阴性结果分析...")
    analyzer = NegativeResultAnalyzer(output_dir="results/dry_eye_negative_analysis")
    
    try:
        results = analyzer.analyze(
            expression_df, 
            sample_df, 
            condition_name="dry_eye_small_study"
        )
        
        # 3. 展示关键结果
        print("\n3. 分析结果摘要:")
        de = results.get('differential_expression', {})
        if de:
            print(f"   测试基因数: {de.get('genes_tested', 'N/A')}")
            print(f"   显著基因数: {de.get('significant_genes', 0)}")
        
        power = results.get('statistical_power', {})
        if power.get('median_effect_size') is not None:
            print(f"   中位效应大小: {power['median_effect_size']:.3f}")
            if power.get('estimated_power') is not None:
                print(f"   估计统计功效: {power['estimated_power']:.1%}")
        
        trends = results.get('exploratory_trends', {})
        if trends.get('trend_genes_count', 0) > 0:
            print(f"   有趋势基因数: {trends['trend_genes_count']}")
        
        # 4. 生成使用指南
        print("\n4. 使用指南:")
        print("""
   在你的研究中使用阴性结果分析模块:
   
   步骤1: 准备数据
      expression_df = pd.read_csv('your_expression_data.csv', index_col=0)
      sample_df = pd.read_csv('your_sample_info.csv')
   
   步骤2: 运行分析
      from src.negative_result_analysis import NegativeResultAnalyzer
      analyzer = NegativeResultAnalyzer(output_dir='your_results')
      results = analyzer.analyze(expression_df, sample_df, 'your_study_name')
   
   步骤3: 解释结果
      - 查看生成的报告文件
      - 关注统计功效分析
      - 考虑探索性趋势基因
      - 基于结果设计后续研究
   
   步骤4: 论文报告
      - 在方法部分描述分析流程
      - 在结果部分报告阴性发现
      - 在讨论部分解释统计局限性
      - 提供完整的分析报告作为补充材料
        """)
        
    except Exception as e:
        print(f"分析失败: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*80)
    print("模块价值总结:")
    print("="*80)
    print("""
即使没有显著基因，阴性结果分析模块提供:

1. **科学严谨性保障**
   - 完整的数据处理记录
   - 适当的统计方法
   - 多重检验校正

2. **研究设计指导**
   - 统计功效评估
   - 样本量建议
   - 效应大小估计

3. **结果解释框架**
   - 阴性结果的科学解释
   - 探索性趋势识别
   - 临床意义分析

4. **可重现性标准**
   - 完整的分析流程
   - 透明的报告生成
   - 标准化的输出格式
    """)
    
    print("\n" + "="*80)
    print("运行示例完成!")
    print("查看生成的文件:")
    print("- results/dry_eye_negative_analysis/comprehensive_negative_analysis_report.json")
    print("- results/dry_eye_negative_analysis/executive_summary.md")
    print("="*80)

if __name__ == "__main__":
    main()
