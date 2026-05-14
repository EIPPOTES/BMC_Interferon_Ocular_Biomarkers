"""
优化分析模块使用示例
展示基于多数据集验证的优化分析
"""

import pandas as pd
import numpy as np
import sys
import os

# 添加src目录到路径
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

try:
    from src.optimized_analysis import OptimizedAnalyzer
    print("成功导入优化分析模块")
except ImportError as e:
    print(f"导入失败: {e}")
    print("请确保在BMC_minimal目录下运行此示例")
    sys.exit(1)

def create_test_dataset():
    """创建测试数据集（模拟小样本眼科研究）"""
    
    np.random.seed(20260301)
    
    # 小样本设计：6对照 + 6病例
    n_samples = 12
    n_genes = 1500
    
    # 基础表达
    expression = np.random.lognormal(mean=5, sigma=0.7, size=(n_samples, n_genes))
    
    # 添加差异表达（模拟真实小样本情况）
    n_de = int(n_genes * 0.04)  # 4%的基因
    de_indices = np.random.choice(n_genes, n_de, replace=False)
    
    for idx in de_indices:
        fold_change = np.random.uniform(1.6, 2.5)
        expression[6:, idx] *= fold_change  # 后6个样本为病例组
    
    # 基因名
    gene_names = [f"GENE_{i:06d}" for i in range(n_genes)]
    
    # 样本名
    sample_names = [f"Control_{i}" for i in range(6)] + [f"Case_{i}" for i in range(6)]
    
    # 创建DataFrame
    expr_df = pd.DataFrame(expression, index=sample_names, columns=gene_names)
    
    # 样本信息
    sample_info = pd.DataFrame({
        'sample_id': sample_names,
        'group': ['control']*6 + ['case']*6,
        'batch': ['A']*4 + ['B']*4 + ['C']*4
    })
    
    return expr_df, sample_info

def main():
    """主函数：演示优化分析"""
    
    print("="*80)
    print("优化分析模块使用示例")
    print("基于多数据集验证的小样本眼科研究优化")
    print("="*80)
    
    # 1. 准备数据
    print("\n1. 准备测试数据...")
    expression_df, sample_df = create_test_dataset()
    
    print(f"   表达矩阵: {expression_df.shape}")
    print(f"   样本信息: {sample_df.shape}")
    print(f"   分组: {sample_df['group'].value_counts().to_dict()}")
    
    # 2. 运行优化分析
    print("\n2. 运行优化分析...")
    analyzer = OptimizedAnalyzer(output_dir="results/optimized_analysis_example")
    
    try:
        results = analyzer.analyze_with_optimization(
            expression_df, 
            sample_df, 
            condition_name="ophthalmology_small_sample_study"
        )
        
        # 3. 展示关键结果
        print("\n3. 分析结果摘要:")
        
        de = results.get("differential_expression", {})
        if de.get("success"):
            print(f"   差异表达分析:")
            print(f"     测试基因: {de['genes_tested']}")
            print(f"     显著基因: {de['significant_genes']} ({de['significance_rate']:.1%})")
            print(f"     中位效应大小: {de.get('median_effect_size', 0):.3f}")
        
        power = results.get("power_analysis", {})
        if power:
            print(f"   统计功效分析:")
            print(f"     当前功效: {power['estimated_current_power']:.1%}")
            print(f"     建议样本量: {power['required_samples_per_group_80_power']}每组")
            print(f"     充足性: {power['power_adequacy']}")
        
        trends = results.get("exploratory_trends", {})
        if trends.get("trend_genes_count", 0) > 0:
            print(f"   探索性趋势:")
            print(f"     有趋势基因: {trends['trend_genes_count']}")
        
        # 4. 生成使用指南
        print("\n4. 在你的研究中使用:")
        print("""
   步骤1: 导入模块
      from src.optimized_analysis import OptimizedAnalyzer
      
   步骤2: 准备数据
      expression_df = pd.read_csv('your_data.csv', index_col=0)
      sample_df = pd.read_csv('your_samples.csv')
      
   步骤3: 运行分析
      analyzer = OptimizedAnalyzer(output_dir='your_results')
      results = analyzer.analyze_with_optimization(
          expression_df, sample_df, 'your_study_name'
      )
      
   步骤4: 解释结果
      - 查看生成的优化报告
      - 关注统计功效分析
      - 考虑探索性趋势
      - 基于验证结果理解局限性
      
   步骤5: 论文报告
      - 描述基于验证的优化策略
      - 报告统计功效和样本量建议
      - 提供完整的优化分析报告
        """)
        
    except Exception as e:
        print(f"分析失败: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*80)
    print("优化分析的价值:")
    print("="*80)
    print("""
基于多数据集验证的优化分析提供:

1. **实证基础**
   - 基于4个眼科数据集的系统验证
   - 针对小样本特性的专门优化
   - 数据驱动的参数选择

2. **科学严谨性**
   - 自动适应样本量的统计方法
   - 详细的统计功效分析
   - 透明的优化决策过程

3. **临床实用性**
   - 针对眼科研究的专门优化
   - 提供具体的样本量建议
   - 连接生物信息学与临床实践

4. **方法学创新**
   - 验证驱动的优化框架
   - 小样本研究的专门工具
   - 可扩展的优化策略
    """)
    
    print("\n" + "="*80)
    print("查看生成的文件:")
    print("- results/optimized_analysis_example/optimized_analysis_report.json")
    print("- results/optimized_analysis_example/optimized_analysis_summary.md")
    print("="*80)

if __name__ == "__main__":
    main()
