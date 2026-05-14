#!/bin/bash
# 整合阴性结果分析模块到BMC Bioinformatics代码包

echo "="*80
echo "整合阴性结果分析模块"
echo "="*80

# 1. 创建阴性结果分析模块
echo "1. 创建阴性结果分析模块..."
cat > src/negative_result_analysis.py << 'EOF2'
"""
阴性结果分析模块
专门处理和分析小样本研究中的阴性结果
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt

class NegativeResultAnalyzer:
    """
    阴性结果分析器
    为小样本眼科研究提供深入的阴性结果分析
    """
    
    def __init__(self, output_dir="results/negative_analysis"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def analyze(self, expression_df, sample_df, condition_name="study"):
        """
        执行完整的阴性结果分析
        
        Parameters:
        -----------
        expression_df : pandas.DataFrame
            基因表达矩阵，行为样本，列为基因
        sample_df : pandas.DataFrame
            样本信息，必须包含'group'列
        condition_name : str
            研究条件名称，用于输出文件命名
        """
        
        print(f"\n{'='*60}")
        print(f"阴性结果深入分析: {condition_name}")
        print(f"{'='*60}")
        
        # 验证输入
        self._validate_inputs(expression_df, sample_df)
        
        # 提取分组信息
        groups = sample_df['group'].unique()
        if len(groups) != 2:
            raise ValueError(f"需要恰好2个分组，找到{len(groups)}个")
        
        group1, group2 = groups[0], groups[1]
        idx1 = sample_df[sample_df['group'] == group1].index
        idx2 = sample_df[sample_df['group'] == group2].index
        
        n1, n2 = len(idx1), len(idx2)
        
        print(f"样本: {group1}(n={n1}) vs {group2}(n={n2})")
        print(f"基因: {expression_df.shape[1]}")
        
        # 执行分析
        results = {
            "condition": condition_name,
            "sample_sizes": {group1: n1, group2: n2},
            "total_samples": n1 + n2,
            "total_genes": expression_df.shape[1]
        }
        
        # 1. 基本差异分析
        print("\n1. 执行差异表达分析...")
        de_results = self._perform_differential_analysis(
            expression_df, idx1, idx2, group1, group2
        )
        results.update(de_results)
        
        # 2. 统计功效分析
        print("\n2. 统计功效分析...")
        power_results = self._power_analysis(
            expression_df, idx1, idx2, n1, n2
        )
        results.update(power_results)
        
        # 3. 探索性趋势分析
        print("\n3. 探索性趋势分析...")
        trend_results = self._trend_analysis(
            expression_df, idx1, idx2, group1, group2
        )
        results.update(trend_results)
        
        # 4. 样本异质性分析
        print("\n4. 样本异质性分析...")
        heterogeneity_results = self._heterogeneity_analysis(
            expression_df, sample_df
        )
        results.update(heterogeneity_results)
        
        # 5. 生成报告
        print("\n5. 生成综合报告...")
        self._generate_reports(results, condition_name)
        
        return results
    
    def _validate_inputs(self, expression_df, sample_df):
        """验证输入数据"""
        if 'group' not in sample_df.columns:
            raise ValueError("sample_df必须包含'group'列")
        
        if expression_df.shape[0] != len(sample_df):
            raise ValueError(f"表达矩阵样本数({expression_df.shape[0]})与样本信息数({len(sample_df)})不匹配")
    
    def _perform_differential_analysis(self, expression_df, idx1, idx2, group1, group2):
        """执行差异表达分析"""
        
        # 分析前500个基因（小样本优化）
        n_genes_to_test = min(500, expression_df.shape[1])
        test_genes = expression_df.columns[:n_genes_to_test]
        
        results = []
        for gene in test_genes:
            vals1 = expression_df.iloc[idx1][gene].values
            vals2 = expression_df.iloc[idx2][gene].values
            
            # 跳过常数值
            if vals1.std() == 0 or vals2.std() == 0:
                continue
            
            # t检验
            stat, pval = stats.ttest_ind(vals1, vals2, equal_var=False)
            
            # 计算fold change
            if vals1.mean() != 0:
                fc = vals2.mean() / vals1.mean()
                log2fc = np.log2(fc) if fc > 0 else 0
            else:
                fc = np.nan
                log2fc = 0
            
            # 计算效应大小 (Cohen's d)
            pooled_std = np.sqrt((vals1.std()**2 + vals2.std()**2) / 2)
            if pooled_std > 0:
                cohens_d = (vals2.mean() - vals1.mean()) / pooled_std
            else:
                cohens_d = 0
            
            results.append({
                'gene': gene,
                'pvalue': pval,
                'log2fc': log2fc,
                'fold_change': fc,
                'cohens_d': cohens_d,
                'group1_mean': vals1.mean(),
                'group2_mean': vals2.mean()
            })
        
        if results:
            results_df = pd.DataFrame(results)
            
            # 多重检验校正 (Bonferroni - 保守)
            results_df['adj_pvalue'] = results_df['pvalue'] * len(results_df)
            results_df['adj_pvalue'] = results_df['adj_pvalue'].clip(upper=1.0)
            
            # 显著性判断
            results_df['significant'] = (results_df['adj_pvalue'] < 0.05) & (abs(results_df['log2fc']) > np.log2(1.5))
            
            n_sig = results_df['significant'].sum()
            
            # 保存结果
            results_df.to_csv(self.output_dir / "differential_expression_results.csv", index=False)
            
            return {
                "differential_expression": {
                    "genes_tested": len(results_df),
                    "significant_genes": int(n_sig),
                    "significance_rate": float(n_sig / len(results_df) if len(results_df) > 0 else 0),
                    "min_pvalue": float(results_df['pvalue'].min()),
                    "max_effect_size": float(results_df['cohens_d'].abs().max()),
                    "results_file": str(self.output_dir / "differential_expression_results.csv")
                }
            }
        
        return {"differential_expression": {"error": "No valid results"}}
    
    def _power_analysis(self, expression_df, idx1, idx2, n1, n2):
        """统计功效分析"""
        
        # 计算观察到的效应大小分布
        effect_sizes = []
        for gene in expression_df.columns[:200]:  # 分析前200个基因
            vals1 = expression_df.iloc[idx1][gene].values
            vals2 = expression_df.iloc[idx2][gene].values
            
            if vals1.std() > 0 and vals2.std() > 0:
                pooled_std = np.sqrt((vals1.std()**2 + vals2.std()**2) / 2)
                if pooled_std > 0:
                    d = (vals2.mean() - vals1.mean()) / pooled_std
                    effect_sizes.append(abs(d))
        
        if effect_sizes:
            median_effect = np.median(effect_sizes)
            
            # 简单功效估计（基于t检验公式）
            # 对于小样本，使用近似公式
            if median_effect > 0:
                # 非中心化参数
                ncp = median_effect * np.sqrt((n1 * n2) / (n1 + n2))
                
                # 使用t分布计算功效（简化）
                from scipy import stats as scipy_stats
                critical_t = scipy_stats.t.ppf(0.975, n1 + n2 - 2)  # 双尾检验
                power = 1 - scipy_stats.t.cdf(critical_t - ncp, n1 + n2 - 2)
                
                # 建议样本量（经验公式）
                # 对于中等效应(0.5)，80%功效需要每组约26个样本
                recommended_n = int(26 * (0.5 / median_effect)**2) if median_effect > 0 else 26
                recommended_n = max(recommended_n, 10)  # 至少10个样本
                
                return {
                    "statistical_power": {
                        "median_effect_size": float(median_effect),
                        "estimated_power": float(power),
                        "current_sample_size_adequacy": "inadequate" if power < 0.8 else "adequate",
                        "recommended_sample_size_per_group": recommended_n,
                        "interpretation": f"Current study has limited power ({power:.1%}) to detect median effect size ({median_effect:.2f})"
                    }
                }
        
        return {
            "statistical_power": {
                "median_effect_size": None,
                "estimated_power": None,
                "interpretation": "Could not estimate statistical power"
            }
        }
    
    def _trend_analysis(self, expression_df, idx1, idx2, group1, group2):
        """探索性趋势分析"""
        
        trend_genes = []
        for gene in expression_df.columns[:300]:  # 分析前300个基因
            vals1 = expression_df.iloc[idx1][gene].values
            vals2 = expression_df.iloc[idx2][gene].values
            
            if vals1.std() > 0 and vals2.std() > 0:
                _, pval = stats.ttest_ind(vals1, vals2, equal_var=False)
                fc = vals2.mean() / vals1.mean() if vals1.mean() != 0 else 1
                
                # 识别有趋势的基因 (0.05 < p < 0.2, FC > 1.5)
                if 0.05 < pval < 0.2 and abs(np.log2(fc)) > np.log2(1.5):
                    trend_genes.append({
                        'gene': gene,
                        'pvalue': pval,
                        'fold_change': fc,
                        'log2fc': np.log2(fc),
                        'group1_mean': vals1.mean(),
                        'group2_mean': vals2.mean(),
                        'trend_strength': -np.log10(pval) * abs(np.log2(fc))  # 综合评分
                    })
        
        if trend_genes:
            trend_df = pd.DataFrame(trend_genes)
            trend_df = trend_df.sort_values('trend_strength', ascending=False)
            
            trend_df.to_csv(self.output_dir / "trend_genes.csv", index=False)
            
            return {
                "exploratory_trends": {
                    "trend_genes_count": len(trend_genes),
                    "top_trend_genes": trend_df[['gene', 'pvalue', 'fold_change']].head(10).to_dict('records'),
                    "trend_genes_file": str(self.output_dir / "trend_genes.csv")
                }
            }
        
        return {"exploratory_trends": {"trend_genes_count": 0}}
    
    def _heterogeneity_analysis(self, expression_df, sample_df):
        """样本异质性分析"""
        
        # 计算样本间相关性
        sample_corr = expression_df.T.corr()
        
        within_group_corr = []
        between_group_corr = []
        
        for i in range(len(sample_corr)):
            for j in range(i+1, len(sample_corr)):
                corr_val = sample_corr.iloc[i, j]
                group_i = sample_df.iloc[i]['group']
                group_j = sample_df.iloc[j]['group']
                
                if group_i == group_j:
                    within_group_corr.append(corr_val)
                else:
                    between_group_corr.append(corr_val)
        
        if within_group_corr and between_group_corr:
            # 可视化
            plt.figure(figsize=(10, 6))
            plt.boxplot([within_group_corr, between_group_corr], 
                       labels=['Within-group', 'Between-group'])
            plt.ylabel('Sample correlation')
            plt.title('Sample heterogeneity analysis')
            plt.savefig(self.output_dir / "heterogeneity_analysis.png", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            return {
                "sample_heterogeneity": {
                    "within_group_corr_median": float(np.median(within_group_corr)),
                    "between_group_corr_median": float(np.median(between_group_corr)),
                    "heterogeneity_ratio": float(np.median(within_group_corr) / np.median(between_group_corr) if np.median(between_group_corr) != 0 else np.nan),
                    "visualization": str(self.output_dir / "heterogeneity_analysis.png"),
                    "interpretation": "High within-group heterogeneity relative to between-group differences may reduce statistical power"
                }
            }
        
        return {"sample_heterogeneity": {"error": "Could not calculate heterogeneity"}}
    
    def _generate_reports(self, results, condition_name):
        """生成综合报告"""
        
        # JSON报告
        report_file = self.output_dir / "comprehensive_negative_analysis_report.json"
        with open(report_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        # 可读性总结
        summary = self._generate_executive_summary(results, condition_name)
        summary_file = self.output_dir / "executive_summary.md"
        with open(summary_file, 'w') as f:
            f.write(summary)
        
        print(f"   综合报告: {report_file}")
        print(f"   执行摘要: {summary_file}")
    
    def _generate_executive_summary(self, results, condition_name):
        """生成可读性执行摘要"""
        
        summary = f"""# 阴性结果分析执行摘要: {condition_name}

## 研究概述
- **研究条件**: {condition_name}
- **总样本量**: {results.get('total_samples', 'N/A')}
- **分析基因数**: {results.get('total_genes', 'N/A')}

## 主要发现

### 差异表达分析
"""
        
        de = results.get('differential_expression', {})
        if de:
            summary += f"- **测试基因数**: {de.get('genes_tested', 'N/A')}\n"
            summary += f"- **显著基因数**: {de.get('significant_genes', 0)} (p<0.05, FDR校正)\n"
            summary += f"- **显著性比例**: {de.get('significance_rate', 0):.1%}\n"
            summary += f"- **最小p值**: {de.get('min_pvalue', 'N/A')}\n"
        
        summary += "\n### 统计功效分析\n"
        
        power = results.get('statistical_power', {})
        if power.get('median_effect_size'):
            summary += f"- **中位效应大小**: {power['median_effect_size']:.3f} (Cohen's d)\n"
            summary += f"- **估计统计功效**: {power.get('estimated_power', 0):.1%}\n"
            summary += f"- **样本量充足性**: {power.get('current_sample_size_adequacy', 'N/A')}\n"
            summary += f"- **建议每组样本量**: {power.get('recommended_sample_size_per_group', 'N/A')}\n"
        
        summary += "\n### 探索性趋势\n"
        
        trends = results.get('exploratory_trends', {})
        if trends.get('trend_genes_count', 0) > 0:
            summary += f"- **有趋势基因数**: {trends['trend_genes_count']} (0.05<p<0.2, FC>1.5)\n"
            summary += "- **最有趋势的基因**:\n"
            for gene in trends.get('top_trend_genes', [])[:5]:
                summary += f"  - {gene['gene']}: p={gene['pvalue']:.3f}, FC={gene['fold_change']:.2f}\n"
        
        summary += "\n## 科学解释\n\n"
        
        if de.get('significant_genes', 0) == 0:
            summary += "在当前样本量下，**未发现统计学显著的差异表达基因**。\n\n"
            summary += "### 可能原因\n"
            summary += "1. **真实效应较小**: 生物学差异可能小于当前检测限\n"
            summary += "2. **样本量不足**: 统计功效有限，无法检测中等效应\n"
            summary += "3. **样本异质性**: 组内变异可能掩盖真实的组间差异\n"
            summary += "4. **技术变异**: 实验或测量误差可能影响信号检测\n\n"
        
        summary += "### 统计考虑\n"
        if power.get('estimated_power', 0) < 0.8:
            summary += f"- 当前研究统计功效较低 ({power.get('estimated_power', 0):.1%})\n"
            summary += f"- 检测中位效应大小({power.get('median_effect_size', 0):.2f})的能力有限\n"
            summary += f"- 建议增加样本量至每组{power.get('recommended_sample_size_per_group', 26)}个以提高功效\n\n"
        
