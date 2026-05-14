"""
整合阴性结果分析模块到BMC Bioinformatics代码包
"""

import os
from pathlib import Path

def create_negative_analysis_module():
    """创建阴性结果分析模块"""
    
    module_content = '''"""
阴性结果分析模块
专门处理和分析小样本研究中的阴性结果
BMC Bioinformatics代码包增强组件
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
    
    主要功能:
    1. 差异表达分析（小样本优化）
    2. 统计功效评估
    3. 探索性趋势识别
    4. 样本异质性分析
    5. 综合报告生成
    """
    
    def __init__(self, output_dir="results/negative_analysis"):
        """
        初始化分析器
        
        Parameters:
        -----------
        output_dir : str or Path
            输出目录路径
        """
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
            
        Returns:
        --------
        dict : 包含所有分析结果的字典
        """
        
        print(f"\\n{'='*60}")
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
            "total_genes": expression_df.shape[1],
            "analysis_date": pd.Timestamp.now().isoformat()
        }
        
        # 1. 基本差异分析
        print("\\n1. 执行差异表达分析...")
        de_results = self._perform_differential_analysis(
            expression_df, idx1, idx2, group1, group2
        )
        results.update(de_results)
        
        # 2. 统计功效分析
        print("2. 统计功效分析...")
        power_results = self._power_analysis(
            expression_df, idx1, idx2, n1, n2
        )
        results.update(power_results)
        
        # 3. 探索性趋势分析
        print("3. 探索性趋势分析...")
        trend_results = self._trend_analysis(
            expression_df, idx1, idx2, group1, group2
        )
        results.update(trend_results)
        
        # 4. 样本异质性分析
        print("4. 样本异质性分析...")
        heterogeneity_results = self._heterogeneity_analysis(
            expression_df, sample_df
        )
        results.update(heterogeneity_results)
        
        # 5. 生成报告
        print("5. 生成综合报告...")
        self._generate_reports(results, condition_name)
        
        print(f"\\n分析完成! 报告保存在: {self.output_dir}")
        return results
    
    def _validate_inputs(self, expression_df, sample_df):
        """验证输入数据"""
        if 'group' not in sample_df.columns:
            raise ValueError("sample_df必须包含'group'列")
        
        if expression_df.shape[0] != len(sample_df):
            raise ValueError(f"表达矩阵样本数({expression_df.shape[0]})与样本信息数({len(sample_df)})不匹配")
    
    def _perform_differential_analysis(self, expression_df, idx1, idx2, group1, group2):
        """执行差异表达分析（小样本优化）"""
        
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
            
            # Welch's t检验（方差不齐）
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
            
            # 多重检验校正 (Bonferroni - 保守，适合小样本)
            results_df['adj_pvalue'] = results_df['pvalue'] * len(results_df)
            results_df['adj_pvalue'] = results_df['adj_pvalue'].clip(upper=1.0)
            
            # 显著性判断（宽松阈值适应小样本）
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
            
            # 简单功效估计
            if median_effect > 0:
                # 使用t分布近似计算功效
                from scipy import stats as scipy_stats
                
                # 非中心化参数
                ncp = median_effect * np.sqrt((n1 * n2) / (n1 + n2))
                
                # 临界t值（双尾检验，α=0.05）
                df = n1 + n2 - 2
                critical_t = scipy_stats.t.ppf(0.975, df)
                
                # 计算功效
                power = 1 - scipy_stats.t.cdf(critical_t - ncp, df)
                power = max(0, min(1, power))  # 限制在0-1之间
                
                # 建议样本量（经验公式）
                # 对于中等效应(0.5)，80%功效需要每组约26个样本
                if median_effect > 0:
                    recommended_n = int(26 * (0.5 / median_effect)**2)
                    recommended_n = max(recommended_n, 10)  # 至少10个样本
                else:
                    recommended_n = 26
                
                adequacy = "adequate" if power >= 0.8 else "inadequate"
                
                return {
                    "statistical_power": {
                        "median_effect_size": float(median_effect),
                        "estimated_power": float(power),
                        "current_sample_size_adequacy": adequacy,
                        "recommended_sample_size_per_group": recommended_n,
                        "interpretation": f"Current study has {power:.1%} power to detect median effect size ({median_effect:.2f})"
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
                    trend_strength = -np.log10(pval) * abs(np.log2(fc))  # 综合评分
                    
                    trend_genes.append({
                        'gene': gene,
                        'pvalue': pval,
                        'fold_change': fc,
                        'log2fc': np.log2(fc),
                        'group1_mean': vals1.mean(),
                        'group2_mean': vals2.mean(),
                        'trend_strength': trend_strength
                    })
        
        if trend_genes:
            trend_df = pd.DataFrame(trend_genes)
            trend_df = trend_df.sort_values('trend_strength', ascending=False)
            
            trend_df.to_csv(self.output_dir / "trend_genes.csv", index=False)
            
            # 取前10个最有趋势的基因
            top_trends = trend_df.head(10)[['gene', 'pvalue', 'fold_change', 'trend_strength']].to_dict('records')
            
            return {
                "exploratory_trends": {
                    "trend_genes_count": len(trend_genes),
                    "top_trend_genes": top_trends,
                    "trend_genes_file": str(self.output_dir / "trend_genes.csv")
                }
            }
        
        return {"exploratory_trends": {"trend_genes_count": 0}}
    
    def _heterogeneity_analysis(self, expression_df, sample_df):
        """样本异质性分析"""
        
        try:
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
                
                within_median = np.median(within_group_corr)
                between_median = np.median(between_group_corr)
                
                if between_median != 0:
                    heterogeneity_ratio = within_median / between_median
                else:
                    heterogeneity_ratio = np.nan
                
                return {
                    "sample_heterogeneity": {
                        "within_group_corr_median": float(within_median),
                        "between_group_corr_median": float(between_median),
                        "heterogeneity_ratio": float(heterogeneity_ratio),
                        "visualization": str(self.output_dir / "heterogeneity_analysis.png"),
                        "interpretation": "High within-group heterogeneity may reduce statistical power for detecting between-group differences"
                    }
                }
        except Exception as e:
            print(f"异质性分析错误: {e}")
        
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
- **分析日期**: {results.get('analysis_date', 'N/A')}

## 主要发现

### 差异表达分析
"""
        
        de = results.get('differential_expression', {})
        if de:
            summary += f"- **测试基因数**: {de.get('genes_tested', 'N/A')}\\n"
            summary += f"- **显著基因数**: {de.get('significant_genes', 0)} (p<0.05, FDR校正)\\n"
            summary += f"- **显著性比例**: {de.get('significance_rate', 0):.1%}\\n"
            if de.get('min_pvalue'):
                summary += f"- **最小p值**: {de['min_pvalue']:.3e}\\n"
            if de.get('max_effect_size'):
                summary += f"- **最大效应大小**: {de['max_effect_size']:.3f} (Cohen's d)\\n"
        
        summary += "\\n### 统计功效分析\\n"
        
        power = results.get('statistical_power', {})
        if power.get('median_effect_size') is not None:
            summary += f"- **中位效应大小**: {power['median_effect_size']:.3f} (Cohen's d)\\n"
            if power.get('estimated_power') is not None:
                summary += f"- **估计统计功效**: {power['estimated_power']:.1%}\\n"
            summary += f"- **样本量充足性**: {power.get('current_sample_size_adequacy', 'N/A')}\\n"
            if power.get('recommended_sample_size_per_group'):
                summary += f"- **建议每组样本量**: {power['recommended_sample_size_per_group']}\\n"
        
        summary += "\\n### 探索性趋势\\n"
        
        trends = results.get('exploratory_trends', {})
        if trends.get('trend_genes_count', 0) > 0:
            summary += f"- **有趋势基因数**: {trends['trend_genes_count']} (0.05<p<0.2, FC>1.5)\\n"
            summary += "- **最有趋势的基因**:\\n"
            for gene in trends.get('top_trend_genes', [])[:5