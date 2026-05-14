"""
眼科GEO数据库搜索和使用指南
基于用户提供的专业指导
"""

import json
from pathlib import Path

def create_geo_search_guide():
    """创建GEO搜索指南"""
    
    guide = {
        "title": "眼科GEO数据库全流程使用指南",
        "author": "Based on user-provided professional guidance",
        "date": "2026-03-01",
        "purpose": "帮助眼科研究者高效使用NCBI GEO数据库",
        
        "sections": [
            {
                "title": "GEO核心编号系统",
                "content": [
                    "GSE (Series - 系列/项目): 最重要！整个实验项目的总称",
                    "GSM (Sample - 样本): 单个样本的数据",
                    "GPL (Platform - 平台): 测序或芯片的技术平台",
                    "GDS (DataSet - 数据集): 官方整理的精选数据集"
                ]
            },
            {
                "title": "高级检索策略",
                "content": [
                    "布尔逻辑: AND, OR, NOT",
                    "字段限定: [Title], [Title/Abstract], [Organism], [Entry Type]",
                    "眼科检索公式示例:",
                    "  '\"Glaucoma\"[Title/Abstract] AND \"Trabecular meshwork\"[All Fields] AND \"Homo sapiens\"[Organism] AND \"gse\"[Entry Type]'"
                ]
            },
            {
                "title": "GSE页面关键信息",
                "content": [
                    "Title & Summary: 快速判断研究内容",
                    "Overall design: 实验分组情况",
                    "Contributor(s): 数据上传作者，可找原始论文",
                    "Platforms: 技术平台信息",
                    "Samples: 所有GSM样本列表"
                ]
            },
            {
                "title": "数据下载指南",
                "content": [
                    "Series Matrix File(s): 最推荐！已处理的表达矩阵",
                    "SOFT formatted family file(s): 元数据，较少使用",
                    "Supplementary file(s): 原始数据或计数矩阵"
                ]
            },
            {
                "title": "GEO2R使用指南",
                "content": [
                    "Define groups: 创建分组（如Disease和Control）",
                    "Assign samples: 分配样本到对应组",
                    "Analyze: 执行差异分析",
                    "查看结果: 关注P.Value和logFC列"
                ]
            },
            {
                "title": "眼科经典数据集推荐",
                "content": [
                    "AMD: GSE29801, GSE115828",
                    "Glaucoma: GSE27276, GSE138114",
                    "Diabetic Retinopathy: GSE60436",
                    "Keratoconus: GSE112155, GSE15020",
                    "Dry Eye: GSE43671"
                ]
            },
            {
                "title": "实用技巧",
                "content": [
                    "使用具体疾病英文全称搜索，比\"Ophthalmology\"更精准",
                    "务必加上\"Homo sapiens\"[Organism]限定人类数据",
                    "查看原始论文了解实验设计和样本信息",
                    "下载Series Matrix File进行初步分析"
                ]
            }
        ],
        
        "example_searches": [
            {
                "disease": "年龄相关性黄斑变性 (AMD)",
                "query": "\"Age-related Macular Degeneration\"[Title] OR \"AMD\"[Title] AND \"Homo sapiens\"[Organism] AND \"gse\"[Entry Type]",
                "expected_results": "GSE29801, GSE115828等"
            },
            {
                "disease": "青光眼",
                "query": "\"Glaucoma\"[Title/Abstract] AND \"Trabecular meshwork\"[All Fields] AND \"Homo sapiens\"[Organism]",
                "expected_results": "GSE27276等"
            },
            {
                "disease": "糖尿病视网膜病变",
                "query": "\"Diabetic retinopathy\"[Title] AND \"fibrovascular membrane\"[All Fields] AND \"Homo sapiens\"[Organism]",
                "expected_results": "GSE60436等"
            }
        ],
        
        "integration_with_framework": {
            "purpose": "将GEO数据与本分析框架结合使用",
            "steps": [
                "1. 使用本指南搜索合适的眼科GEO数据集",
                "2. 下载Series Matrix File或计数矩阵",
                "3. 使用本框架的`process_real_data.py`处理数据",
                "4. 运行`analysis_pipeline.py`进行分析",
                "5. 使用`optimized_analysis.py`进行优化分析",
                "6. 生成完整的分析报告"
            ],
            "benefits": [
                "标准化的工作流程",
                "优化的分析方法",
                "透明的结果报告",
                "可重现的分析过程"
            ]
        }
    }
    
    return guide

def create_search_scripts():
    """创建搜索脚本示例"""
    
    scripts = {
        "python_search_example": """
# Python搜索GEO的示例代码
import requests

def search_geo(query):
    \"\"\"搜索GEO数据库\"\"\"
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "gds",
        "term": query,
        "retmode": "json",
        "retmax": 50
    }
    
    response = requests.get(base_url, params=params)
    return response.json()

# 搜索AMD相关数据集
amd_query = '"Age-related Macular Degeneration"[Title] AND "Homo sapiens"[Organism]'
results = search_geo(amd_query)
print(f"找到{results.get('esearchresult', {}).get('count', 0)}个AMD数据集")
""",
        
        "r_analysis_example": """
# R分析GEO数据的示例代码
library(GEOquery)
library(limma)

# 下载GSE数据
gse <- getGEO("GSE29801", GSEMatrix = TRUE)

# 获取表达矩阵和样本信息
exprs <- exprs(gse[[1]])
pdata <- pData(gse[[1]])

# 差异表达分析
design <- model.matrix(~0 + factor(pdata$group))
colnames(design) <- c("control", "case")
fit <- lmFit(exprs, design)
contrast.matrix <- makeContrasts(case-control, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取结果
results <- topTable(fit2, number=1000, adjust.method="fdr")
""",
        
        "bash_download_example": """
# Bash下载GEO数据的示例
# 下载Series Matrix File
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29801/matrix/GSE29801_series_matrix.txt.gz

# 解压
gunzip GSE29801_series_matrix.txt.gz

# 查看前几行
head -n 20 GSE29801_series_matrix.txt
"""
    }
    
    return scripts

def main():
    """主函数"""
    
    print("="*80)
    print("眼科GEO数据库使用指南生成器")
    print("基于用户提供的专业指导")
    print("="*80)
    
    # 创建指南
    guide = create_geo_search_guide()
    
    # 保存指南
    output_dir = Path("geo_search_guides")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # JSON格式
    with open(output_dir / "geo_search_guide.json", 'w') as f:
        json.dump(guide, f, indent=2)
    
    # Markdown格式
    with open(output_dir / "geo_search_guide.md", 'w') as f:
        f.write(f"# {guide['title']}\n\n")
        f.write(f"**作者**: {guide['author']}\n")
        f.write(f"**日期**: {guide['date']}\n")
        f.write(f"**目的**: {guide['purpose']}\n\n")
        
        for section in guide['sections']:
            f.write(f"## {section['title']}\n\n")
            for item in section['content']:
                f.write(f"- {item}\n")
            f.write("\n")
        
        f.write("## 搜索示例\n\n")
        for example in guide['example_searches']:
            f.write(f"### {example['disease']}\n")
            f.write(f"- **搜索式**: `{example['query']}`\n")
            f.write(f"- **预期结果**: {example['expected_results']}\n\n")
        
        f.write("## 与本框架集成\n\n")
        f.write(f"**目的**: {guide['integration_with_framework']['purpose']}\n\n")
        f.write("**步骤**:\n")
        for i, step in enumerate(guide['integration_with_framework']['steps'], 1):
            f.write(f"{i}. {step}\n")
        
        f.write("\n**优势**:\n")
        for benefit in guide['integration_with_framework']['benefits']:
            f.write(f"- {benefit}\n")
    
    # 创建搜索脚本
    scripts = create_search_scripts()
    with open(output_dir / "search_scripts.py", 'w') as f:
        f.write(scripts['python_search_example'])
    
    with open(output_dir / "search_scripts.R", 'w') as f:
        f.write(scripts['r_analysis_example'])
    
    with open(output_dir / "search_scripts.sh", 'w') as f:
        f.write(scripts['bash_download_example'])
    
    print(f"\n指南生成完成!")
    print(f"文件保存至: {output_dir}/")
    print(f"  - geo_search_guide.json (JSON格式)")
    print(f"  - geo_search_guide.md (Markdown格式)")
    print(f"  - search_scripts.py (Python示例)")
    print(f"  - search_scripts.R (R示例)")
    print(f"  - search_scripts.sh (Bash示例)")
    
    print(f"\n{'='*80}")
    print("使用建议:")
    print("1. 使用本指南搜索合适的眼科GEO数据集")
    print("2. 下载数据后使用本框架进行分析")
    print("3. 结合用户提供的专业数据集信息")
    print("4. 生成完整的分析报告")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()