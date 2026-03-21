#!/usr/bin/env python3

def main():
    """主函数，包装原有执行代码"""
    """
    生成最新版统计分析报告
    """

    import os
    import pandas as pd
    from datetime import datetime
    import re

    print("=" * 80)
    print("最新版统计分析结果报告")
    print("生成时间: 2026-03-15 15:23")
    print("=" * 80)

    # 定义分析类型
    analysis_types = [
        ("1. 描述性统计分析", ["描述性", "Descriptive"]),
        ("2. 组间比较（MDD vs Control）", ["组间比较", "Group_Comparison", "Group Comparison"]),
        ("3. 相关性分析", ["相关性", "Correlation"]),
        ("4. ROC分析", ["ROC", "AUC"]),
        ("5. 多变量回归分析", ["多变量回归", "Multivariate", "回归分析"]),
        ("6. 亚组分析", ["亚组分析", "Subgroup"])
    ]

    # 检查文件并获取信息
    def check_file_info(filepath, analysis_name):
        info = {
            "文件": os.path.basename(filepath),
            "路径": filepath,
            "修改时间": datetime.fromtimestamp(os.path.getmtime(filepath)).strftime('%Y-%m-%d %H:%M'),
            "大小_KB": round(os.path.getsize(filepath) / 1024, 1),
            "分析类型": analysis_name,
            "样本信息": "未知",
            "数据版本": "未知",
            "状态": "待验证"
        }

        try:
            # 尝试读取Excel文件
            xl = pd.ExcelFile(filepath)
            sheets = xl.sheet_names
            info["工作表"] = sheets

            # 读取第一个工作表
            if sheets:
                df = xl.parse(sheets[0])
                info["行数"] = df.shape[0]
                info["列数"] = df.shape[1]

                # 尝试推断样本信息
                sample_keywords = ["样本量", "样本", "n", "N", "患者", "眼", "eyes", "participants"]
                for col in df.columns:
                    if any(keyword in str(col) for keyword in sample_keywords):
                        if df.shape[0] > 0:
                            sample_val = df[col].iloc[0]
                            if pd.notna(sample_val):
                                info["样本信息"] = f"{col}: {sample_val}"
                                break

                # 检查是否包含499或485样本信息
                if "样本信息" in info and info["样本信息"] != "未知":
                    if "499" in str(info["样本信息"]):
                        info["数据版本"] = "499眼版本"
                    elif "485" in str(info["样本信息"]):
                        info["数据版本"] = "485眼版本"
                    elif "463" in str(info["样本信息"]):
                        info["数据版本"] = "463眼(499眼版本子集)"

                # 特殊检查：文件名中包含版本信息
                filename = info["文件"]
                if "485" in filename:
                    info["数据版本"] = "485眼版本"
                elif "499" in filename:
                    info["数据版本"] = "499眼版本"
                elif "20260315" in filename:
                    info["数据版本"] = "今日分析(2026-03-15)"
                    info["状态"] = "最新"

        except Exception as e:
            info["错误"] = str(e)[:100]

        return info

    # 收集文件信息
    print("\n🔍 正在收集文件信息...")

    # 主要目录
    base_dirs = [
        "/mnt/c/Users/CUI/Desktop/投稿、数据修改/03_Tables",
        "/mnt/c/Users/CUI/Desktop/投稿、数据修改/08_Unified_Analysis_20260315/unified_analysis_results_20260315"
    ]

    all_files = []
    for base_dir in base_dirs:
        if os.path.exists(base_dir):
            for root, dirs, files in os.walk(base_dir):
                for file in files:
                    if file.endswith('.xlsx'):
                        all_files.append(os.path.join(root, file))

    print(f"发现 {len(all_files)} 个Excel文件")

    # 为每个分析类型查找文件
    analysis_results = {}
    for analysis_name, keywords in analysis_types:
        matched_files = []
        for filepath in all_files:
            filename = os.path.basename(filepath)
            for keyword in keywords:
                if keyword in filename:
                    matched_files.append(filepath)
                    break

        # 按修改时间排序
        if matched_files:
            matched_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            latest_file = matched_files[0]
            info = check_file_info(latest_file, analysis_name)
            analysis_results[analysis_name] = {
                "最新文件": info,
                "所有文件": [check_file_info(f, analysis_name) for f in matched_files[:3]],  # 只保留前3个
                "文件总数": len(matched_files)
            }

    # 生成报告
    print("\n" + "=" * 80)
    print("📊 最新版统计分析结果清单")
    print("=" * 80)

    for analysis_name, keywords in analysis_types:
        if analysis_name in analysis_results:
            result = analysis_results[analysis_name]
            latest = result["最新文件"]

            print(f"\n{analysis_name}")
            print("-" * 40)
            print(f"✅ 最新文件: {latest['文件']}")
            print(f"   修改时间: {latest['修改时间']}")
            print(f"   文件大小: {latest['大小_KB']} KB")
            print(f"   数据版本: {latest['数据版本']}")
            print(f"   状态: {latest['状态']}")

            if "行数" in latest:
                print(f"   数据规模: {latest['行数']} 行 × {latest['列数']} 列")

            if "样本信息" in latest and latest["样本信息"] != "未知":
                print(f"   样本信息: {latest['样本信息']}")

            if "工作表" in latest:
                print(f"   包含工作表: {latest['工作表'][:3]}")  # 只显示前3个工作表

            # 显示关键结果摘要
            if analysis_name == "1. 描述性统计分析":
                print(f"   内容: 73个OCT指标的描述性统计（均值、标准差等）")
            elif analysis_name == "2. 组间比较（MDD vs Control）":
                print(f"   内容: MDD组与对照组的OCT指标比较（t检验/非参数检验）")
            elif analysis_name == "3. 相关性分析":
                print(f"   内容: OCT指标与PHQ-9评分的相关性分析")
            elif analysis_name == "4. ROC分析":
                print(f"   内容: OCT指标诊断抑郁的ROC曲线分析（AUC、敏感度、特异度）")
            elif analysis_name == "5. 多变量回归分析":
                print(f"   内容: 控制年龄、性别后的多变量回归分析")
            elif analysis_name == "6. 亚组分析":
                print(f"   内容: 性别、年龄、PHQ-9严重程度的亚组分析")

            # 如果有多个文件
            if result["文件总数"] > 1:
                print(f"   该分析共有 {result['文件总数']} 个版本，其他版本:")
                for file_info in result["所有文件"][1:]:  # 跳过第一个（最新）
                    print(f"     • {file_info['文件']} ({file_info['修改时间']}, {file_info['数据版本']})")

    print("\n" + "=" * 80)
    print("📈 基于数据版本的分析建议")
    print("=" * 80)

    # 分析数据版本一致性
    version_counts = {}
    for analysis_name in analysis_results:
        latest = analysis_results[analysis_name]["最新文件"]
        version = latest["数据版本"]
        version_counts[version] = version_counts.get(version, 0) + 1

    print(f"\n数据版本分布:")
    for version, count in version_counts.items():
        print(f"   {version}: {count} 个分析")

    # 检查一致性
    print(f"\n🔍 一致性检查:")
    if len(version_counts) == 1:
        print(f"   ✅ 所有分析使用相同数据版本: {list(version_counts.keys())[0]}")
    else:
        print(f"   ⚠️ 存在数据版本不一致:")
        for version, count in version_counts.items():
            print(f"      {version}: {count} 个分析")

    print(f"\n🎯 使用建议:")
    print(f"   1. 优先使用带有'20260315'或'今日分析'标记的文件")
    print(f"   2. 描述性统计和组间比较可能需要基于499眼版本更新")
    print(f"   3. 统一分析结果汇总.xlsx包含最全面的最新分析")
    print(f"   4. 注意ROC分析可能存在方法问题（AUC<0.5）")

    print("\n" + "=" * 80)
    print("💾 核心推荐文件")
    print("=" * 80)

    # 推荐核心文件
    core_files = [
        ("统一分析结果汇总.xlsx", 
         "08_Unified_Analysis_20260315/unified_analysis_results_20260315/统一分析结果汇总.xlsx",
         "包含多变量回归、亚组分析、ROC分析、相关性分析（基于463眼样本）"),

        ("All_OCT_Indicators_Descriptive_Statistics_Complete.xlsx",
         "03_Tables/All_OCT_Indicators_Descriptive_Statistics_Complete.xlsx",
         "73个OCT指标的完整描述性统计（基于485眼版本，需要验证）"),

        ("All_OCT_Indicators_Group_Comparison.xlsx",
         "03_Tables/All_OCT_Indicators_Group_Comparison.xlsx",
         "MDD与对照组的组间比较（基于485眼版本，需要验证）")
    ]

    for filename, rel_path, description in core_files:
        full_path = f"/mnt/c/Users/CUI/Desktop/投稿、数据修改/{rel_path}"
        if os.path.exists(full_path):
            mtime = datetime.fromtimestamp(os.path.getmtime(full_path)).strftime('%Y-%m-%d %H:%M')
            print(f"\n📄 {filename}")
            print(f"   位置: {rel_path}")
            print(f"   修改: {mtime}")
            print(f"   说明: {description}")

    print("\n" + "=" * 80)
    print("✅ 报告生成完成")
    print("=" * 80)


if __name__ == "__main__":
    main()