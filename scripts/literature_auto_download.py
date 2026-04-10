#!/usr/bin/env python3
"""
文献自动下载与分析脚本
目标: 自动检索、下载、分析与OCT-MDD相关的最新文献
作者: 小cui科研助手
"""

import os
import json
import time
import requests
from datetime import datetime

# 配置
DOWNLOAD_DIR = "/root/.openclaw/workspace/downloads"
PROXY = "http://127.0.0.1:7890"
LOG_FILE = "/root/.openclaw/workspace/logs/literature_auto.log"

# 关键文献列表 (DOI, 文件名模板, 优先级)
TARGET_PAPERS = [
    ("10.1038/s41467-021-27604-x", "Vascular_BBB_depression_2022", "high"),
    ("10.1038/s41467-024-50309-w", "Eye_brain_connections_2024", "high"),
    ("10.1016/j.preteyeres.2025.101350", "Oculomics_2025", "medium"),
    ("10.3390/jcm13113081", "Depression_Eye_Disease_2024", "medium"),
    ("10.1093/schbul/sbae102", "Schizophrenia_Bipolar_Retina_2024", "medium"),
    ("10.1016/j.biopsych.2024.04.014", "Schizophrenia_Retinal_Microstructures_2024", "medium"),
    ("10.1038/s41392-024-01738-y", "MDD_mechanism_2024", "high"),  # 已下载
]

def log(msg):
    """日志记录"""
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    with open(LOG_FILE, "a") as f:
        f.write(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
    print(msg)

def download_pdf(doi, filename, max_retries=3):
    """下载PDF"""
    url = f"https://doi.org/{doi}"
    
    # 尝试多个下载源
    sources = [
        f"https://doi.org/{doi}",
        f"https://www.nature.com/articles/{doi}.pdf",
        f"https://www.ncbi.nlm.nih.gov/pmc/articles/{doi}/pdf",
    ]
    
    for attempt in range(max_retries):
        try:
            # 先解析DOI获取真实URL
            resp = requests.get(url, proxies={"http": PROXY, "https": PROXY}, 
                              allow_redirects=True, timeout=30)
            actual_url = resp.url
            
            # 检查是否为PDF链接
            if not actual_url.endswith('.pdf'):
                # 尝试直接添加.pdf
                actual_url = actual_url.rstrip('/') + '.pdf'
            
            # 下载
            filepath = os.path.join(DOWNLOAD_DIR, f"{filename}.pdf")
            if os.path.exists(filepath):
                log(f"跳过(已存在): {filename}")
                return True
                
            resp = requests.get(actual_url, proxies={"http": PROXY, "https": PROXY},
                              timeout=120, stream=True)
            
            if resp.status_code == 200:
                with open(filepath, 'wb') as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        f.write(chunk)
                size = os.path.getsize(filepath)
                log(f"✅ 下载成功: {filename} ({size/1024/1024:.1f}MB)")
                return True
            else:
                log(f"⚠️ 尝试{attempt+1}失败 HTTP {resp.status_code}: {filename}")
                
        except Exception as e:
            log(f"⚠️ 错误 {attempt+1}: {filename} - {str(e)[:50]}")
            time.sleep(2)
    
    log(f"❌ 下载失败: {filename}")
    return False

def search_and_download(keywords, max_results=5):
    """搜索并下载最新文献"""
    query = " AND ".join(keywords)
    url = "https://api.openalex.org/works"
    params = {
        "search": query,
        "filter": "publication_date:2020-01-01|2026-04-07",
        "sort": "cited_by_count:desc",
        "per-page": max_results
    }
    
    try:
        resp = requests.get(url, params=params, timeout=30)
        results = resp.json().get("results", [])
        
        for work in results:
            doi = work.get("doi", "").replace("https://doi.org/", "")
            title = work.get("title", "")[:50]
            citations = work.get("cited_by_count", 0)
            
            # 检查是否已下载
            safe_name = "".join(c for c in title if c.isalnum())[:30]
            
            if citations > 10 and doi:
                log(f"🔍 发现高引用文献: {title}... ({citations}引用)")
                # 这里只记录，不自动下载以避免版权问题
                
    except Exception as e:
        log(f"搜索错误: {e}")

def main():
    """主函数"""
    log("=" * 50)
    log("📚 文献自动下载任务开始")
    log("=" * 50)
    
    # 确保下载目录存在
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)
    
    # 1. 下载目标文献
    downloaded = 0
    for doi, filename, priority in TARGET_PAPERS:
        if download_pdf(doi, filename):
            downloaded += 1
        time.sleep(1)  # 避免请求过快
    
    # 2. 搜索新文献
    search_keywords = [
        ["retinal OCT", "depression", "biomarker"],
        ["eye brain", "cognitive", "imaging"],
        ["neurovascular", "depression", "blood brain barrier"]
    ]
    
    for keywords in search_keywords:
        search_and_download(keywords)
    
    log(f"📊 任务完成: 下载 {downloaded}/{len(TARGET_PAPERS)} 篇文献")
    log("=" * 50)

if __name__ == "__main__":
    main()