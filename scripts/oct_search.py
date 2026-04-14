#!/usr/bin/env python3
import json, urllib.request, urllib.parse

# 搜索眼科关键词
keywords = "retina OR glaucoma OR macular OR ophthalmology OR eye disease"
params = urllib.parse.urlencode({
    "search": keywords,
    "sort": "cited_by_count:desc",
    "per-page": 25
})

url = f"https://api.openalex.org/works?{params}"
print(f"Requesting: {url[:80]}...")

try:
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)
    
    results = data.get('results', [])
    meta = data.get('meta', {})
    
    print(f"\n=== 眼科文献检索结果 ===")
    print(f"匹配总数: {meta.get('count', 0)}")
    print(f"返回结果: {len(results)}")
    print()
    
    # 筛选引用>50
    high_impact = [w for w in results if w.get('cited_by_count', 0) > 50 and w.get('doi')]
    print(f"引用>50高影响力文献: {len(high_impact)}篇")
    print()
    
    for w in high_impact[:10]:
        title = w.get('title', '')[:70]
        doi = w.get('doi', '')
        cited = w.get('cited_by_count', 0)
        date = w.get('publication_date', 'N/A')
        print(f"[{date}] 引用:{cited}")
        print(f"  {title}")
        print(f"  DOI: {doi}")
        print()
        
except Exception as e:
    print(f"Error: {e}")