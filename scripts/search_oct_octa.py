#!/usr/bin/env python3
"""搜索OCT+OCTA联合分析研究"""

import requests

all_results = []

queries = [
    "OCT OCTA combined analysis retina",
    "multimodal retinal imaging OCT OCTA",
    "retinal capillary density OCTA neurological",
    "optical coherence tomography angiography multimodal",
    "OCTA retinal vessel density analysis",
]

for query in queries:
    url = "https://api.openalex.org/works"
    params = {
        "search": query,
        "filter": "publication_date:2015-01-01|2026-04-06",
        "sort": "cited_by_count:desc",
        "per-page": 10,
        "select": "id,title,publication_year,publication_date,cited_by_count,doi,primary_location"
    }
    
    try:
        resp = requests.get(url, params=params, timeout=20)
        data = resp.json()
        results = data.get("results", [])
        
        for work in results:
            doi = work.get("doi", "")
            title = work.get("title", "N/A")[:80]
            citations = work.get("cited_by_count", 0)
            date = work.get("publication_date", "N/A")
            primary_loc = work.get("primary_location", {})
            source = primary_loc.get("source", {})
            journal = source.get("display_name", "N/A") if source else "N/A"
            
            if citations > 0 and title != "N/A":
                all_results.append({
                    "title": title,
                    "citations": citations,
                    "date": date,
                    "journal": journal,
                    "doi": doi
                })
    except Exception as e:
        print(f"Query error: {e}")
        continue

# 去重并按引用排序
seen = set()
unique_results = []
for r in all_results:
    key = r["title"][:50]
    if key not in seen:
        seen.add(key)
        unique_results.append(r)

unique_results.sort(key=lambda x: x["citations"], reverse=True)

print("="*70)
print("🔬 OCT + OCTA 联合分析 - 高引用研究")
print("="*70)
print(f"共找到 {len(unique_results)} 篇独特高引用文献\n")

for i, r in enumerate(unique_results[:20], 1):
    print(f"{i}. {r['title']}")
    print(f"   📅 {r['date']} | ⭐ {r['citations']} | 📖 {r['journal'][:40]}")
    if r['doi']:
        print(f"   🔗 {r['doi']}")
    print()