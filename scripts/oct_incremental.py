#!/usr/bin/env python3
"""眼科文献增量更新脚本"""
import json
import urllib.request
import urllib.parse

# 眼科关键词
keywords = 'OCT retina glaucoma macular AMD diabetic retinopathy'
params = {
    'search': keywords,
    'from_publication_date': '2025-10-01',
    'sort': 'cited_by_count:desc',
    'per-page': 30
}

url = 'https://api.openalex.org/works?' + urllib.parse.urlencode(params)
print(f"URL: {url[:80]}...")

try:
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)
except Exception as e:
    print(f"Error: {e}")
    exit(1)

results = data.get('results', [])
meta = data.get('meta', {})

print(f'\n=== 眼科文献检索 ===')
print(f'总匹配: {meta.get("count", 0)}')
print(f'返回: {len(results)}')

# 筛选引用>50
high = [w for w in results if w.get('cited_by_count', 0) > 50 and w.get('doi')]
print(f'引用>50: {len(high)}')
print()

# 输出文献
papers = []
for w in high[:10]:
    paper = {
        'title': w.get('title', ''),
        'doi': w.get('doi', ''),
        'cited': w.get('cited_by_count', 0),
        'date': w.get('publication_date', ''),
        'journal': w.get('primary_location', {}).get('source', {}).get('display_name', '')
    }
    papers.append(paper)
    print(f"[{paper['date']}] 引用:{paper['cited']}")
    print(f"  {paper['title'][:70]}")
    print(f"  {paper['journal'][:50]}")
    print(f"  DOI: {paper['doi']}")
    print()

# 保存到知识库
output_file = '/root/.openclaw/workspace/knowledge/oct_incremental_update.md'
output = f"""# 眼科文献增量更新 - 2026-04-14

## 检索摘要
- 检索关键词: OCT, retina, glaucoma, macular, AMD, diabetic retinopathy  
- 时间范围: 2025-10-01 至今
- 总匹配: {meta.get('count', 0)}
- 返回结果: {len(results)}
- 高影响力(引用>50): {len(high)}篇

## 高影响力文献

"""
for i, p in enumerate(papers, 1):
    output += f"""### {i}. {p['title'][:80]}
- **日期**: {p['date']}
- **引用**: {p['cited']}
- **期刊**: {p['journal']}
- **DOI**: {p['doi']}

"""

with open(output_file, 'w', encoding='utf-8') as f:
    f.write(output)

print(f"已保存到: {output_file}")