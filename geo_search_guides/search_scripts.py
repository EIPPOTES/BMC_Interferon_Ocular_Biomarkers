
# Python搜索GEO的示例代码
import requests

def search_geo(query):
    """搜索GEO数据库"""
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
