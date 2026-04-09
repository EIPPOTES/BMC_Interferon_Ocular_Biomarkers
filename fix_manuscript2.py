import re

# Read the original file
with open('/root/.openclaw/workspace/temp_draft_fixed.txt', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix: m → μm (but not in words like "mm", "mm2", "mm3")
# Only for specific patterns like "271.45 ± 16.91 m," or "vs. 278.19 ± 14.89 m"
content = re.sub(r'(\d+\.\d+\s*±\s*\d+\.\d+)\s+m,(?!\s*\d)', r'\1 μm,', content)
content = re.sub(r'(\d+\.\d+\s*±\s*\d+\.\d+)\s+m(?!\u03bc)', r'\1 μm', content)

# Fix: P<0.008 → P=0.008 (these are specific values, not "less than")
content = re.sub(r'P<0\.008\)', r'P=0.008)', content)
content = re.sub(r'P<0\.009\)', r'P=0.009)', content)

# Fix: Drug naive → drug-naïve consistently
content = re.sub(r'drug-naive', r'drug-naïve', content)

# Fix remaining number patterns like "28.813.5" → "28.8 ± 13.5"
content = re.sub(r'(\d+\.\d)(\d+\.\d+)', r'\1 ± \2', content)

# Fix: "P<0.58" → "P=0.58" 
content = re.sub(r'P<0\.58\)', r'P=0.58)', content)

# Fix: "P<0.48" → "P=0.48"
content = re.sub(r'P<0\.48\)', r'P=0.48)', content)

# Fix: Add space after period before new sentence
content = re.sub(r'\)([A-Z])', r') \1', content)

# Write the fixed file
with open('/root/.openclaw/workspace/temp_draft_fixed2.txt', 'w', encoding='utf-8') as f:
    f.write(content)

print("Done!")
