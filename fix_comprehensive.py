import re

with open('/root/.openclaw/workspace/temp_draft_final.txt', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix 1: Replace m with μm in specific contexts (thickness measurements)
# Pattern: number ± number m, or number vs number m,
content = re.sub(r'(\d+\.\d+\s*±\s*\d+\.\d+)\s+m,', r'\1 μm,', content)
content = re.sub(r'(\d+\.\d+\s*±\s*\d+\.\d+)\s+m\)', r'\1 μm)', content)
content = re.sub(r'vs\.\s*(\d+\.\d+\s*±\s*\d+\.\d+)\s+m', r'vs. \1 μm', content)

# Fix 2: These P values are actual values, not "less than"
content = re.sub(r'P<0\.008\)', r'P=0.008)', content)
content = re.sub(r'P<0\.009\)', r'P=0.009)', content)

# Fix 3: drug-naive → drug-naïve
content = re.sub(r'drug-naive', r'drug-naïve', content)

# Write output
with open('/root/.openclaw/workspace/temp_draft_comprehensive.txt', 'w', encoding='utf-8') as f:
    f.write(content)

print("Done!")
