import re

with open('/root/.openclaw/workspace/temp_draft_fixed.txt', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix: "μmm2" → "mm²" and "μmm3" → "mm³"
content = re.sub(r'μmm(\d)', r'mm\1', content)

# Fix: Missing β in regression coefficients description
content = re.sub(r'Regression coefficients \(\)', r'Regression coefficients (β)', content)

# Fix: R² = → R²= (remove space)
content = re.sub(r'R² =', r'R²=', content)

# Fix: mean MTC → mean macular thickness
content = re.sub(r'\bmean MTC\b', r'mean macular thickness', content)

# Write output
with open('/root/.openclaw/workspace/temp_draft_final.txt', 'w', encoding='utf-8') as f:
    f.write(content)

print("Done!")
