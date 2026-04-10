import re

with open('/root/.openclaw/workspace/temp_draft_fixed2.txt', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix: "μmm2" → "mm²" and "μmm3" → "mm³"
content = re.sub(r'μmm(\d)', r'mm\1', content)

# Fix: Missing β in regression coefficients description
content = re.sub(r'Regression coefficients \(\)', r'Regression coefficients (β)', content)

# Fix: Scientific notation P values like P=5.89104 → P=5.89×10⁻⁴
content = re.sub(r'P=(\d+\.\d+)(\d{2,3})\)', 
                 lambda m: f'P={float(m.group(1))*10**(-int(m.group(2))):.2e})', content)

# Fix: Missing = sign in regression (line 447 area)
content = re.sub(r'were significant \(= (\d+\.\d+)', r'were significant (β = \1', content)

# Fix: R² = → R²= (remove space)
content = re.sub(r'R² =', r'R²=', content)

# Write output
with open('/root/.openclaw/workspace/temp_draft_fixed3.txt', 'w', encoding='utf-8') as f:
    f.write(content)

print("Done!")
