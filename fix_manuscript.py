import re

# Read the original file
with open('/root/.openclaw/workspace/temp_draft.txt', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix 1: P0.001 → P<0.001 (and similar patterns)
content = re.sub(r'P(\d+\.\d+)', r'P<\1', content)

# Fix 2: Numerical format - add ± between mean and SD
# Pattern: number + number (like 271.4516.91) → number ± number
content = re.sub(r'(\d+\.\d{2})(\d+\.\d{2})', r'\1 ± \2', content)

# Fix 3: β symbol - add β before = in regression equations
content = re.sub(r'\(=-(\d+\.\d+)', r'(β = -\1', content)
content = re.sub(r', β=-(\d+\.\d+)', r', β = -\1', content)

# Fix 4: CI error 0.0597 → 0.597
content = re.sub(r'0\.0597', r'0.597', content)

# Fix 5: mean MTC → mean macular thickness (or MT)
content = re.sub(r'\bmean MTC\b', r'mean macular thickness', content)

# Fix 6: R2 → R²
content = re.sub(r'\bR2\b', r'R²', content)

# Fix 7: Remove duplicate "Results:" in Abstract
content = re.sub(r'(Results:.*?)(Results:)', r'\1', content, flags=re.DOTALL)

# Fix 8: Cross sectional → cross-sectional
content = re.sub(r'cross sectional', r'cross-sectional', content)

# Fix 9: drug naive → drug-naïve (and standardize)
content = re.sub(r'drug naive', r'drug-naïve', content)
content = re.sub(r'drug -naïve', r'drug-naïve', content)

# Fix 10: μm or uM → μm (fix missing multiplication sign)
content = re.sub(r'(\d+\.\d+)m(?!\u03bc)', r'\1 μm', content)

# Fix 11: Fix extra spaces around ±
content = re.sub(r'(\d+\.\d+)\s*±\s*(\d+\.\d+)', r'\1 ± \2', content)

# Fix 12: Fix patterns like "271.45 16.91" → "271.45 ± 16.91" in specific contexts
# This handles cases where the above pattern missed
lines = content.split('\n')
fixed_lines = []
for line in lines:
    # If line has "vs." and numbers, try to fix the ±
    if 'vs.' in line and re.search(r'\d+\.\d+\s+\d+\.\d+', line):
        line = re.sub(r'(\d+\.\d+)\s+(\d+\.\d+)\s+(?:vs\.)\s+(\d+\.\d+)\s+(\d+\.\d+)', 
                      r'\1 ± \2 vs. \3 ± \4', line)
    fixed_lines.append(line)
content = '\n'.join(fixed_lines)

# Write the fixed file
with open('/root/.openclaw/workspace/temp_draft_fixed.txt', 'w', encoding='utf-8') as f:
    f.write(content)

print("Done! Fixed file saved.")
