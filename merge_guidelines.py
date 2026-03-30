import sys

# Read original file
with open('ophthalmology_knowledge_base.md', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Read new guidelines
with open('new_guidelines.md', 'r', encoding='utf-8') as f:
    new_content = f.read()

# Find the line with "## 真实性验证方法"
for i, line in enumerate(lines):
    if line.strip() == '## 真实性验证方法':
        break
else:
    i = len(lines)

# Insert new content before that line
# Ensure there is an empty line before new section
if i > 0 and lines[i-1].strip() != '':
    lines.insert(i, '\n')
    i += 1

lines.insert(i, new_content + '\n\n')

# Write back
with open('ophthalmology_knowledge_base.md', 'w', encoding='utf-8') as f:
    f.writelines(lines)

print('Merged successfully')