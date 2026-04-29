#!/bin/bash
# 提取所有考试资料内容并整合

echo "开始提取考试资料..."

# 提取PDF内容
echo "=== 提取打印-合辑-硕士眼科学试题 ==="
pdftotext "/root/.openclaw/workspace/exam_materials/打印-合辑-硕士眼科学试题 - 副本.pdf" - 2>/dev/null > /tmp/exam_full.txt

echo "=== 提取专业课考核100题 ==="
pdftotext "/root/.openclaw/workspace/exam_materials/专业课考核100题.pdf" - 2>/dev/null >> /tmp/exam_full.txt

echo "=== 提取眼科问答题 ==="
pdftotext "/root/.openclaw/workspace/exam_materials/眼科问答题.pdf" - 2>/dev/null >> /tmp/exam_full.txt

echo "=== 提取眼科十大急症 ==="
pdftotext "/root/.openclaw/workspace/exam_materials/眼科十大急症  limeng.pdf" - 2>/dev/null >> /tmp/exam_full.txt

echo "=== 提取名词解释 ==="
pdftotext "/root/.openclaw/workspace/exam_materials/打印-名词解释-打印1.pdf" - 2>/dev/null >> /tmp/exam_full.txt

echo "=== 提取历年题答案总结 ==="
pdftotext "/root/.openclaw/workspace/exam_materials/打印-硕士眼科学-历年题答案总结.pdf" - 2>/dev/null >> /tmp/exam_full.txt

echo "=== 提取分章节总结 ==="
pdftotext "/root/.openclaw/workspace/exam_materials/眼科分章节总结资料.pdf" - 2>/dev/null >> /tmp/exam_full.txt

echo "提取完成，总行数:"
wc -l /tmp/exam_full.txt