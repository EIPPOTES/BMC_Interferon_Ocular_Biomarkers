#!/usr/bin/env python3
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.style import WD_STYLE_TYPE
import re

def markdown_to_docx(md_file, docx_file):
    """Convert Markdown file to DOCX with basic formatting."""
    
    # Read Markdown content
    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Create document
    doc = Document()
    
    # Set default style
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Times New Roman'
    font.size = Pt(12)
    style.paragraph_format.line_spacing = 1.5
    style.paragraph_format.space_after = Pt(6)
    
    # Split content into lines
    lines = content.split('\n')
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Skip empty lines
        if not line:
            i += 1
            continue
            
        # Check for title
        if line.startswith('# Title'):
            # Get the actual title from next line
            if i + 1 < len(lines):
                title_line = lines[i + 1].strip()
                while not title_line and i + 2 < len(lines):
                    i += 1
                    title_line = lines[i].strip()
                
                if title_line:
                    title = doc.add_heading(title_line, level=0)
                    for run in title.runs:
                        run.font.name = 'Times New Roman'
                        run.font.size = Pt(14)
                        run.font.bold = True
                    i += 1
            else:
                # Fallback title
                title = doc.add_heading('OCT-MDD Manuscript', level=0)
                for run in title.runs:
                    run.font.name = 'Times New Roman'
                    run.font.size = Pt(14)
                    run.font.bold = True
        
        # Check for Abstract section
        elif line.startswith('# Abstract'):
            doc.add_heading('Abstract', level=1)
            # Collect abstract content
            i += 1
            abstract_lines = []
            while i < len(lines) and not lines[i].strip().startswith('# '):
                if lines[i].strip():
                    abstract_lines.append(lines[i].strip())
                i += 1
            continue  # Don't increment i again
            
        # Check for other sections (level 1 headers)
        elif line.startswith('# ') and not line.startswith('# Title'):
            section_title = line[2:].strip()
            doc.add_heading(section_title, level=1)
            
        # Check for subsections (level 2 headers)
        elif line.startswith('## '):
            subsection_title = line[3:].strip()
            doc.add_heading(subsection_title, level=2)
            
        # Check for sub-subsections (level 3 headers)
        elif line.startswith('### '):
            subsubsection_title = line[4:].strip()
            doc.add_heading(subsubsection_title, level=3)
            
        # Regular paragraph
        else:
            # Skip separator lines
            if line.startswith('---') or line.startswith('***') or line.startswith('___'):
                i += 1
                continue
                
            # Handle bold formatting
            line_text = line
            # Convert **bold** to bold
            line_text = re.sub(r'\*\*(.*?)\*\*', r'\1', line_text)
            
            # Create paragraph
            p = doc.add_paragraph(line_text)
            
            # Check if it contains keywords like Background:, Methods:, etc.
            if line_text.startswith('Background:') or line_text.startswith('Methods:') or \
               line_text.startswith('Results:') or line_text.startswith('Conclusions:'):
                for run in p.runs:
                    if 'Background:' in run.text or 'Methods:' in run.text or \
                       'Results:' in run.text or 'Conclusions:' in run.text:
                        run.bold = True
        
        i += 1
    
    # Save document
    doc.save(docx_file)
    print(f"Converted {md_file} to {docx_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python convert_to_docx.py <input.md> <output.docx>")
        sys.exit(1)
    
    markdown_to_docx(sys.argv[1], sys.argv[2])
