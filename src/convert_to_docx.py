"""
Program: convert_to_docx.py
Deskripsi: Mengonversi TK1_Bagian_Yosua_Technical_Report.md ke berkas Microsoft Word (.docx)
           agar dapat diunggah ke Google Drive dan dibuka via Google Docs tanpa berantakan.
"""

import os
import re
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
md_path = os.path.join(base_dir, 'TK1_Bagian_Yosua_Technical_Report.md')
docx_path = os.path.join(base_dir, 'TK1_Technical_Report.docx')

doc = Document()

# Set standard A4 margins (1 inch = 2.54 cm)
sections = doc.sections
for s in sections:
    s.top_margin = Inches(1.0)
    s.bottom_margin = Inches(1.0)
    s.left_margin = Inches(1.0)
    s.right_margin = Inches(1.0)
    s.page_width = Inches(8.27)   # A4
    s.page_height = Inches(11.69)

# Base styling
style_normal = doc.styles['Normal']
style_normal.font.name = 'Times New Roman'
style_normal.font.size = Pt(11)
style_normal.font.color.rgb = RGBColor(0x22, 0x22, 0x22)
style_normal.paragraph_format.line_spacing = 1.15
style_normal.paragraph_format.space_after = Pt(6)

with open(md_path, 'r', encoding='utf-8') as f:
    text = f.read()

lines = text.split('\n')
i = 0
in_code = False
code_lines = []

def set_cell_background(cell, fill_hex):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{fill_hex}"/>')
    tcPr.append(shd)

def add_styled_paragraph(doc, text):
    p = doc.add_paragraph()
    # Simple inline parser for bold and math
    parts = re.split(r'(\*\*.*?\*\*|\*.*?\*|`.*?`)', text)
    for part in parts:
        if part.startswith('**') and part.endswith('**'):
            run = p.add_run(part[2:-2])
            run.bold = True
        elif part.startswith('*') and part.endswith('*'):
            run = p.add_run(part[1:-1])
            run.italic = True
        elif part.startswith('`') and part.endswith('`'):
            run = p.add_run(part[1:-1])
            run.font.name = 'Courier New'
            run.font.size = Pt(10)
        else:
            p.add_run(part)
    return p

while i < len(lines):
    line = lines[i].strip()
    
    # 1. Code block handling
    if line.startswith('```'):
        if not in_code:
            in_code = True
            code_lines = []
        else:
            in_code = False
            # Render code block in table or shaded box
            tbl = doc.add_table(rows=1, cols=1)
            tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
            cell = tbl.cell(0, 0)
            set_cell_background(cell, "F5F6F8")
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.0
            run = p.add_run('\n'.join(code_lines))
            run.font.name = 'Courier New'
            run.font.size = Pt(8.5)
            run.font.color.rgb = RGBColor(0x2C, 0x3E, 0x50)
            doc.add_paragraph()  # spacing
        i += 1
        continue
        
    if in_code:
        code_lines.append(lines[i])
        i += 1
        continue
        
    if not line:
        i += 1
        continue
        
    # 2. Headings
    if line.startswith('# '):
        h = doc.add_heading(line[2:].strip(), level=1)
        h.paragraph_format.space_before = Pt(14)
        h.paragraph_format.space_after = Pt(8)
        i += 1
        continue
    elif line.startswith('## '):
        h = doc.add_heading(line[3:].strip(), level=2)
        h.paragraph_format.space_before = Pt(12)
        h.paragraph_format.space_after = Pt(6)
        i += 1
        continue
    elif line.startswith('### '):
        h = doc.add_heading(line[4:].strip(), level=3)
        h.paragraph_format.space_before = Pt(10)
        h.paragraph_format.space_after = Pt(4)
        i += 1
        continue
    elif line.startswith('#### '):
        p = doc.add_paragraph()
        run = p.add_run(line[5:].strip())
        run.bold = True
        run.font.size = Pt(11)
        p.paragraph_format.space_before = Pt(8)
        p.paragraph_format.space_after = Pt(4)
        i += 1
        continue
        
    # 3. Horizontal rules
    if line == '---':
        i += 1
        continue
        
    # 4. Images
    img_match = re.match(r'!\[(.*?)\]\((.*?)\)', line)
    if img_match:
        caption = img_match.group(1)
        img_url = img_match.group(2)
        # Extract relative or absolute local path
        local_img = img_url.replace('file:///', '').replace('%20', ' ')
        # If relative to figures
        if not os.path.exists(local_img):
            fname = os.path.basename(local_img)
            local_img = os.path.join(base_dir, 'figures', fname)
            
        if os.path.exists(local_img):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(8)
            p.paragraph_format.space_after = Pt(2)
            run = p.add_run()
            run.add_picture(local_img, width=Inches(6.2))
            
            p_cap = doc.add_paragraph()
            p_cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_cap.paragraph_format.space_after = Pt(10)
            run_cap = p_cap.add_run(f"{caption}")
            run_cap.italic = True
            run_cap.font.size = Pt(9.5)
            run_cap.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
        i += 1
        continue
        
    # 5. Tables
    if line.startswith('|') and '|' in line[1:]:
        table_lines = []
        while i < len(lines) and lines[i].strip().startswith('|'):
            table_lines.append(lines[i].strip())
            i += 1
            
        if len(table_lines) >= 2:
            # Parse header
            header = [c.strip() for c in table_lines[0].split('|')[1:-1]]
            # Row 1 is divider (|:---|:---|)
            data_rows = []
            for t_line in table_lines[2:]:
                data_rows.append([c.strip() for c in t_line.split('|')[1:-1]])
                
            n_cols = len(header)
            table = doc.add_table(rows=len(data_rows) + 1, cols=n_cols)
            table.alignment = WD_TABLE_ALIGNMENT.CENTER
            table.style = 'Table Grid'
            
            # Format header
            for col_idx, h_text in enumerate(header):
                cell = table.cell(0, col_idx)
                cell.paragraphs[0].paragraph_format.space_after = Pt(2)
                cell.paragraphs[0].paragraph_format.space_before = Pt(2)
                run = cell.paragraphs[0].add_run(h_text.replace('$', ''))
                run.bold = True
                run.font.size = Pt(9.5)
                set_cell_background(cell, "EAECEE")
                
            # Format data rows
            for row_idx, r_data in enumerate(data_rows):
                for col_idx in range(min(n_cols, len(r_data))):
                    cell = table.cell(row_idx + 1, col_idx)
                    cell.paragraphs[0].paragraph_format.space_after = Pt(2)
                    cell.paragraphs[0].paragraph_format.space_before = Pt(2)
                    cell_text = r_data[col_idx].replace('**', '')
                    run = cell.paragraphs[0].add_run(cell_text)
                    run.font.size = Pt(9)
                    if row_idx % 2 == 1:
                        set_cell_background(cell, "F8F9FA")
            
            doc.add_paragraph()  # spacing after table
        continue
        
    # 6. Quotes / Callouts
    if line.startswith('> '):
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.4)
        run = p.add_run(line[2:].strip().replace('*', ''))
        run.italic = True
        run.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
        i += 1
        continue
        
    # 7. Lists
    if line.startswith('- ') or line.startswith('* '):
        p = doc.add_paragraph(style='List Bullet')
        p.paragraph_format.space_after = Pt(3)
        add_styled_paragraph(doc, line[2:].strip())
        i += 1
        continue
        
    # Regular paragraph
    add_styled_paragraph(doc, line)
    i += 1

doc.save(docx_path)
print(f"File Word berhasil dibuat di: {docx_path}")
