# -*- coding: utf-8 -*-
"""md2docx.py — 提交版 Markdown → .docx（无第三方依赖）

docx = OOXML zip 包：手写 document.xml + 打包。支持标题/表格/列表/粗体/引用/分隔线。
h2 起新页（四页结构）；正文 微软雅黑 10.5pt。

用法：python scripts/md2docx.py <input.md> <output.docx>
"""
import re
import sys
import zipfile
from pathlib import Path

NS = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'

RPR_BASE = '<w:rFonts w:ascii="Calibri" w:eastAsia="微软雅黑"/><w:sz w:val="21"/><w:szCs w:val="21"/>'
BOLD = "<w:b/><w:bCs/>"


def esc(s: str) -> str:
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def run(text: str, extra_rpr: str = "") -> str:
    return f'<w:r><w:rPr>{RPR_BASE}{extra_rpr}</w:rPr><w:t xml:space="preserve">{text}</w:t></w:r>'


def inline_runs(s: str, extra_rpr: str = "") -> str:
    """**bold** 与 `code` → run 序列（全部内容都在 run 内）"""
    out, pos = [], 0
    for m in re.finditer(r"\*\*(.+?)\*\*|`(.+?)`", s):
        out.append(run(esc(s[pos:m.start()]), extra_rpr))
        if m.group(1) is not None:
            out.append(run(esc(m.group(1)), BOLD + extra_rpr))
        else:
            out.append(run(esc(m.group(2)), extra_rpr))
        pos = m.end()
    out.append(run(esc(s[pos:]), extra_rpr))
    return "".join(out)


def para(text: str, size=21, bold=False, page_break=False, style=None) -> str:
    ppr = ""
    if style:
        ppr += f'<w:pStyle w:val="{style}"/>'
    if page_break:
        ppr += "<w:pageBreakBefore/>"
    rpr = f'<w:rFonts w:ascii="Calibri" w:eastAsia="微软雅黑"/><w:sz w:val="{size}"/><w:szCs w:val="{size}"/>'
    if bold:
        rpr += BOLD
    ppr_full = f"<w:pPr>{ppr}</w:pPr>" if ppr else ""
    return f"<w:p>{ppr_full}<w:r><w:rPr>{rpr}</w:rPr><w:t xml:space=\"preserve\">{esc(text)}</w:t></w:r></w:p>"


def rich_para(text: str, size=21, bold=False, page_break=False) -> str:
    """带行内粗体/代码的段落"""
    ppr = "<w:pageBreakBefore/>" if page_break else ""
    lead = f"<w:p><w:pPr>{ppr}</w:pPr>" if ppr else "<w:p>"
    if bold:
        lead += f"<w:r><w:rPr>{RPR_BASE}{BOLD}</w:rPr><w:t xml:space=\"preserve\">{esc(text)}</w:t></w:r>"
    else:
        lead += inline_runs(text)
    return lead + "</w:p>"


def table(rows) -> str:
    ncol = max(len(r) for r in rows)
    borders = "".join(
        f'<w:{b} w:val="single" w:sz="4" w:space="0" w:color="999999"/>'
        for b in ["top", "left", "bottom", "right", "insideH", "insideV"])
    out = [f'<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/><w:tblBorders>{borders}</w:tblBorders></w:tblPr>']
    for ri, row in enumerate(rows):
        out.append("<w:tr>")
        for ci in range(ncol):
            cell = row[ci] if ci < len(row) else ""
            bold = "<w:b/><w:bCs/>" if ri == 0 else ""  # 表头加粗
            out.append(
                f'<w:tc><w:tcPr><w:tcW w:w="0" w:type="auto"/></w:tcPr>'
                f'<w:p><w:pPr><w:pStyle w:val="CellText"/></w:pPr>'
                f'{inline_runs(cell, bold)}'
                f'</w:p></w:tc>')
        out.append("</w:tr>")
    out.append("</w:tbl>")
    return "".join(out)


def md2docx(text: str) -> bytes:
    lines = text.splitlines()
    body, i = [], 0
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            i += 1; continue
        if re.match(r"^---+$", ln.strip()):
            i += 1; continue
        if ln.startswith(">"):
            body.append(rich_para(ln[1:].strip(), size=19))
            i += 1; continue
        # 表格
        if ln.strip().startswith("|"):
            rows, started = [], False
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not started and len(cells) >= 2 and all(re.match(r"^:?-+:?$", c) for c in cells[1:]):
                    started = True  # 分隔行：跳过
                else:
                    rows.append(cells)
                    started = True
                i += 1
            body.append(table(rows)); continue
        m = re.match(r"^(#{1,4})\s+(.*)$", ln)
        if m:
            lvl, t = len(m.group(1)), m.group(2)
            size = {1: 30, 2: 26, 3: 23, 4: 21}[lvl]
            body.append(para(t, size=size, bold=True, page_break=(lvl == 2)))
            i += 1; continue
        if re.match(r"^\s*[-*]\s+", ln):
            items, indent = [], None
            while i < len(lines) and re.match(r"^\s*[-*]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*]\s+", "", lines[i])); i += 1
            for it in items:
                body.append(rich_para("• " + it, size=21))
            continue
        if re.match(r"^\s*\d+\.\s+", ln):
            items = []
            while i < len(lines) and re.match(r"^\s*\d+\.\s+", lines[i]):
                items.append(re.sub(r"^\s*\d+\.\s+", "", lines[i])); i += 1
            for it in items:
                body.append(rich_para("• " + it, size=21))
            continue
        body.append(rich_para(ln))
        i += 1

    document = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document {NS}>
<w:body>
{''.join(body)}
<w:sectPr><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr>
</w:body>
</w:document>"""

    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/></Types>"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>"""
    doc_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"/>"""

    buf = Path("md2docx_tmp.docx")
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("word/_rels/document.xml.rels", doc_rels)
        z.writestr("word/document.xml", document)
    data = buf.read_bytes()
    buf.unlink()
    return data


if __name__ == "__main__":
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    dst.write_bytes(md2docx(src.read_text(encoding="utf-8")))
    print("docx 已生成:", dst, f"({dst.stat().st_size // 1024} KB)")
