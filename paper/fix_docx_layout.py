#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Post-process submission docx: dedup Table-5 caption + explicit fonts.

Fixes found during full layout audit (2026-08-28):
1. Markdown source had a duplicated "Table 5." caption (fixed in md); the
   regenerated docx still carries both -> remove one paragraph.
2. Fonts fell back to theme fonts (Aptos / 等线) because the last docx
   regeneration ran without the font post-processing step -> replace all
   theme-based rFonts in styles.xml with explicit fonts:
     EN: Times New Roman everywhere
     CN: eastAsia = SimSun, latin = Times New Roman
"""
import zipfile, shutil, re, sys, os

BASE = "/Users/luoxiaowen/Desktop/LKDock/LKina论文"

FONT_MAP = {
    "LKina_manuscript_EN_submission.docx": {
        "ascii": "Times New Roman", "hAnsi": "Times New Roman",
        "eastAsia": "Times New Roman", "cs": "Times New Roman",
    },
    "LKina_manuscript_CN_submission.docx": {
        "ascii": "Times New Roman", "hAnsi": "Times New Roman",
        "eastAsia": "SimSun", "cs": "Times New Roman",
    },
}

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def fix_styles(xml: str, fonts: dict) -> str:
    """Replace every theme-based rFonts with explicit font choices."""
    repl = (f'<w:rFonts w:ascii="{fonts["ascii"]}" w:hAnsi="{fonts["hAnsi"]}" '
            f'w:eastAsia="{fonts["eastAsia"]}" w:cs="{fonts["cs"]}"/>')
    # single-tag theme rFonts
    xml = re.sub(r'<w:rFonts[^>]*w:asciiTheme[^>]*/>', repl, xml)
    # any remaining rFonts without ascii (rare) -> normalize too
    xml = re.sub(r'<w:rFonts[^>]*w:eastAsiaTheme[^>]*/>', repl, xml)
    return xml


def dedup_captions(doc_xml: str) -> tuple[str, int]:
    """Remove consecutive duplicated caption paragraphs (Table N. / 表 N.)."""
    removed = 0
    # match a <w:p>...</w:p> whose visible text is a caption; then check the
    # immediately following paragraph for identical visible text
    para_re = re.compile(r'<w:p\b[^>]*>(?:(?!</w:p>).)*?</w:p>', re.S)

    def visible_text(p_xml: str) -> str:
        return "".join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', p_xml))

    out = []
    pos = 0
    paras = list(para_re.finditer(doc_xml))
    skip_until = -1
    for i, m in enumerate(paras):
        if m.start() < skip_until:
            continue
        txt = visible_text(m.group(0)).strip()
        if re.match(r'^(Table|表)\s*\d+[.:：]', txt) and i + 1 < len(paras):
            nxt = paras[i + 1]
            if visible_text(nxt.group(0)).strip() == txt:
                out.append(doc_xml[pos:m.end()])   # keep first
                pos = nxt.start()                   # skip duplicate
                skip_until = nxt.end()
                removed += 1
                continue
    out.append(doc_xml[pos:])
    return "".join(out), removed


def process(path: str, fonts: dict):
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin:
        names = zin.namelist()
        styles = zin.read("word/styles.xml").decode("utf-8")
        doc = zin.read("word/document.xml").decode("utf-8")
        styles2 = fix_styles(styles, fonts)
        doc2, removed = dedup_captions(doc)
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
            for n in names:
                if n == "word/styles.xml":
                    zout.writestr(n, styles2)
                elif n == "word/document.xml":
                    zout.writestr(n, doc2)
                else:
                    zout.writestr(n, zin.read(n))
    shutil.move(tmp, path)
    return removed


if __name__ == "__main__":
    for name, fonts in FONT_MAP.items():
        p = os.path.join(BASE, name)
        n = process(p, fonts)
        print(f"{name}: removed {n} duplicate caption(s); fonts -> "
              f"{fonts['ascii']} / ea:{fonts['eastAsia']}")
