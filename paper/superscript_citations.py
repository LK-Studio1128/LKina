#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Superscript in-text citations [n] / [n,m] / [n-m] in submission docx.

Body-only: after the References / 参考文献 heading, citations stay plain.
Pattern proven in submission v5 (commit d272569): vertAlign superscript + sz 18.
"""
import zipfile, shutil, re, os, sys

BASE = "/Users/luoxiaowen/Desktop/LKDock/LKina论文"
FILES = ["LKina_manuscript_EN_submission.docx", "LKina_manuscript_CN_submission.docx"]

CIT = re.compile(r'\[(\d{1,2}(?:[,，]\s*\d{1,2})*(?:-\d{1,2})?)\]')


def find_refs_split(doc_xml: str) -> int:
    """Return char index where the references section starts (last occurrence
    of a heading whose text is References / 参考文献)."""
    last = -1
    for m in re.finditer(r'<w:p\b[^>]*>(?:(?!</w:p>).)*?</w:p>', doc_xml, re.S):
        txt = "".join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', m.group(0))).strip()
        if txt in ("References", "参考文献") or txt.startswith("References") and len(txt) < 15:
            last = m.end()
    return last


def process(path):
    with zipfile.ZipFile(path) as zin:
        names = zin.namelist()
        doc = zin.read("word/document.xml").decode("utf-8")
    split = find_refs_split(doc)
    head, tail = doc[:split], doc[split:]

    def repl(mm):
        inner = mm.group(1)
        rpr = '<w:rPr><w:vertAlign w:val="superscript"/><w:sz w:val="18"/></w:rPr>'
        return f'<w:r>{rpr}<w:t>[{inner}]</w:t></w:r>'

    # Citations live inside runs: <w:r>...<w:t>text [6,7] more</w:t></w:r>
    # Split run text around citations, superscripting the bracket parts only.
    n_marks = 0
    out_parts = []
    pos = 0
    # iterate over runs in head
    for rm in re.finditer(r'<w:r\b[^>]*>(?:(?!</w:r>).)*?</w:r>', head, re.S):
        run = rm.group(0)
        tm = re.search(r'(<w:t[^>]*>)([^<]*)(</w:t>)', run)
        if not tm or '[' not in tm.group(2):
            continue
        text = tm.group(2)
        if not CIT.search(text):
            continue
        # find rPr in run (copy into each new run, preserving existing props)
        rprm = re.search(r'<w:rPr>.*?</w:rPr>', run, re.S)
        base_rpr = rprm.group(0) if rprm else ''
        sup_rpr = ('<w:rPr><w:vertAlign w:val="superscript"/><w:sz w:val="18"/>'
                   + ''.join(re.findall(r'<w:rFonts[^>]*/>', base_rpr))
                   + '</w:rPr>')
        pieces = []
        last = 0
        for cm in CIT.finditer(text):
            pre = text[last:cm.start()]
            if pre:
                pieces.append(f'<w:r>{base_rpr}<w:t xml:space="preserve">{pre}</w:t></w:r>')
            pieces.append(f'<w:r>{sup_rpr}<w:t>[{cm.group(1)}]</w:t></w:r>')
            n_marks += 1
            last = cm.end()
        rest = text[last:]
        if rest:
            pieces.append(f'<w:r>{base_rpr}<w:t xml:space="preserve">{rest}</w:t></w:r>')
        out_parts.append((rm.start(), rm.end(), ''.join(pieces)))

    # apply replacements from the end
    for start, end, repl_txt in reversed(out_parts):
        head = head[:start] + repl_txt + head[end:]

    doc2 = head + tail
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin:
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
            for name in names:
                zout.writestr(name, doc2 if name == "word/document.xml" else zin.read(name))
    shutil.move(tmp, path)
    print(f"{os.path.basename(path)}: {n_marks} citations superscripted "
          f"(refs split at {split})")


if __name__ == "__main__":
    for f in FILES:
        process(os.path.join(BASE, f))
