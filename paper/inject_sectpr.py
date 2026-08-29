#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Inject Letter+1in page setup into submission docx files (post-pandoc).

pandoc's default reference doc has an empty sectPr; every regeneration must
re-apply this. (Superscript citations are added by a separate step.)
"""
import zipfile, shutil, re, os

BASE = "/Users/luoxiaowen/Desktop/LKDock/LKina论文"
FILES = ["LKina_manuscript_EN_submission.docx", "LKina_manuscript_CN_submission.docx"]

SECTPR = ('<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>'
          '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" '
          'w:header="720" w:footer="720" w:gutter="0"/>'
          '<w:footnotePr><w:numRestart w:val="eachSect"/></w:footnotePr></w:sectPr>')


def inject(path):
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin:
        names = zin.namelist()
        doc = zin.read("word/document.xml").decode("utf-8")
        n = doc.count("<w:sectPr")
        if 'w:w="12240"' in doc:
            print(f"{os.path.basename(path)}: sectPr already present, skipped")
            return
        # replace the (single) empty sectPr
        doc2 = re.sub(r"<w:sectPr.*?</w:sectPr>", SECTPR, doc, count=1, flags=re.S)
        if doc2 == doc:  # no sectPr existed -> insert before </w:body>
            doc2 = doc.replace("</w:body>", SECTPR + "</w:body>")
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
            for name in names:
                zout.writestr(name, doc2 if name == "word/document.xml" else zin.read(name))
    shutil.move(tmp, path)
    print(f"{os.path.basename(path)}: sectPr injected (was {n})")


if __name__ == "__main__":
    for f in FILES:
        inject(os.path.join(BASE, f))
