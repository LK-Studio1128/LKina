#!/usr/bin/env python
"""Post-process pandoc LKina HTML (EN/ZH) into the LKDock SCI-manuscript style. v3

v3 adds IN-TEXT FIGURE PLACEMENT (per user request: figures must sit inside the
related subsection, not stacked in a trailing legend section):
  - every <figure class="fig"> is extracted and re-inserted at the END of its
    mapped subsection (anchor ids below), like the hand-made LKina论文.html;
  - S-figures stay where they were (supplementary area);
  - the figure style adopts the old card look: surface bg + rounded border +
    caption inside the card.

Caption harvesting decoded from pandoc's messy merge (unchanged from v2):
  - <p><img N/> - <strong>图 N+1.</strong> caption(N+1)</p>   (merged paragraphs)
  - <ul><li><strong>图 N.</strong> caption(N)</li></ul>        (surviving bullets)
"""
import re, sys, os, html as H

LANG_CSS_FONT = {
    'zh': "'Songti SC','SimSun',serif",
    'en': "'Times New Roman', Times, serif",
}
LABEL_WORD = {'zh': '图', 'en': 'Figure'}
# Graphical-abstract label per language (unnumbered figure, captioned inline).
GA_LABEL = {'zh': '图形摘要', 'en': 'Graphical Abstract'}

# Subsection anchor -> figures that must be placed at the END of that section.
# (zh_anchor, en_anchor): ids exist under the same heading in both files.
ANCHOR_FIGS = [
    ('金属蛋白重对接基准实测数据集',
     'metalloprotein-redocking-benchmarks-measured-datasets', ['6', '7']),
    ('系统性金属覆盖基准110-种金属模式',
     'systematic-metal-coverage-benchmark-110-metal-modes', ['1', '2']),
    ('功能族实测伪原子几何-bvs-水桥-金属作配体',
     'feature-family-measurements-pseudoatoms-bvs-water-bridge-metal-as-ligand', ['3', '4']),
    ('响应式共价预设66-全部通过',
     'reactive-covalent-presets-66', ['5']),
    ('评分函数对比参数扫描与热点靶点案例',
     'scoring-function-comparison-parameter-sweeps-and-hotspot-target-cases', ['S5', 'S6', 'S7', 'S8', 'S9', 'S10']),
]
H_RE = re.compile(r'<h[23][ >]')

CSS = """<style>
:root {{
  --ink:#1a1a1a; --muted:#5a5a5a; --line:#c9c9c9; --accent:#8b0000; --bg:#ffffff;
}}
* {{ box-sizing:border-box; }}
body {{
  font-family:{font};
  color:var(--ink); background:var(--bg);
  max-width:820px; margin:2rem auto; padding:0 1.4rem 4rem;
  line-height:1.75; font-size:1.02rem;
}}
header#title-block-header {{ text-align:center; margin-bottom:1rem; }}
h1.title {{ font-size:1.7rem; line-height:1.35; text-align:center; margin:1.2rem 0 .8rem; }}
h2 {{
  font-size:1.25rem; border-bottom:2px solid var(--ink); padding-bottom:.25rem;
  margin-top:2.2rem;
}}
h3 {{ font-size:1.12rem; margin-top:1.8rem; color:var(--accent); }}
h4 {{ font-size:1.03rem; }}
p {{ text-align:justify; margin:.85rem 0; hyphens:auto; }}
strong {{ font-weight:700; }}
code {{
  font-family:"SF Mono",Consolas,monospace; font-size:.88em;
  background:#f4f4f4; padding:.05em .3em; border-radius:3px;
}}
pre {{ background:#f7f6f3; border:1px solid #e2e0da; border-radius:4px;
  padding:.7rem .9rem; overflow-x:auto; }}
pre code {{ background:none; padding:0; font-size:.85em; }}
table {{
  width:100%; border-collapse:collapse; margin:1.2rem 0; font-size:.93rem;
  display:table; table-layout:auto;
}}
th, td {{
  border:1px solid var(--line); padding:.42rem .6rem;
  text-align:left; vertical-align:top;
}}
th {{ background:#f5f2ec; font-weight:700; }}
tbody tr:nth-child(even) {{ background:#fafaf8; }}
nav#TOC {{
  background:#faf9f7; border:1px solid #e2e0da; border-radius:4px;
  padding:.9rem 1.2rem; margin:1.4rem 0 2rem; font-size:.92rem;
}}
nav#TOC::before {{
  content:"{toc_title}"; display:block; font-weight:700; margin-bottom:.4rem;
  color:var(--accent);
}}
nav#TOC ul {{ padding-left:1.1em; margin:.2em 0; list-style:none; }}
nav#TOC > ul {{ padding-left:0; }}
nav#TOC a {{ text-decoration:none; color:var(--ink); }}
nav#TOC a:hover {{ color:var(--accent); }}
figure.fig {{
  margin:1.6rem 0; text-align:center; break-inside:avoid;
  background:#faf9f7; border:1px solid #e2e0da; border-radius:10px;
  padding:14px 18px;
}}
figure.fig img {{ max-width:100%; height:auto; border:1px solid #e2e0da; padding:4px; }}
figure.fig figcaption {{
  font-size:.9rem; color:var(--muted); margin-top:.5rem;
  text-align:justify; line-height:1.5;
}}
hr {{ border:none; border-top:1px solid var(--line); margin:2rem auto; width:60%; }}
ul, ol {{ padding-left:1.6rem; }}
ol li, ul li {{ margin:.45rem 0; text-align:justify; font-size:.95rem; }}
p.note {{ color:var(--muted); font-style:italic; }}
@media print {{
  body {{ max-width:100%; margin:0; font-size:10.5pt; }}
  h2 {{ break-after:avoid; }}
  table {{ break-inside:avoid; }}
  figure.fig {{ break-inside:avoid; }}
  nav#TOC {{ display:none; }}
}}
</style>"""

IMG_RE = re.compile(
    r'<img role="img" aria-label="(?P<label>[^"]+)" '
    r'src="(?P<uri>data:image/png;base64,[^"]+)" alt="[^"]*" */?>')

def process(src, dst, lang):
    w = LABEL_WORD[lang]
    s = open(src, encoding='utf-8').read()

    # ---------- 1. stylesheet ----------
    toc_title = '目录' if lang == 'zh' else 'Contents'
    s = re.sub(r'<style>.*?</style>',
               CSS.format(font=LANG_CSS_FONT[lang], toc_title=toc_title),
               s, count=1, flags=re.S)

    # ---------- 1b. fix MathJax font path (pandoc 3.9 bakes /output/... into bundle) ----------
    # Inject window.MathJax config BEFORE the inline webpack bundle runs; MathJax 3
    # merges user config at startup. Images stay embedded; only webfonts hit CDN
    # (and MathJax falls back to system fonts when offline).
    s = re.sub(
        r'(<script type="text/javascript">\(function\(\)\{"use strict";)',
        '<script>window.MathJax={"chtml":{"fontURL":'
        '"https://cdn.jsdelivr.net/npm/mathjax@3/es5/output/chtml/fonts/woff-v2"}};</script>'
        '\\1',
        s, count=1)

    # ---------- 2. harvest captions by label (from <li> bullets AND merged <p>) ----------
    caps = {}
    # 2a li bullets
    for m in re.finditer(
            rf'<li><strong>({re.escape(w)} [^<]+?)\.</strong>(.*?)</li>', s, re.S):
        caps[m.group(1).strip()] = m.group(2).strip()
    # 2a' graphical-abstract bullet (label does not start with the figure word)
    for m in re.finditer(
            rf'<li><strong>({re.escape(GA_LABEL[lang])})\.</strong>(.*?)</li>', s, re.S):
        caps[m.group(1).strip()] = m.group(2).strip()
    # 2b merged <p><img A/> - <strong>B.</strong> text</p>
    for m in re.finditer(
            rf'<strong>({re.escape(w)}[^<]+?)\.</strong>(.*?)(?=</p>|<strong>|$)', s, re.S):
        label = m.group(1).strip()
        # normalize whitespace-broken label e.g. "图\n8"
        label_norm = re.sub(r'\s+', ' ', label)
        caps[label_norm] = m.group(2).strip()

    # ---------- 3. rebuild every image as a proper figure ----------
    def make_fig(label, uri):
        cap = caps.get(label)
        if cap is None:
            cap_n = re.sub(r'\s+', ' ', label)
            cap = caps.get(cap_n)
        if cap is None:
            return None
        plain = untag(H.unescape(re.sub(r'<[^>]+>', '', cap)))
        alt = H.escape(plain)[:180]
        return (f'<figure class="fig"><img src="{uri}" alt="{alt}">'
                f'<figcaption><strong>{label}.</strong> {cap}</figcaption></figure>')

    # 3a. merged paragraphs: <p><img .../> - [stuff] </p>  -> keep only image as figure;
    #     the trailing caption fragment belongs to the NEXT label (already harvested).
    def fix_para(mm):
        inner = mm.group(1)
        im = IMG_RE.search(inner)
        if not im:
            return mm.group(0)
        fig = make_fig(im.group('label'), im.group('uri'))
        return fig if fig else mm.group(0)
    s = re.sub(r'<p>(<img role="img"[^>]+/>.*?)</p>', fix_para, s, flags=re.S)

    # 3b. old-style <figure> with dummy figcaption
    def fix_old_fig(mm):
        im = IMG_RE.search(mm.group(0))
        if not im:
            return mm.group(0)
        fig = make_fig(im.group('label'), im.group('uri'))
        return fig if fig else mm.group(0)
    s = re.sub(r'<figure>\s*<img role="img"[^>]+/>\s*<figcaption[^>]*>[^<]*</figcaption>\s*</figure>',
               fix_old_fig, s, flags=re.S)

    # 3c. remove now-orphaned caption bullets/fragments (they are inside figcaptions now)
    s = re.sub(rf'<li><strong>{re.escape(w)} [^<]+\.</strong>.*?</li>', '', s, flags=re.S)
    # graphical-abstract orphan bullet too
    s = re.sub(rf'<li><strong>{re.escape(GA_LABEL[lang])}\.</strong>.*?</li>', '', s, flags=re.S)
    # stray "- <strong>图 X.</strong> text" fragments directly before figures (paranoia)
    s = re.sub(rf'\s*-\s*<p class="note"><strong>{re.escape(w)}', '<p class="note"><strong>' + w, s)

    # 3d. drop empty ULs left over
    s = re.sub(r'<ul>\s*</ul>', '', s)

    # 3e. remove residual old-style <figure> blocks whose label already has a new figure
    new_labels = set(m.group(1) for m in re.finditer(
        r'<figure class="fig">.*?<figcaption><strong>([^<]+)</strong>', s, re.S))
    def drop_dup_fig(mm):
        im = IMG_RE.search(mm.group(0))
        if im:
            lbl = re.sub(r'\s+', ' ', im.group('label')).rstrip('.') + '.'
            if lbl in new_labels:
                return ''
        return mm.group(0)
    s = re.sub(r'<figure>\s*<img role="img"[^>]+/>\s*<figcaption[^>]*>[^<]*</figcaption>\s*</figure>',
               drop_dup_fig, s, flags=re.S)

    # ---------- 4. move figures into their related subsection (in-text placement) ----------
    # Extract every figure block; key by caption number ("图 1" / "Figure 1" -> '1').
    fig_blocks = {}
    def extract(mm):
        fm = re.search(r'<figcaption><strong>' + re.escape(w) + r'\s*([0-9S]+)\.',
                       mm.group(0))
        if fm:
            fig_blocks[fm.group(1)] = mm.group(0)
            return ''
        # Unnumbered figures (e.g. graphical abstract "图形摘要." / "Graphical
        # abstract.") must be kept inline, not dropped.
        return mm.group(0)
    s = re.sub(r'<figure class="fig">.*?</figure>', extract, s, flags=re.S)

    # Section body = between its heading and the next h2/h3 heading.
    for zh_anchor, en_anchor, nums in ANCHOR_FIGS:
        anchor = zh_anchor if lang == 'zh' else en_anchor
        hm = re.search(rf'(<h[23] id="{re.escape(anchor)}"[^>]*>.*?</h[23]>)', s, flags=re.S)
        if not hm:
            print(f'  WARN: anchor not found: {anchor}')
            continue
        nxt = H_RE.search(s, hm.end())
        insert_at = nxt.start() if nxt else len(s)
        insert_txt = '\n'.join(fig_blocks.pop(n, '') for n in nums if n in fig_blocks)
        if insert_txt:
            s = s[:insert_at] + '\n' + insert_txt + '\n' + s[insert_at:]
            print(f'  placed {" & ".join(n for n in nums if n)} -> #{anchor[:40]}')

    # Any remaining numbered figure that missed mapping: restore at previous position
    if fig_blocks:
        leftovers = [k for k in fig_blocks if k.isdigit()]
        print('  WARN leftover figures:', list(fig_blocks))
        for k, blk in fig_blocks.items():
            s += '\n' + blk

    # ---------- 5. drop the now-empty "figure legends" trailing section ----------
    legend_id = '图注' if lang == 'zh' else 'figure-legends'
    m = re.search(rf'<h2 id="{re.escape(legend_id)}">.*?(?=<h2)', s, flags=re.S)
    if m and '<figure class="fig">' not in m.group(0) \
            and len(re.sub(r'<[^>]+>|\s', '', m.group(0))) <= len(legend_id) + 1:
        s = s[:m.start()] + s[m.end():]
        print('  removed empty legend section:', legend_id)

    open(dst, 'w', encoding='utf-8').write(s)
    print('saved:', dst, '%.0f KB' % (os.path.getsize(dst) / 1024))
    n_fig = s.count('<figure class="fig">')
    n_old = len(re.findall(r'<figure>\s*<img', s))
    n_pimg = len(re.findall(r'<p><img', s))
    print(f'  figures: {n_fig} | old-style left: {n_old} | p-img left: {n_pimg}')

def untag(s):
    return (s.replace('&lt;', '<').replace('&gt;', '>')
             .replace('&amp;', '&').replace('&quot;', '"'))

if __name__ == '__main__':
    src, dst, lang = sys.argv[1], sys.argv[2], sys.argv[3]
    process(src, dst, lang)
