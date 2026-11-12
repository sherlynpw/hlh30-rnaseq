"""Measure printed text size of each manifest panel.

For each panel: parse every content stream, collect |d| from `a b c d e f Tm`
times the current Tf size, scale by placed_width / canvas_width, and report
the smallest, modal and largest printed point sizes.
"""
import re
import subprocess
import sys
from collections import Counter

import pikepdf

MANIFEST = "manifest.tsv"
NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
TM = re.compile(rf"({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s+Tm\b")
TF = re.compile(rf"/\S+\s+({NUM})\s+Tf\b")
TOK = re.compile(rf"(?:{NUM}\s+){{6}}Tm\b|/\S+\s+{NUM}\s+Tf\b")


def sizes(path):
    pdf = pikepdf.open(path)
    page = pdf.pages[0]
    w = float(page.mediabox[2]) - float(page.mediabox[0])
    out = []
    streams = []
    contents = page.obj.get("/Contents")
    if isinstance(contents, pikepdf.Array):
        streams = [s.read_bytes() for s in contents]
    elif contents is not None:
        streams = [contents.read_bytes()]
    for raw in streams:
        text = raw.decode("latin-1")
        tf = 1.0
        for m in TOK.finditer(text):
            s = m.group(0)
            if s.endswith("Tf"):
                tf = float(TF.match(s).group(1))
            else:
                mm = TM.match(s)
                a, b, c, d = (float(mm.group(i)) for i in range(1, 5))
                # handle rotated text: size is the column norm
                size = (c * c + d * d) ** 0.5 * tf
                if size > 0.5:
                    out.append(round(size, 2))
    return w / 72.0, out


rows = []
for line in open(MANIFEST):
    if line.startswith("#") or not line.strip():
        continue
    panel, src, placed, *_ = line.rstrip("\n").split("\t")
    path = f"out/{panel}.pdf"
    try:
        canvas, s = sizes(path)
    except Exception as e:  # noqa
        rows.append((panel, f"ERR {e}"))
        continue
    if not s:
        rows.append((panel, "no Tm text"))
        continue
    k = float(placed) / canvas
    printed = [round(x * k, 1) for x in s]
    cnt = Counter(printed)
    mode = cnt.most_common(1)[0][0]
    rows.append((panel, f"{float(placed):5.2f} {canvas:5.2f} {k:5.3f} min {min(printed):4.1f} mode {mode:4.1f} max {max(printed):4.1f}  n={len(printed)}  dist={dict(sorted(cnt.items()))}"))

print("panel        placed canvas ratio  printed pt")
for p, r in rows:
    print(f"{p:12s} {r}")
