"""Draw build_contrast_analysis() as a factory for man/figures/prolfqua_factory.png.

Usage:
  python3 data-raw/make_factory_svg.py /tmp/prolfqua_factory.svg
  rsvg-convert -w 2920 /tmp/prolfqua_factory.svg -o man/figures/prolfqua_factory.png
The README and the vignettes prolfqua and ModellingMethods use the PNG; pkgdown
copies man/figures to the site, but not SVG figures.
The lanes, their models and colours come from inst/extdata/facade_styles.tsv
(one row per facade, in lane order); add a row there when a facade is added to
the registry.
"""
import csv
import os
import sys

out = sys.argv[1]

styles = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inst", "extdata", "facade_styles.tsv")
lanes = []
with open(styles, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if not lanes or lanes[-1][0] != row["family"]:
            lanes.append((row["family"], [], row["fill"], row["colour"]))
        lanes[-1][1].append(row["method"])

W, H = 1460, 720
MONO = "Menlo, Consolas, 'DejaVu Sans Mono', monospace"
SANS = "Helvetica, Arial, sans-serif"
CHAR = 7.4  # approx. width of a 12px monospace glyph

el = []
add = el.append
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="{SANS}">')
add(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')

# smoke and chimney (drawn before the roof)
for cx, cy, r in [(501, 58, 12), (522, 40, 15), (551, 26, 18)]:
    add(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="#d9d9d9"/>')
add('<rect x="488" y="66" width="26" height="110" fill="#7d4436" stroke="#4a2820" stroke-width="2"/>')

# factory body and sawtooth roof
add('<polygon points="290,178 290,128 380,178 380,128 470,178 470,128 560,178" fill="#9c4433" stroke="#4a2820" stroke-width="2"/>')
add('<rect x="290" y="178" width="270" height="470" fill="#b5523b" stroke="#4a2820" stroke-width="2"/>')
add('<rect x="300" y="196" width="228" height="40" rx="4" fill="#ffffff" stroke="#4a2820" stroke-width="2"/>')
add(f'<text x="414" y="221" text-anchor="middle" font-family="{MONO}" font-size="13" font-weight="bold" fill="#222">build_contrast_analysis()</text>')
for row_y in (262, 318):
    for x in (308, 358, 408, 458):
        add(f'<rect x="{x}" y="{row_y}" width="32" height="34" fill="#f2d16b" stroke="#4a2820" stroke-width="1.5"/>')
add('<rect x="312" y="578" width="44" height="70" fill="#5a3329" stroke="#4a2820" stroke-width="2"/>')

# input card: the R call, and an arrow into the factory
add('<rect x="20" y="318" width="250" height="174" rx="8" fill="#eef3fb" stroke="#5b7db1" stroke-width="2"/>')
add('<text x="34" y="344" font-size="14" font-weight="bold" fill="#2b4a7a">R</text>')
call = [(34, "build_contrast_analysis("), (50, "lfqdata,"), (50, '"~ group_",'), (50, "contrasts,"), (50, 'method = "lm_impute"'), (34, ")")]
for i, (x, line) in enumerate(call):
    add(f'<text x="{x}" y="{370 + i * 22}" font-family="{MONO}" font-size="13" fill="#222">{line}</text>')
add('<line x1="270" y1="405" x2="280" y2="405" stroke="#5b7db1" stroke-width="3"/>')
add('<polygon points="280,397 292,405 280,413" fill="#5b7db1"/>')

# conveyor lanes with crates
for i, (label, models, fill, stroke) in enumerate(lanes):
    cy = 204 + i * 64
    add(f'<rect x="540" y="{cy - 24}" width="20" height="44" fill="#3b2a25"/>')
    add(f'<rect x="560" y="{cy + 14}" width="{W - 580}" height="7" rx="3" fill="#555555"/>')
    for rx in range(574, W - 20, 30):
        add(f'<circle cx="{rx}" cy="{cy + 26}" r="4.5" fill="#9a9a9a"/>')
    add(f'<text x="572" y="{cy + 4}" font-size="13" font-style="italic" fill="#444">{label}</text>')
    x = 700
    for m in models:
        w = len(m) * CHAR + 18
        add(f'<rect x="{x:.1f}" y="{cy - 16}" width="{w:.1f}" height="29" rx="3" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>')
        add(f'<text x="{x + w / 2:.1f}" y="{cy + 3}" text-anchor="middle" font-family="{MONO}" font-size="12" fill="#222">{m}</text>')
        x += w + 10
    assert x < W - 10, (label, x)

# common interface
add(f'<text x="700" y="672" font-size="14" fill="#333">every model rolls out with the same interface:</text>')
for x, f in zip((700, 840, 960, 1080), ("get_contrasts()", "get_Plotter()", "get_missing()", "to_wide()")):
    add(f'<text x="{x}" y="696" font-family="{MONO}" font-size="13" fill="#222">{f}</text>')
add('</svg>')

with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(el) + "\n")
print("models:", sum(len(l[1]) for l in lanes))
