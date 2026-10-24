#!/usr/bin/env python3
"""Draw the light theme's logo from the logo itself.

The header's logo, docs/assets/amap-logo.png, is drawn for a dark ground: its
lettering, rings and rules are light greys that vanish on a light page. This
writes docs/assets/amap-logo-light.png from .github/AMAP.png with the same
shapes and alpha: the greys become the site's ink, and the teals deepen so
they hold against a light ground.

    python3 tools/logo_light.py [.github/AMAP.png] [docs/assets/amap-logo-light.png]

Needs Pillow. The header shows it at 57 x 40, so it is written at 137 x 96,
the size of the dark one.
"""
import colorsys
import sys
from pathlib import Path

from PIL import Image

SIZE = (137, 96)


def recolour(src: Path, out: Path) -> None:
    im = Image.open(src).convert("RGBA").resize(SIZE, Image.LANCZOS)
    px = im.load()
    for y in range(im.height):
        for x in range(im.width):
            r, g, b, a = px[x, y]
            if a == 0:
                continue
            h, l, s = colorsys.rgb_to_hls(r / 255, g / 255, b / 255)
            if s < 0.18 or max(r, g, b) - min(r, g, b) < 24:
                # A grey: the lighter it was, the darker its ink, in the
                # site's ink hue.
                h2, l2, s2 = 160 / 360, 0.10 + max(0.0, 0.88 - l) * 0.70, 0.14
            else:
                # A teal or the slate: the same hue, deeper and a little richer.
                h2, l2, s2 = h, min(max(0.15 + (l - 0.35) * 0.55, 0.14), 0.40), min(1.0, s * 1.08 + 0.05)
            r2, g2, b2 = colorsys.hls_to_rgb(h2, l2, s2)
            px[x, y] = (round(r2 * 255), round(g2 * 255), round(b2 * 255), a)
    im.save(out, optimize=True)


if __name__ == "__main__":
    src = Path(sys.argv[1] if len(sys.argv) > 1 else ".github/AMAP.png")
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/assets/amap-logo-light.png")
    recolour(src, out)
    print(f"{out}: {SIZE[0]} x {SIZE[1]}")
