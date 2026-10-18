#!/usr/bin/env python3
"""Join PNG frames into a GIF that loops: gif.py OUT.gif FRAME.png:MS [FRAME.png:MS ...]

Each frame gets a palette of its own, 256 colours by median cut, with no
dithering: a palette shared by every frame bands the soft gradient behind the
drawing, and dithering speckles it. tools/pictures.mjs runs this for the
README's animation of the front page's drawing; it needs Pillow, which
requirements-test.txt installs with matplotlib.
"""
import sys

from PIL import Image


def main(argv):
    if len(argv) < 3:
        print(__doc__.splitlines()[0], file=sys.stderr)
        return 2
    frames, durations = [], []
    for arg in argv[2:]:
        path, _, ms = arg.rpartition(":")
        if not path or not ms.isdigit():
            print(f"not FRAME.png:MS: {arg}", file=sys.stderr)
            return 2
        with Image.open(path) as im:
            frames.append(im.convert("RGB").quantize(
                colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))
        durations.append(int(ms))
    frames[0].save(argv[1], save_all=True, append_images=frames[1:],
                   duration=durations, loop=0)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
