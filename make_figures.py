"""Regenerate manuscript figures as 600 dpi PNGs from simulation data.

Run from the folder containing chapter3_results:
    python make_figures.py

Requires NumPy 2 or newer and Pillow 10.1 or newer.
Output images are 4224 x 1488 pixels with 600 dpi metadata.
All plot elements are rendered directly at the output resolution.
"""
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path("chapter3_results")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

DPI = 600
SCALE = 2.4
IMAGE_SIZE = (4224, 1488)
PALETTE = ["#21577b", "#ab4e22", "#397548", "#6a407e", "#202020"]
FONT = "C:/Windows/Fonts/arial.ttf"
BOLD = "C:/Windows/Fonts/arialbd.ttf"
FIGURE_NUMBERS = {
    "baseline": 1, "kinetic": 2, "kinetic09": 3, "kinetic08": 4,
    "near": 5, "joint1809": 6, "joint1609": 7, "joint1808": 8,
    "joint1608": 9, "joint": 10, "jointnear": 11, "outlet": 12,
    "inventory": 13, "variant": 14, "refinement": 15,
}


def load(alpha, beta, tag="fine"):
    with np.load(ROOT / f"{tag}_{alpha:g}_{beta:g}.npz") as data:
        return dict(data)


def font(size, bold=False):
    pixel_size = round(size * SCALE)
    for name in (BOLD if bold else FONT,
                 "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, pixel_size)
        except OSError:
            continue
    return ImageFont.load_default(size=pixel_size)


def chart(name, panels):
    image = Image.new("RGB", IMAGE_SIZE, "white")
    draw = ImageDraw.Draw(image)
    panel_width = 880 * SCALE
    grid_width = round(2 * SCALE)
    curve_width = round(5 * SCALE)

    for index, (title, xlabel, ylabel, xmax, ymax, series) in enumerate(panels):
        offset = index * panel_width
        left, right = offset + 108 * SCALE, offset + 845 * SCALE
        bottom, top = 485 * SCALE, 70 * SCALE
        draw.text((offset + 104 * SCALE, 22 * SCALE),
                  f"({'ab'[index]}) {title}", fill="#111111", font=font(28, True))

        for j in range(5):
            fraction = j / 4
            x = left + (right - left) * fraction
            y = bottom - (bottom - top) * fraction
            draw.line((x, bottom, x, top), fill="#dddddd", width=grid_width)
            draw.line((left, y, right, y), fill="#dddddd", width=grid_width)
            draw.text((x, bottom + 10 * SCALE), f"{xmax * fraction:g}",
                      anchor="mt", fill="#202020", font=font(23))
            draw.text((left - 11 * SCALE, y), f"{ymax * fraction:.3g}",
                      anchor="rm", fill="#202020", font=font(23))

        draw.rectangle((left, top, right, bottom), outline="#444444", width=grid_width)
        draw.text(((left + right) / 2, 535 * SCALE), xlabel,
                  anchor="mt", fill="#222222", font=font(26))

        vertical_label = Image.new("RGBA", (round(390 * SCALE), round(34 * SCALE)),
                                   (255, 255, 255, 0))
        label_draw = ImageDraw.Draw(vertical_label)
        label_draw.text((195 * SCALE, 17 * SCALE), ylabel,
                        anchor="mm", fill="#222222", font=font(25))
        vertical_label = vertical_label.rotate(90, expand=True)
        image.paste(vertical_label, (round(offset + 6 * SCALE), round(83 * SCALE)),
                    vertical_label)

        for j, (label, x_values, y_values) in enumerate(series):
            points = [
                (float(left + (right - left) * x / xmax),
                 float(bottom - (bottom - top) * y / ymax))
                for x, y in zip(x_values, y_values)
            ]
            draw.line(points, fill=PALETTE[j % len(PALETTE)], width=curve_width)

        legend_font = font(22)
        widths = [draw.textlength(label, font=legend_font) + 60 * SCALE
                  for label, _, _ in series]
        x = offset + (panel_width - sum(widths)) / 2
        y = 587 * SCALE
        if x < offset or x + sum(widths) > offset + panel_width:
            raise ValueError(f"Legend exceeds the panel width in {name}")
        for j, (label, _, _) in enumerate(series):
            draw.line((x, y, x + 25 * SCALE, y),
                      fill=PALETTE[j % len(PALETTE)], width=curve_width)
            draw.text((x + 33 * SCALE, y), label, anchor="lm",
                      fill="#202020", font=legend_font)
            x += widths[j]

    number = FIGURE_NUMBERS[name]
    output_path = OUTPUT_DIR / f"Figure_{number:02d}_{name}.png"
    image.save(output_path, dpi=(DPI, DPI))
    print(f"Saved {output_path}")

def profiles(name,cases,t=4800):
    records=[(label,load(a,b,tag)) for label,a,b,tag in cases]
    def lines(key,scale):
        return [(label,r['x'],r[key][np.searchsorted(r['t'],t)]/scale) for label,r in records]
    chart(name,[
        ('Suspension','Distance x (m)','Normalized concentration c/c0',0.4,1,lines('c',.05)),
        ('Deposit','Distance x (m)','Deposit fraction ρ',0.4,.104,lines('deposit',1))])
def histories(name,cases,points=(.02,.04)):
    records=[(label,load(a,b,tag)) for label,a,b,tag in cases]
    chart(name,[(f'x = {x:g} m','Time (s)','Deposit fraction ρ',4800,.104,
        [(label,r['t'],r['deposit'][:,np.searchsorted(r['x'],x-1e-9)]) for label,r in records]) for x in points])
def transient(name,a,b,times):
    r=load(a,b)
    chart(name,[('Suspension','Distance x (m)','Normalized concentration c/c0',.4,1,
        [(f'{t} s',r['x'],r['c'][np.searchsorted(r['t'],t)]/.05) for t in times]),
        ('Deposit','Distance x (m)','Deposit fraction ρ',.4,.104,
        [(f'{t} s',r['x'],r['deposit'][np.searchsorted(r['t'],t)]) for t in times])])
def outlet(name,groups):
    chart(name,[(title,'Time (s)','Outlet c/c0',4800,.45,[(label,r['t'],r['c'][:,-1]/.05)
        for label,a,b,tag in cases for r in [load(a,b,tag)]]) for title,cases in groups])
def inventory(name,groups):
    chart(name,[(title,'Time (s)','Integrated deposit (m)',4800,.028,[(label,r['t'],np.trapezoid(r['deposit'],r['x'],axis=1))
        for label,a,b,tag in cases for r in [load(a,b,tag)]]) for title,cases in groups])

kin=[(f'β = {b:g}',2,b,'fine') for b in (1,.9,.8)]
base=[(f'α = {a:g}',a,1,'fine') for a in (2,1.8,1.6)]
joint=[(f'{a:g} / {b:g}',a,b,'fine') for b in (.9,.8) for a in (1.8,1.6)]
profiles('baseline',base)
profiles('kinetic',kin)
histories('near',kin)
transient('kinetic09',2,.9,(1600,3200,4800))
transient('kinetic08',2,.8,(1600,3200,4800))
transient('joint1809',1.8,.9,(450,900,1350))
transient('joint1609',1.6,.9,(450,900,1350))
transient('joint1808',1.8,.8,(450,900,1350))
transient('joint1608',1.6,.8,(450,900,1350))
profiles('joint',joint)
histories('jointnear',joint)
outlet('outlet',[('Kinetic sweep',kin),('Coupled cases',joint)])
inventory('inventory',[('Kinetic sweep',kin),('Coupled cases',joint)])
profiles('variant',[('Fractional balance',1.6,.8,'fine'),('Ordinary balance',1.6,.8,'ordinary'),('Higher threshold',1.6,.8,'threshold')])
profiles('refinement',[(s,1.6,.8,s) for s in ('coarse','medium','fine')])
