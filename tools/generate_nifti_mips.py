import gzip
import struct
import zlib
from pathlib import Path


PX_PER_MM = 3.0
MARGIN = 12
GAP = 16
LABEL_PAD = 20

BG = 18
BORDER = 96
TEXT = 210
ON = 255

FONT_5X7 = {
    "X": [
        "10001",
        "01010",
        "00100",
        "00100",
        "00100",
        "01010",
        "10001",
    ],
    "Y": [
        "10001",
        "01010",
        "00100",
        "00100",
        "00100",
        "00100",
        "00100",
    ],
    "Z": [
        "11111",
        "00010",
        "00100",
        "01000",
        "10000",
        "10000",
        "11111",
    ],
}


def png_chunk(tag, data):
    crc = zlib.crc32(tag + data) & 0xFFFFFFFF
    return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)


def write_grayscale_png(path, width, height, pixels):
    raw = bytearray()
    for y in range(height):
        raw.append(0)
        start = y * width
        raw.extend(pixels[start:start + width])

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0)
    idat = zlib.compress(bytes(raw), level=9)

    png = bytearray(b"\x89PNG\r\n\x1a\n")
    png.extend(png_chunk(b"IHDR", ihdr))
    png.extend(png_chunk(b"IDAT", idat))
    png.extend(png_chunk(b"IEND", b""))
    path.write_bytes(png)


def read_nifti_u8(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as f:
        hdr = f.read(348)
        sizeof_hdr = struct.unpack("<I", hdr[0:4])[0]
        if sizeof_hdr != 348:
            raise ValueError(f"{path}: unexpected header size {sizeof_hdr}")

        dim = struct.unpack("<8h", hdr[40:56])
        datatype = struct.unpack("<h", hdr[70:72])[0]
        bitpix = struct.unpack("<h", hdr[72:74])[0]
        pixdim = struct.unpack("<8f", hdr[76:108])
        vox_offset = int(struct.unpack("<f", hdr[108:112])[0])

        shape = tuple(dim[1:1 + dim[0]])
        spacing = tuple(pixdim[1:1 + dim[0]])
        if datatype != 2 or bitpix != 8:
            raise ValueError(
                f"{path}: only uint8 masks supported, got datatype={datatype}, bitpix={bitpix}"
            )

        f.seek(vox_offset)
        raw = f.read(shape[0] * shape[1] * shape[2])
    return shape, spacing, raw


def project_mask(shape, raw):
    xdim, ydim, zdim = shape
    slice_size = xdim * ydim

    xy = bytearray(xdim * ydim)
    xz = bytearray(xdim * zdim)
    yz = bytearray(ydim * zdim)

    min_x = min_y = min_z = None
    max_x = max_y = max_z = None

    for z in range(zdim):
        slab = raw[z * slice_size:(z + 1) * slice_size]
        xz_row = z * xdim
        yz_row = z * ydim

        for offset, value in enumerate(slab):
            if not value:
                continue

            x = offset % xdim
            y = offset // xdim

            xy[y * xdim + x] = ON
            xz[xz_row + x] = ON
            yz[yz_row + y] = ON

            if min_x is None or x < min_x:
                min_x = x
            if max_x is None or x > max_x:
                max_x = x
            if min_y is None or y < min_y:
                min_y = y
            if max_y is None or y > max_y:
                max_y = y
            if min_z is None or z < min_z:
                min_z = z
            if max_z is None or z > max_z:
                max_z = z

    bbox = (min_x, min_y, min_z, max_x, max_y, max_z)
    return xy, xz, yz, bbox


def crop_plane(data, src_w, box):
    x0, y0, x1, y1 = box
    out_w = x1 - x0
    out_h = y1 - y0
    cropped = bytearray(out_w * out_h)

    for row in range(out_h):
        src_start = (y0 + row) * src_w + x0
        dst_start = row * out_w
        cropped[dst_start:dst_start + out_w] = data[src_start:src_start + out_w]

    return out_w, out_h, cropped


def resize_nn(data, src_w, src_h, dst_w, dst_h):
    if src_w == dst_w and src_h == dst_h:
        return bytearray(data)

    xmap = [min((x * src_w) // dst_w, src_w - 1) for x in range(dst_w)]
    ymap = [min((y * src_h) // dst_h, src_h - 1) for y in range(dst_h)]

    out = bytearray(dst_w * dst_h)
    for y in range(dst_h):
        src_row = ymap[y] * src_w
        dst_row = y * dst_w
        for x in range(dst_w):
            out[dst_row + x] = data[src_row + xmap[x]]
    return out


def mm_to_px(length_voxels, spacing):
    return max(1, round(length_voxels * spacing * PX_PER_MM))


def draw_rect(buf, canvas_w, canvas_h, x, y, w, h, color):
    x2 = min(x + w - 1, canvas_w - 1)
    y2 = min(y + h - 1, canvas_h - 1)
    if x < 0 or y < 0 or x >= canvas_w or y >= canvas_h:
        return

    for xx in range(x, x2 + 1):
        buf[y * canvas_w + xx] = color
        buf[y2 * canvas_w + xx] = color
    for yy in range(y, y2 + 1):
        buf[yy * canvas_w + x] = color
        buf[yy * canvas_w + x2] = color


def blit(buf, canvas_w, x, y, img_w, img_h, img):
    for row in range(img_h):
        dst = (y + row) * canvas_w + x
        src = row * img_w
        buf[dst:dst + img_w] = img[src:src + img_w]


def draw_char(buf, canvas_w, canvas_h, x, y, ch, scale=2, color=TEXT):
    pattern = FONT_5X7.get(ch)
    if not pattern:
        return
    for py, row in enumerate(pattern):
        for px, bit in enumerate(row):
            if bit != "1":
                continue
            for dy in range(scale):
                for dx in range(scale):
                    xx = x + px * scale + dx
                    yy = y + py * scale + dy
                    if 0 <= xx < canvas_w and 0 <= yy < canvas_h:
                        buf[yy * canvas_w + xx] = color


def draw_text(buf, canvas_w, canvas_h, x, y, text, scale=2, color=TEXT):
    cursor = x
    for ch in text:
        if ch == " ":
            cursor += 4 * scale
            continue
        draw_char(buf, canvas_w, canvas_h, cursor, y, ch, scale=scale, color=color)
        cursor += 6 * scale


def render_triptych(xy_img, xy_size, xz_img, xz_size, yz_img, yz_size):
    panels = [
        ("XY", xy_img, xy_size[0], xy_size[1]),
        ("XZ", xz_img, xz_size[0], xz_size[1]),
        ("YZ", yz_img, yz_size[0], yz_size[1]),
    ]
    max_h = max(h for _, _, _, h in panels)
    canvas_w = GAP * (len(panels) + 1) + sum(w for _, _, w, _ in panels)
    canvas_h = GAP * 2 + LABEL_PAD + max_h

    canvas = bytearray([BG] * (canvas_w * canvas_h))

    x = GAP
    for label, img, w, h in panels:
        y = GAP + LABEL_PAD + (max_h - h) // 2
        blit(canvas, canvas_w, x, y, w, h, img)
        draw_rect(canvas, canvas_w, canvas_h, x - 1, y - 1, w + 2, h + 2, BORDER)
        draw_text(canvas, canvas_w, canvas_h, x, GAP, label, scale=2, color=TEXT)
        x += w + GAP

    return canvas_w, canvas_h, canvas


def crop_bounds(bbox, axis_a, axis_b, dims):
    starts = {
        "x": bbox[0],
        "y": bbox[1],
        "z": bbox[2],
    }
    ends = {
        "x": bbox[3],
        "y": bbox[4],
        "z": bbox[5],
    }
    limits = {
        "x": dims[0],
        "y": dims[1],
        "z": dims[2],
    }
    a0 = max(starts[axis_a] - MARGIN, 0)
    b0 = max(starts[axis_b] - MARGIN, 0)
    a1 = min(ends[axis_a] + MARGIN + 1, limits[axis_a])
    b1 = min(ends[axis_b] + MARGIN + 1, limits[axis_b])
    return a0, b0, a1, b1


def stem_from_path(path):
    parts = path.parts
    return f"{parts[-3]}_{parts[-2]}"


def main():
    root = Path.cwd()
    inputs = sorted(root.glob("case*/**/mask.nii.gz"))
    if not inputs:
        raise SystemExit("No mask.nii.gz files found under case*/")

    out_dir = root / "nifti_mips"
    out_dir.mkdir(exist_ok=True)

    for path in inputs:
        shape, spacing, raw = read_nifti_u8(path)
        xy, xz, yz, bbox = project_mask(shape, raw)
        if bbox[0] is None:
            raise ValueError(f"{path}: empty mask")

        xy_box = crop_bounds(bbox, "x", "y", shape)
        xz_box = crop_bounds(bbox, "x", "z", shape)
        yz_box = crop_bounds(bbox, "y", "z", shape)

        xy_w, xy_h, xy_crop = crop_plane(xy, shape[0], xy_box)
        xz_w, xz_h, xz_crop = crop_plane(xz, shape[0], xz_box)
        yz_w, yz_h, yz_crop = crop_plane(yz, shape[1], yz_box)

        xy_img = resize_nn(
            xy_crop,
            xy_w,
            xy_h,
            mm_to_px(xy_w, spacing[0]),
            mm_to_px(xy_h, spacing[1]),
        )
        xz_img = resize_nn(
            xz_crop,
            xz_w,
            xz_h,
            mm_to_px(xz_w, spacing[0]),
            mm_to_px(xz_h, spacing[2]),
        )
        yz_img = resize_nn(
            yz_crop,
            yz_w,
            yz_h,
            mm_to_px(yz_w, spacing[1]),
            mm_to_px(yz_h, spacing[2]),
        )

        xy_size = (mm_to_px(xy_w, spacing[0]), mm_to_px(xy_h, spacing[1]))
        xz_size = (mm_to_px(xz_w, spacing[0]), mm_to_px(xz_h, spacing[2]))
        yz_size = (mm_to_px(yz_w, spacing[1]), mm_to_px(yz_h, spacing[2]))

        canvas_w, canvas_h, canvas = render_triptych(
            xy_img, xy_size, xz_img, xz_size, yz_img, yz_size
        )

        out_path = out_dir / f"{stem_from_path(path)}_mip_triptych.png"
        write_grayscale_png(out_path, canvas_w, canvas_h, canvas)

        print(
            f"{out_path.name}\tshape={shape}\tspacing="
            f"{tuple(round(v, 6) for v in spacing)}\tbbox={bbox}"
        )


if __name__ == "__main__":
    main()
