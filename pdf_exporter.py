from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas

from alignment_core import (
    ANALYSIS_MODE_CONSERVATIVE,
    ANALYSIS_MODE_MUTATION,
    classify_substitution,
    get_aa_primary_class,
)

try:
    from html_exporter import (
        DISPLAY_STYLE_AA_CLASS,
        DISPLAY_STYLE_ESPRIPT,
        DEFAULT_ESPRIPT_COLORS,
    )
except Exception:
    DISPLAY_STYLE_AA_CLASS = "aa_class"
    DISPLAY_STYLE_ESPRIPT = "espript"
    DEFAULT_ESPRIPT_COLORS = {
        "conserved_bg": "#B85C5C",
        "conserved_text": "#FFFFFF",
        "similar_text": "#9E4F4F",
        "similar_border": "#5D7FA8",
        "gap_bg": "#E9E9E9",
        "gap_text": "#777777",
    }


AA_COLORS = {
    "hydrophobic": colors.Color(1.0, 0.95, 0.64),
    "aromatic": colors.Color(1.0, 0.82, 0.60),
    "positive": colors.Color(0.72, 0.86, 1.0),
    "negative": colors.Color(1.0, 0.75, 0.84),
    "polar": colors.Color(0.78, 1.0, 0.75),
    "special": colors.Color(0.90, 0.90, 0.90),
    "gap": colors.Color(0.82, 0.82, 0.82),
    "unknown": colors.white,
}


MUTATION_COLORS = {
    "identical_column": colors.white,
    "mutation_column": colors.Color(0.784, 0.353, 0.353),
    "gap_column": colors.white,
    "identical": colors.white,
    "mutation": colors.Color(0.784, 0.353, 0.353),
    "gap": colors.white,
}


BORDER_COLORS = {
    "identical": colors.Color(0.10, 0.45, 0.22),
    "strong_conservative": colors.Color(0.25, 0.55, 0.35),
    "weak_similar": colors.grey,
    "neutral": colors.grey,
    "identical_column": colors.Color(0.10, 0.45, 0.22),
    "variable_column": None,
    "mutation": None,
    "mutation_column": None,
    "gap": colors.grey,
    "gap_column": colors.grey,
    "espript_conserved": None,
    "espript_strong_similar": None,
    "espript_weak_similar": None,
    "espript_nonconserved": None,
    "espript_gap": colors.grey,
}


ESPRIPT_GROUPS = [
    set("AVLIM"),
    set("FWY"),
    set("KRH"),
    set("DE"),
    set("STNQ"),
    set("NQ"),
    set("ST"),
    set("AGS"),
    set("GAS"),
    set("GP"),
    set("C"),
]


def hex_to_color(hex_color, default=colors.white):
    if not hex_color:
        return default

    value = hex_color.strip().lstrip("#")

    if len(value) != 6:
        return default

    try:
        r = int(value[0:2], 16) / 255
        g = int(value[2:4], 16) / 255
        b = int(value[4:6], 16) / 255
        return colors.Color(r, g, b)
    except Exception:
        return default


def merged_espript_colors(color_config=None):
    merged = DEFAULT_ESPRIPT_COLORS.copy()
    if color_config:
        merged.update(color_config)

    return {
        "conserved_bg": hex_to_color(merged["conserved_bg"], colors.Color(0.72, 0.36, 0.36)),
        "conserved_text": hex_to_color(merged["conserved_text"], colors.white),
        "similar_text": hex_to_color(merged["similar_text"], colors.Color(0.62, 0.31, 0.31)),
        "similar_border": hex_to_color(merged["similar_border"], colors.Color(0.36, 0.50, 0.66)),
        "gap_bg": hex_to_color(merged["gap_bg"], colors.Color(0.91, 0.91, 0.91)),
        "gap_text": hex_to_color(merged["gap_text"], colors.grey),
    }


def get_column_classes_for_msa(records, analysis_mode):
    if not records:
        return []

    aln_len = len(records[0]["sequence"])
    column_classes = []

    for col_index in range(aln_len):
        column = [record["sequence"][col_index] for record in records]

        if "-" in column:
            column_classes.append("gap_column")
        elif len(set(column)) == 1:
            column_classes.append("identical_column")
        else:
            if analysis_mode == ANALYSIS_MODE_MUTATION:
                column_classes.append("mutation_column")
            else:
                column_classes.append("variable_column")

    return column_classes


def espript_column_classes(records):
    if not records:
        return []

    aln_len = len(records[0]["sequence"])
    classes = []

    for col_index in range(aln_len):
        column = [record["sequence"][col_index] for record in records]
        non_gap = [aa for aa in column if aa != "-"]

        if not non_gap or len(non_gap) != len(column):
            classes.append("espript_gap")
            continue

        unique = set(non_gap)

        if len(unique) == 1:
            classes.append("espript_conserved")
            continue

        if residues_share_group(non_gap):
            classes.append("espript_strong_similar")
            continue

        if majority_share_group(non_gap, cutoff=0.70):
            classes.append("espript_weak_similar")
            continue

        classes.append("espript_nonconserved")

    return classes


def residues_share_group(residues):
    residue_set = set(residues)
    for group in ESPRIPT_GROUPS:
        if residue_set.issubset(group):
            return True
    return False


def majority_share_group(residues, cutoff=0.70):
    n = len(residues)
    if n == 0:
        return False

    for group in ESPRIPT_GROUPS:
        count = sum(1 for aa in residues if aa in group)
        if count / n >= cutoff:
            return True

    return False


def get_runs(classes):
    if not classes:
        return []

    runs = []
    start = 0
    current = classes[0]

    for i in range(1, len(classes)):
        if classes[i] != current:
            runs.append((start, i, current))
            start = i
            current = classes[i]

    runs.append((start, len(classes), current))
    return runs


def draw_header(c, title, subtitle, page_num, width, height):
    c.setFillColor(colors.black)
    c.setFont("Helvetica-Bold", 17)
    c.drawString(16 * mm, height - 14 * mm, title)

    c.setFont("Helvetica", 10)
    c.setFillColor(colors.Color(0.30, 0.30, 0.30))
    c.drawString(16 * mm, height - 20 * mm, subtitle)

    c.setFont("Helvetica", 9)
    c.drawRightString(width - 16 * mm, height - 14 * mm, f"Page {page_num}")


def draw_summary_box(c, x, y, width, lines):
    line_h = 5.0 * mm
    box_h = (len(lines) + 1) * line_h

    c.setFillColor(colors.Color(0.97, 0.97, 0.97))
    c.setStrokeColor(colors.Color(0.82, 0.82, 0.82))
    c.roundRect(x, y - box_h, width, box_h, 3 * mm, fill=1, stroke=1)

    c.setFillColor(colors.black)
    c.setFont("Helvetica", 9)

    text_y = y - 4.5 * mm
    for line in lines:
        c.drawString(x + 4 * mm, text_y, line)
        text_y -= line_h

    return y - box_h - 6 * mm


def draw_legend_conservative(c, x, y):
    c.setFont("Helvetica-Bold", 9)
    c.setFillColor(colors.black)
    c.drawString(x, y, "Legend: amino acid classes")
    y -= 5 * mm

    legend_items = [
        ("A V L I M", "Hydrophobic", "hydrophobic"),
        ("F W Y", "Aromatic", "aromatic"),
        ("K R H", "Positive", "positive"),
        ("D E", "Negative", "negative"),
        ("S T N Q", "Polar", "polar"),
        ("G P C", "Special", "special"),
        ("-", "Gap", "gap"),
    ]

    c.setFont("Helvetica", 8)

    cursor_x = x
    for text, label, aa_class in legend_items:
        bg = AA_COLORS.get(aa_class, colors.white)
        c.setFillColor(bg)
        c.rect(cursor_x, y - 3.2 * mm, 18 * mm, 4.2 * mm, fill=1, stroke=0)

        c.setFillColor(colors.black)
        c.setFont("Courier-Bold", 7.5)
        c.drawCentredString(cursor_x + 9 * mm, y - 2.1 * mm, text)

        c.setFont("Helvetica", 8)
        c.drawString(cursor_x + 20 * mm, y - 2.1 * mm, label)

        cursor_x += 47 * mm

    return y - 10 * mm


def draw_legend_mutation(c, x, y):
    c.setFont("Helvetica-Bold", 9)
    c.setFillColor(colors.black)
    c.drawString(x, y, "Legend: mutation site comparison")
    y -= 5 * mm

    legend_items = [
        ("A", "Same / unchanged", "identical_column"),
        ("A", "Different / mutation site", "mutation_column"),
        ("-", "Gap / insertion-deletion", "gap_column"),
    ]

    c.setFont("Helvetica", 8)

    cursor_x = x
    for text, label, cls in legend_items:
        bg = MUTATION_COLORS.get(cls, colors.white)

        c.setFillColor(bg)
        c.rect(cursor_x, y - 3.2 * mm, 7 * mm, 4.2 * mm, fill=1, stroke=0)

        if cls == "mutation_column":
            c.setFillColor(colors.white)
            c.setFont("Courier-Bold", 8)
        else:
            c.setFillColor(colors.black)
            c.setFont("Courier", 8)

        c.drawCentredString(cursor_x + 3.5 * mm, y - 2.1 * mm, text)

        c.setFillColor(colors.black)
        c.setFont("Helvetica", 8)
        c.drawString(cursor_x + 9 * mm, y - 2.1 * mm, label)

        cursor_x += 58 * mm

    return y - 10 * mm


def draw_legend_espript(c, x, y, espript_colors):
    c.setFont("Helvetica-Bold", 9)
    c.setFillColor(colors.black)
    c.drawString(x, y, "Legend: ESPript-like conservation style")
    y -= 5 * mm

    legend_items = [
        ("A", "Fully conserved", "espript_conserved"),
        ("A", "Strongly similar", "espript_strong_similar"),
        ("A", "Weakly similar", "espript_weak_similar"),
        ("A", "Non-conserved", "espript_nonconserved"),
        (".", "Gap-containing", "espript_gap"),
    ]

    cursor_x = x
    for text, label, cls in legend_items:
        bg, text_color, border_color, bold = get_espript_cell_style(cls, espript_colors)

        c.setFillColor(bg)
        c.rect(cursor_x, y - 3.2 * mm, 7 * mm, 4.2 * mm, fill=1, stroke=0)

        if border_color is not None:
            c.setStrokeColor(border_color)
            c.setLineWidth(0.7)
            c.rect(cursor_x, y - 3.2 * mm, 7 * mm, 4.2 * mm, fill=0, stroke=1)

        c.setFillColor(text_color)
        c.setFont("Courier-Bold" if bold else "Courier", 8)
        c.drawCentredString(cursor_x + 3.5 * mm, y - 2.1 * mm, text)

        c.setFillColor(colors.black)
        c.setFont("Helvetica", 8)
        c.drawString(cursor_x + 9 * mm, y - 2.1 * mm, label)

        cursor_x += 45 * mm

    return y - 10 * mm


def get_espript_cell_style(cls, espript_colors):
    if cls == "espript_conserved":
        return (
            espript_colors["conserved_bg"],
            espript_colors["conserved_text"],
            None,
            True,
        )

    if cls == "espript_strong_similar":
        return (
            colors.white,
            espript_colors["similar_text"],
            espript_colors["similar_border"],
            True,
        )

    if cls == "espript_weak_similar":
        return (
            colors.white,
            espript_colors["similar_text"],
            None,
            True,
        )

    if cls == "espript_gap":
        return (
            espript_colors["gap_bg"],
            espript_colors["gap_text"],
            colors.grey,
            False,
        )

    return (
        colors.white,
        colors.black,
        None,
        False,
    )


def draw_residue_cell(c, x, y, aa, bg_color, cell_w, cell_h, text_color=colors.black, bold=False):
    c.setFillColor(bg_color)
    c.rect(x, y, cell_w, cell_h, fill=1, stroke=0)

    c.setFillColor(text_color)
    c.setFont("Courier-Bold" if bold else "Courier", 10)
    c.drawCentredString(x + cell_w / 2, y + 1.25 * mm, aa)


def draw_run_borders(c, x0, y, classes, cell_w, cell_h, color_overrides=None):
    runs = get_runs(classes)

    for start, end, cls in runs:
        border = None

        if color_overrides and cls in color_overrides:
            border = color_overrides[cls]
        else:
            border = BORDER_COLORS.get(cls)

        if border is None:
            continue

        c.setStrokeColor(border)
        c.setLineWidth(0.75)
        c.rect(
            x0 + start * cell_w,
            y,
            (end - start) * cell_w,
            cell_h,
            fill=0,
            stroke=1,
        )


def export_pairwise_pdf(
    result,
    file_path,
    title="Pairwise Protein Sequence Alignment",
    seq_a_name="SeqA",
    seq_b_name="SeqB",
    line_width=85,
):
    analysis_mode = result.get("analysis_mode", ANALYSIS_MODE_CONSERVATIVE)

    page_size = landscape(A3)
    width, height = page_size

    c = canvas.Canvas(file_path, pagesize=page_size)

    left = 16 * mm
    top = height - 28 * mm
    bottom = 16 * mm

    label_w = 42 * mm
    cell_w = 4.0 * mm
    cell_h = 5.0 * mm
    row_gap = 1.2 * mm
    block_gap = 8 * mm

    page_num = 1

    subtitle = (
        "Analysis mode: Mutation site comparison"
        if analysis_mode == ANALYSIS_MODE_MUTATION
        else "Analysis mode: Conservative analysis"
    )

    draw_header(c, title, subtitle, page_num, width, height)

    y = top

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        summary_lines = [
            f"Identity: {result['identity']:.2f}%",
            f"Mutation sites: {result['mutation_percent']:.2f}%",
            f"Gaps: {result['gap_percent']:.2f}%",
            f"Alignment score: {result['score']:.2f}",
        ]
        y = draw_summary_box(c, left, y, 95 * mm, summary_lines)
        y = draw_legend_mutation(c, left, y)
    else:
        summary_lines = [
            f"Identity: {result['identity']:.2f}%",
            f"Strong conservative: {result['strong_conservative']:.2f}%",
            f"Weak / neutral similarity: {result['weak_similar']:.2f}%",
            f"Total similarity: {result['similarity']:.2f}%",
            f"Non-conservative: {result['non_conservative']:.2f}%",
            f"Gaps: {result['gap_percent']:.2f}%",
            f"Alignment score: {result['score']:.2f}",
        ]
        y = draw_summary_box(c, left, y, 110 * mm, summary_lines)
        y = draw_legend_conservative(c, left, y)

    aligned_a = result["aligned_seq_a"]
    aligned_b = result["aligned_seq_b"]

    for start in range(0, len(aligned_a), line_width):
        block_a = aligned_a[start:start + line_width]
        block_b = aligned_b[start:start + line_width]

        block_classes = []
        for a, b in zip(block_a, block_b):
            info = classify_substitution(a, b, analysis_mode=analysis_mode)
            block_classes.append(info["type"])

        needed_height = 3 * cell_h + 2 * row_gap + block_gap + 6 * mm

        if y - needed_height < bottom:
            c.showPage()
            page_num += 1
            draw_header(c, title, subtitle, page_num, width, height)
            y = top

        c.setFont("Helvetica", 8)
        c.setFillColor(colors.grey)
        c.drawString(left, y, f"Position {start + 1}-{start + len(block_a)}")
        y -= 5 * mm

        rows = [
            (seq_a_name, block_a),
            ("", result["match_line"][start:start + line_width]),
            (seq_b_name, block_b),
        ]

        for row_index, (label, seq) in enumerate(rows):
            c.setFont("Helvetica-Bold", 9)
            c.setFillColor(colors.black)
            c.drawString(left, y + 1.6 * mm, label[:26])

            x = left + label_w

            if row_index == 1:
                for offset, symbol in enumerate(seq):
                    cls = block_classes[offset]
                    c.setFont("Courier-Bold", 10)

                    if analysis_mode == ANALYSIS_MODE_MUTATION:
                        if cls == "mutation":
                            c.setFillColor(colors.Color(0.78, 0.00, 0.00))
                        elif cls == "gap":
                            c.setFillColor(colors.grey)
                        else:
                            c.setFillColor(colors.black)
                    else:
                        if cls == "identical":
                            c.setFillColor(colors.Color(0.00, 0.35, 0.00))
                        elif cls == "strong_conservative":
                            c.setFillColor(colors.Color(0.10, 0.35, 0.18))
                        else:
                            c.setFillColor(colors.grey)

                    c.drawCentredString(x + cell_w / 2, y + 1.25 * mm, symbol)
                    x += cell_w
            else:
                for offset, aa in enumerate(seq):
                    cls = block_classes[offset]

                    if analysis_mode == ANALYSIS_MODE_MUTATION:
                        bg = MUTATION_COLORS.get(cls, colors.white)
                        text_color = colors.white if cls == "mutation" else colors.black
                        bold = cls in ["mutation", "identical"]
                    else:
                        bg = AA_COLORS.get(get_aa_primary_class(aa), colors.white)
                        text_color = colors.black
                        bold = cls in ["identical", "identical_column"]

                    draw_residue_cell(
                        c,
                        x,
                        y,
                        aa,
                        bg,
                        cell_w,
                        cell_h,
                        text_color=text_color,
                        bold=bold,
                    )
                    x += cell_w

                if analysis_mode != ANALYSIS_MODE_MUTATION:
                    draw_run_borders(c, left + label_w, y, block_classes, cell_w, cell_h)

            y -= cell_h + row_gap

        y -= block_gap

    c.save()


def export_msa_pdf(
    records,
    file_path,
    title="Multiple Sequence Alignment",
    line_width=85,
    analysis_mode=ANALYSIS_MODE_CONSERVATIVE,
    display_style=DISPLAY_STYLE_AA_CLASS,
    color_config=None,
):
    if not records:
        raise ValueError("No MSA records to export.")

    page_size = landscape(A3)
    width, height = page_size

    c = canvas.Canvas(file_path, pagesize=page_size)

    left = 16 * mm
    top = height - 28 * mm
    bottom = 16 * mm

    max_name_len = max(len(r.get("display_name", r.get("id", ""))) for r in records)
    label_w = min(max(42 * mm, max_name_len * 2.0 * mm), 70 * mm)

    cell_w = 4.0 * mm
    cell_h = 5.0 * mm
    row_gap = 1.0 * mm
    block_gap = 8 * mm

    aln_len = len(records[0]["sequence"])

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        column_classes = get_column_classes_for_msa(records, analysis_mode)
        effective_style = "Mutation site comparison"
    elif display_style == DISPLAY_STYLE_ESPRIPT:
        column_classes = espript_column_classes(records)
        effective_style = "Conservative analysis | ESPript-like conservation style"
    else:
        column_classes = get_column_classes_for_msa(records, analysis_mode)
        effective_style = "Conservative analysis | Amino acid class style"

    espript_colors = merged_espript_colors(color_config)

    page_num = 1
    draw_header(c, title, effective_style, page_num, width, height)

    y = top

    summary_lines = [
        f"Number of sequences: {len(records)}",
        f"Alignment length: {aln_len}",
    ]
    y = draw_summary_box(c, left, y, 95 * mm, summary_lines)

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        y = draw_legend_mutation(c, left, y)
    elif display_style == DISPLAY_STYLE_ESPRIPT:
        y = draw_legend_espript(c, left, y, espript_colors)
    else:
        y = draw_legend_conservative(c, left, y)

    for start in range(0, aln_len, line_width):
        end = min(start + line_width, aln_len)
        block_classes = column_classes[start:end]

        needed_height = len(records) * (cell_h + row_gap) + block_gap + 6 * mm

        if y - needed_height < bottom:
            c.showPage()
            page_num += 1
            draw_header(c, title, effective_style, page_num, width, height)
            y = top

        c.setFont("Helvetica", 8)
        c.setFillColor(colors.grey)
        c.drawString(left, y, f"Position {start + 1}-{end}")
        y -= 5 * mm

        for record in records:
            name = record.get("display_name", record.get("id", "sequence"))
            seq = record["sequence"][start:end]

            c.setFont("Helvetica-BoldOblique", 10)
            c.setFillColor(colors.black)
            c.drawString(left, y + 1.2 * mm, name[:32])

            x = left + label_w

            for offset, aa in enumerate(seq):
                cls = block_classes[offset]

                if analysis_mode == ANALYSIS_MODE_MUTATION:
                    bg = MUTATION_COLORS.get(cls, colors.white)
                    text_color = colors.white if cls == "mutation_column" else colors.black
                    bold = cls in ["mutation_column", "identical_column"]
                    display_aa = aa
                elif display_style == DISPLAY_STYLE_ESPRIPT:
                    bg, text_color, border_color, bold = get_espript_cell_style(cls, espript_colors)
                    display_aa = "." if aa == "-" else aa
                else:
                    bg = AA_COLORS.get(get_aa_primary_class(aa), colors.white)
                    text_color = colors.black
                    bold = cls == "identical_column"
                    display_aa = aa

                draw_residue_cell(
                    c,
                    x,
                    y,
                    display_aa,
                    bg,
                    cell_w,
                    cell_h,
                    text_color=text_color,
                    bold=bold,
                )
                x += cell_w

            if analysis_mode == ANALYSIS_MODE_MUTATION:
                pass
            elif display_style == DISPLAY_STYLE_ESPRIPT:
                draw_run_borders(
                    c,
                    left + label_w,
                    y,
                    block_classes,
                    cell_w,
                    cell_h,
                    color_overrides={
                        "espript_conserved": espript_colors["conserved_bg"],
                        "espript_strong_similar": espript_colors["similar_border"],
                        "espript_gap": colors.grey,
                    },
                )
            else:
                draw_run_borders(c, left + label_w, y, block_classes, cell_w, cell_h)

            y -= cell_h + row_gap

        y -= block_gap

    c.save()
