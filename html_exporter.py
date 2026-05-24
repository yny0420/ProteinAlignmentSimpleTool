import html

from alignment_core import (
    ANALYSIS_MODE_CONSERVATIVE,
    ANALYSIS_MODE_MUTATION,
    classify_substitution,
    get_aa_primary_class,
)

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


def merged_colors(color_config=None):
    colors = DEFAULT_ESPRIPT_COLORS.copy()
    if color_config:
        colors.update(color_config)
    return colors


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


def get_run_edges(classes):
    edges = []

    for i, cls in enumerate(classes):
        prev_cls = classes[i - 1] if i > 0 else None
        next_cls = classes[i + 1] if i < len(classes) - 1 else None

        edge = []

        if cls != prev_cls:
            edge.append("run_start")

        if cls != next_cls:
            edge.append("run_end")

        edges.append(edge)

    return edges


def aa_span_conservative(aa: str, column_class: str = "", edge_classes=None) -> str:
    aa_class = get_aa_primary_class(aa)
    classes = ["aa", aa_class]

    if column_class:
        classes.append(column_class)

    if edge_classes:
        classes.extend(edge_classes)

    class_text = " ".join(classes)
    return f'<span class="{class_text}">{html.escape(aa)}</span>'


def aa_span_mutation(aa: str, column_class: str = "", edge_classes=None) -> str:
    classes = ["aa", "mutation_plain"]

    if aa == "-":
        classes.append("gap")

    if column_class:
        classes.append(column_class)

    if edge_classes:
        classes.extend(edge_classes)

    display_aa = "-" if aa == "-" else aa
    class_text = " ".join(classes)
    return f'<span class="{class_text}">{html.escape(display_aa)}</span>'


def aa_span_espript(aa: str, column_class: str = "", edge_classes=None) -> str:
    classes = ["aa", "espript"]

    if column_class:
        classes.append(column_class)

    if edge_classes:
        classes.extend(edge_classes)

    display_aa = "." if aa == "-" else aa
    class_text = " ".join(classes)
    return f'<span class="{class_text}">{html.escape(display_aa)}</span>'


def match_span(symbol: str, sub_type: str) -> str:
    safe_symbol = html.escape(symbol if symbol != " " else "\u00A0")
    return f'<span class="match {sub_type}">{safe_symbol}</span>'


def pairwise_result_to_html(
    result: dict,
    title="Pairwise Protein Sequence Alignment",
    seq_a_name="SeqA",
    seq_b_name="SeqB",
    line_width=80,
) -> str:
    analysis_mode = result.get("analysis_mode", ANALYSIS_MODE_CONSERVATIVE)

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        return pairwise_result_to_html_mutation(
            result,
            title=title,
            seq_a_name=seq_a_name,
            seq_b_name=seq_b_name,
            line_width=line_width,
        )

    return pairwise_result_to_html_conservative(
        result,
        title=title,
        seq_a_name=seq_a_name,
        seq_b_name=seq_b_name,
        line_width=line_width,
    )


def pairwise_result_to_html_conservative(
    result: dict,
    title="Pairwise Protein Sequence Alignment",
    seq_a_name="SeqA",
    seq_b_name="SeqB",
    line_width=80,
) -> str:
    aligned_a = result["aligned_seq_a"]
    aligned_b = result["aligned_seq_b"]
    counts = result.get("counts", {})

    blocks = []

    for start in range(0, len(aligned_a), line_width):
        block_a = aligned_a[start:start + line_width]
        block_b = aligned_b[start:start + line_width]

        col_classes = []
        for a, b in zip(block_a, block_b):
            info = classify_substitution(a, b, analysis_mode=ANALYSIS_MODE_CONSERVATIVE)
            col_classes.append(info["type"])

        edge_classes = get_run_edges(col_classes)

        seq_a_html = []
        seq_b_html = []
        match_html = []

        for i, (a, b) in enumerate(zip(block_a, block_b)):
            info = classify_substitution(a, b, analysis_mode=ANALYSIS_MODE_CONSERVATIVE)
            sub_type = info["type"]

            seq_a_html.append(aa_span_conservative(a, sub_type, edge_classes[i]))
            seq_b_html.append(aa_span_conservative(b, sub_type, edge_classes[i]))
            match_html.append(match_span(info["symbol"], sub_type))

        blocks.append(
            f"""
<div class="alignment-block">
<div class="position">Position {start + 1}-{start + len(block_a)}</div>
<div><span class="label">{html.escape(seq_a_name)}</span> <span class="seq">{''.join(seq_a_html)}</span></div>
<div><span class="label"></span> <span class="seq">{''.join(match_html)}</span></div>
<div><span class="label">{html.escape(seq_b_name)}</span> <span class="seq">{''.join(seq_b_html)}</span></div>
</div>
"""
        )

    body = "\n".join(blocks)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{html.escape(title)}</title>
<style>
{base_css()}
{conservative_css()}
</style>
</head>
<body>

<h1>{html.escape(title)}</h1>
<h2>Analysis mode: Conservative analysis</h2>

<div class="summary">
<div class="summary-grid">
<div><strong>Identity:</strong></div><div>{result["identity"]:.2f}%</div>
<div><strong>Strong conservative:</strong></div><div>{result["strong_conservative"]:.2f}%</div>
<div><strong>Weak / neutral similarity:</strong></div><div>{result["weak_similar"]:.2f}%</div>
<div><strong>Total similarity:</strong></div><div>{result["similarity"]:.2f}%</div>
<div><strong>Non-conservative:</strong></div><div>{result["non_conservative"]:.2f}%</div>
<div><strong>Gaps:</strong></div><div>{result["gap_percent"]:.2f}%</div>
<div><strong>Alignment score:</strong></div><div>{result["score"]:.2f}</div>
</div>

<br>

<div class="summary-grid">
<div>Identical count:</div><div>{counts.get("identical", 0)}</div>
<div>Strong conservative count:</div><div>{counts.get("strong_conservative", 0)}</div>
<div>Weak similar count:</div><div>{counts.get("weak_similar", 0)}</div>
<div>Neutral count:</div><div>{counts.get("neutral", 0)}</div>
<div>Non-conservative count:</div><div>{counts.get("non_conservative", 0)}</div>
<div>Gap count:</div><div>{counts.get("gap", 0)}</div>
</div>
</div>

{body}

{conservative_legend_html()}

</body>
</html>
"""


def pairwise_result_to_html_mutation(
    result: dict,
    title="Pairwise Protein Sequence Alignment",
    seq_a_name="SeqA",
    seq_b_name="SeqB",
    line_width=80,
) -> str:
    aligned_a = result["aligned_seq_a"]
    aligned_b = result["aligned_seq_b"]
    counts = result.get("counts", {})

    blocks = []

    for start in range(0, len(aligned_a), line_width):
        block_a = aligned_a[start:start + line_width]
        block_b = aligned_b[start:start + line_width]

        col_classes = []
        for a, b in zip(block_a, block_b):
            info = classify_substitution(a, b, analysis_mode=ANALYSIS_MODE_MUTATION)
            col_classes.append(info["type"])

        edge_classes = get_run_edges(col_classes)

        seq_a_html = []
        seq_b_html = []
        match_html = []

        for i, (a, b) in enumerate(zip(block_a, block_b)):
            info = classify_substitution(a, b, analysis_mode=ANALYSIS_MODE_MUTATION)
            sub_type = info["type"]

            seq_a_html.append(aa_span_mutation(a, sub_type, edge_classes[i]))
            seq_b_html.append(aa_span_mutation(b, sub_type, edge_classes[i]))
            match_html.append(match_span(info["symbol"], sub_type))

        blocks.append(
            f"""
<div class="alignment-block mutation-block">
<div class="position">Position {start + 1}-{start + len(block_a)}</div>
<div><span class="label">{html.escape(seq_a_name)}</span> <span class="seq">{''.join(seq_a_html)}</span></div>
<div><span class="label"></span> <span class="seq">{''.join(match_html)}</span></div>
<div><span class="label">{html.escape(seq_b_name)}</span> <span class="seq">{''.join(seq_b_html)}</span></div>
</div>
"""
        )

    body = "\n".join(blocks)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{html.escape(title)}</title>
<style>
{base_css()}
{mutation_css()}
</style>
</head>
<body>

<h1>{html.escape(title)}</h1>
<h2>Analysis mode: Mutation site comparison</h2>

<div class="summary">
<div class="summary-grid">
<div><strong>Identity:</strong></div><div>{result["identity"]:.2f}%</div>
<div><strong>Mutation sites:</strong></div><div>{result["mutation_percent"]:.2f}%</div>
<div><strong>Gaps:</strong></div><div>{result["gap_percent"]:.2f}%</div>
<div><strong>Alignment score:</strong></div><div>{result["score"]:.2f}</div>
</div>

<br>

<div class="summary-grid">
<div>Identical count:</div><div>{counts.get("identical", 0)}</div>
<div>Mutation count:</div><div>{counts.get("mutation", 0)}</div>
<div>Gap count:</div><div>{counts.get("gap", 0)}</div>
</div>
</div>

{body}

{mutation_legend_html()}

</body>
</html>
"""


def msa_to_html(
    records,
    title="Multiple Sequence Alignment",
    line_width=80,
    analysis_mode: str = ANALYSIS_MODE_CONSERVATIVE,
    display_style: str = DISPLAY_STYLE_AA_CLASS,
    color_config=None,
):
    if analysis_mode == ANALYSIS_MODE_MUTATION:
        return msa_to_html_mutation(records, title=title, line_width=line_width)

    if display_style == DISPLAY_STYLE_ESPRIPT:
        return msa_to_html_espript(
            records,
            title=title,
            line_width=line_width,
            color_config=color_config,
        )

    return msa_to_html_conservative(records, title=title, line_width=line_width)


def msa_to_html_conservative(records, title="Multiple Sequence Alignment", line_width=80):
    if not records:
        raise ValueError("No MSA records to export.")

    aln_len = len(records[0]["sequence"])
    max_name_len = max(len(record.get("display_name", record.get("id", ""))) for record in records)
    column_classes = get_column_classes_for_msa(records, ANALYSIS_MODE_CONSERVATIVE)

    blocks = []

    for start in range(0, aln_len, line_width):
        end = min(start + line_width, aln_len)
        block_column_classes = column_classes[start:end]
        block_edge_classes = get_run_edges(block_column_classes)
        block_lines = []

        for record in records:
            name = html.escape(record.get("display_name", record.get("id", "sequence")).ljust(max_name_len))
            seq = record["sequence"][start:end]

            seq_html = []

            for offset, aa in enumerate(seq):
                seq_html.append(
                    aa_span_conservative(
                        aa,
                        block_column_classes[offset],
                        block_edge_classes[offset],
                    )
                )

            block_lines.append(
                f'<div><span class="label">{name}</span> <span class="seq">{"".join(seq_html)}</span></div>'
            )

        blocks.append(
            f"""
<div class="alignment-block">
<div class="position">Position {start + 1}-{end}</div>
{''.join(block_lines)}
</div>
"""
        )

    body = "\n".join(blocks)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{html.escape(title)}</title>
<style>
{base_css()}
{conservative_css()}
</style>
</head>
<body>

<h1>{html.escape(title)}</h1>
<h2>Analysis mode: Conservative analysis | Display style: Amino acid class</h2>

<div class="summary">
<div><strong>Number of sequences:</strong> {len(records)}</div>
<div><strong>Alignment length:</strong> {aln_len}</div>
</div>

{body}

{conservative_legend_html()}

</body>
</html>
"""


def msa_to_html_espript(records, title="Multiple Sequence Alignment", line_width=80, color_config=None):
    if not records:
        raise ValueError("No MSA records to export.")

    colors = merged_colors(color_config)

    aln_len = len(records[0]["sequence"])
    max_name_len = max(len(record.get("display_name", record.get("id", ""))) for record in records)
    column_classes = espript_column_classes(records)

    blocks = []

    for start in range(0, aln_len, line_width):
        end = min(start + line_width, aln_len)
        block_column_classes = column_classes[start:end]
        block_edge_classes = get_run_edges(block_column_classes)
        block_lines = []

        for record in records:
            name = html.escape(record.get("display_name", record.get("id", "sequence")).ljust(max_name_len))
            seq = record["sequence"][start:end]

            seq_html = []

            for offset, aa in enumerate(seq):
                seq_html.append(
                    aa_span_espript(
                        aa,
                        block_column_classes[offset],
                        block_edge_classes[offset],
                    )
                )

            block_lines.append(
                f'<div><span class="label espript-name">{name}</span> <span class="seq">{"".join(seq_html)}</span></div>'
            )

        blocks.append(
            f"""
<div class="alignment-block espript-block">
<div class="position">Position {start + 1}-{end}</div>
{''.join(block_lines)}
</div>
"""
        )

    body = "\n".join(blocks)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{html.escape(title)}</title>
<style>
{base_css()}
{espript_css(colors)}
</style>
</head>
<body>

<h1>{html.escape(title)}</h1>
<h2>Analysis mode: Conservative analysis | Display style: ESPript-like conservation</h2>

<div class="summary">
<div><strong>Number of sequences:</strong> {len(records)}</div>
<div><strong>Alignment length:</strong> {aln_len}</div>
<div><strong>Color rule:</strong> Muted conserved background, soft similar residue color, and low-saturation similarity frame.</div>
</div>

{body}

{espript_legend_html(colors)}

</body>
</html>
"""


def msa_to_html_mutation(records, title="Multiple Sequence Alignment", line_width=80):
    if not records:
        raise ValueError("No MSA records to export.")

    aln_len = len(records[0]["sequence"])
    max_name_len = max(len(record.get("display_name", record.get("id", ""))) for record in records)
    column_classes = get_column_classes_for_msa(records, ANALYSIS_MODE_MUTATION)

    blocks = []

    for start in range(0, aln_len, line_width):
        end = min(start + line_width, aln_len)
        block_column_classes = column_classes[start:end]
        block_edge_classes = get_run_edges(block_column_classes)
        block_lines = []

        for record in records:
            name = html.escape(record.get("display_name", record.get("id", "sequence")).ljust(max_name_len))
            seq = record["sequence"][start:end]

            seq_html = []

            for offset, aa in enumerate(seq):
                seq_html.append(
                    aa_span_mutation(
                        aa,
                        block_column_classes[offset],
                        block_edge_classes[offset],
                    )
                )

            block_lines.append(
                f'<div><span class="label">{name}</span> <span class="seq">{"".join(seq_html)}</span></div>'
            )

        blocks.append(
            f"""
<div class="alignment-block mutation-block">
<div class="position">Position {start + 1}-{end}</div>
{''.join(block_lines)}
</div>
"""
        )

    body = "\n".join(blocks)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{html.escape(title)}</title>
<style>
{base_css()}
{mutation_css()}
</style>
</head>
<body>

<h1>{html.escape(title)}</h1>
<h2>Analysis mode: Mutation site comparison</h2>

<div class="summary">
<div><strong>Number of sequences:</strong> {len(records)}</div>
<div><strong>Alignment length:</strong> {aln_len}</div>
</div>

{body}

{mutation_legend_html()}

</body>
</html>
"""


def conservative_legend_html():
    return """
<div class="legend">
<div class="legend-section">
<strong>Amino acid classes:</strong><br><br>
<span class="aa hydrophobic">A V L I M</span> Hydrophobic
<span class="aa aromatic">F W Y</span> Aromatic
<span class="aa positive">K R H</span> Positive
<span class="aa negative">D E</span> Negative
<span class="aa polar">S T N Q</span> Polar
<span class="aa special">G P C</span> Special
<span class="aa gap">-</span> Gap
</div>

<div class="legend-section">
<strong>Column conservation:</strong><br><br>
<span class="aa hydrophobic identical_column run_start run_end">A</span> Fully conserved column
<span class="aa hydrophobic variable_column run_start run_end">A</span> Variable column
<span class="aa gap gap_column run_start run_end">-</span> Gap-containing column
</div>
</div>
"""


def mutation_legend_html():
    return """
<div class="legend">
<div class="legend-section">
<strong>Mutation site comparison:</strong><br><br>
<span class="aa mutation_plain identical_column">A</span> Same / unchanged
<span class="aa mutation_plain mutation_column">A</span> Different / mutation site
<span class="aa mutation_plain gap_column">-</span> Gap / insertion-deletion
</div>
</div>
"""


def espript_legend_html(colors):
    return f"""
<div class="legend">
<div class="legend-section">
<strong>ESPript-like conservation style:</strong><br><br>
<span class="aa espript espript_conserved run_start run_end">A</span> Fully conserved column
<span class="aa espript espript_strong_similar run_start run_end">A</span> Strongly similar column
<span class="aa espript espript_weak_similar">A</span> Weakly similar column
<span class="aa espript espript_nonconserved">A</span> Non-conserved column
<span class="aa espript espript_gap run_start run_end">.</span> Gap-containing column
</div>
<div class="legend-section">
<strong>Current ESPript colors:</strong>
Conserved bg {colors["conserved_bg"]};
Similar text {colors["similar_text"]};
Similar border {colors["similar_border"]};
Gap bg {colors["gap_bg"]}.
</div>
</div>
"""


def base_css():
    return """
body {
    font-family: Arial, sans-serif;
    margin: 30px;
    background: #ffffff;
    color: #222222;
}

h1 {
    font-size: 26px;
    margin-bottom: 5px;
}

h2 {
    font-size: 16px;
    margin-top: 0;
    color: #555555;
}

.summary {
    margin-bottom: 25px;
    padding: 15px;
    border: 1px solid #dddddd;
    border-radius: 8px;
    background: #fafafa;
    font-size: 14px;
}

.summary-grid {
    display: grid;
    grid-template-columns: 260px 120px;
    row-gap: 6px;
}

.alignment-block {
    font-family: "Courier New", Menlo, Monaco, monospace;
    font-size: 18px;
    line-height: 1.45;
    margin-bottom: 22px;
    white-space: pre;
}

.position {
    font-family: Arial, sans-serif;
    font-size: 12px;
    color: #666666;
    margin-bottom: 5px;
}

.label {
    display: inline-block;
    min-width: 18ch;
    font-weight: bold;
    color: #111111;
}

.seq span {
    display: inline-block;
    min-width: 13px;
    height: 22px;
    line-height: 22px;
    text-align: center;
    box-sizing: border-box;
}

.legend {
    margin-top: 30px;
    padding: 12px;
    border: 1px solid #dddddd;
    border-radius: 8px;
    background: #fafafa;
    font-size: 14px;
}

.legend-section {
    margin-top: 12px;
}

.legend span {
    padding: 3px 8px;
    margin-right: 8px;
    font-family: "Courier New", monospace;
    border-radius: 3px;
}
"""


def conservative_css():
    return """
.aa.hydrophobic { background: #fff3a3; }
.aa.aromatic { background: #ffd39a; }
.aa.positive { background: #b8dcff; }
.aa.negative { background: #ffc0d6; }
.aa.polar { background: #c7ffbf; }
.aa.special { background: #e5e5e5; }
.aa.gap { background: #cfcfcf; color: #555555; }

.identical,
.strong_conservative,
.weak_similar,
.neutral,
.identical_column,
.variable_column,
.gap_column {
    border-top: 1.5px solid transparent;
    border-bottom: 1.5px solid transparent;
}

.identical,
.identical_column {
    font-weight: bold;
    border-top-color: #1b7f3a;
    border-bottom-color: #1b7f3a;
}

.strong_conservative {
    border-top-color: #40916c;
    border-bottom-color: #40916c;
}

.weak_similar,
.neutral {
    border-top-color: #999999;
    border-bottom-color: #999999;
}

.gap_column {
    border-top-color: #999999;
    border-bottom-color: #999999;
}

.run_start {
    border-left: 1.5px solid currentColor;
}

.run_end {
    border-right: 1.5px solid currentColor;
}

.identical.run_start,
.identical_column.run_start {
    border-left-color: #1b7f3a;
}

.identical.run_end,
.identical_column.run_end {
    border-right-color: #1b7f3a;
}

.strong_conservative.run_start {
    border-left-color: #40916c;
}

.strong_conservative.run_end {
    border-right-color: #40916c;
}

.weak_similar.run_start,
.neutral.run_start,
.gap_column.run_start {
    border-left-color: #999999;
}

.weak_similar.run_end,
.neutral.run_end,
.gap_column.run_end {
    border-right-color: #999999;
}

.match { font-weight: bold; }
.match.identical { color: #006400; }
.match.strong_conservative { color: #2d6a4f; }
.match.weak_similar,
.match.neutral { color: #555555; }
.match.non_conservative,
.match.gap { color: #aaaaaa; }
"""


def mutation_css():
    return """
.aa.mutation_plain {
    background: #ffffff !important;
    color: #111111 !important;
    font-weight: normal !important;
    border: none !important;
    outline: none !important;
}

.aa.mutation_plain.identical,
.identical_column {
    background: #ffffff !important;
    color: #111111 !important;
    font-weight: normal !important;
    border: none !important;
    outline: none !important;
}

.aa.mutation_plain.mutation,
.mutation_column {
    background: #C85A5A !important;
    color: #ffffff !important;
    font-weight: bold !important;
    border: none !important;
    outline: none !important;
}

.aa.mutation_plain.gap,
.gap_column {
    background: #ffffff !important;
    color: #111111 !important;
    font-weight: normal !important;
    border: none !important;
    outline: none !important;
}

.run_start,
.run_end {
    border: none !important;
    outline: none !important;
}

.match {
    font-weight: normal;
}

.match.identical {
    color: #BBBBBB;
}

.match.mutation {
    color: #C85A5A;
    font-weight: bold;
}

.match.gap {
    color: #111111;
}
"""


def espript_css(colors):
    return f"""
.espript-block {{
    font-size: 18px;
}}

.espript-name {{
    font-style: italic;
    font-weight: bold;
}}

.aa.espript {{
    background: #ffffff;
    color: #111111;
    font-weight: normal;
    border-top: 1.5px solid transparent;
    border-bottom: 1.5px solid transparent;
}}

.aa.espript_conserved {{
    background: {colors["conserved_bg"]};
    color: {colors["conserved_text"]};
    font-weight: bold;
}}

.aa.espript_strong_similar {{
    background: #ffffff;
    color: {colors["similar_text"]};
    font-weight: bold;
    border-top-color: {colors["similar_border"]};
    border-bottom-color: {colors["similar_border"]};
}}

.aa.espript_weak_similar {{
    background: #ffffff;
    color: {colors["similar_text"]};
    font-weight: bold;
}}

.aa.espript_nonconserved {{
    background: #ffffff;
    color: #111111;
}}

.aa.espript_gap {{
    background: {colors["gap_bg"]};
    color: {colors["gap_text"]};
}}

.run_start.espript_strong_similar,
.run_start.espript_conserved,
.run_start.espript_gap {{
    border-left: 1.5px solid currentColor;
}}

.run_end.espript_strong_similar,
.run_end.espript_conserved,
.run_end.espript_gap {{
    border-right: 1.5px solid currentColor;
}}

.run_start.espript_strong_similar {{
    border-left-color: {colors["similar_border"]};
}}

.run_end.espript_strong_similar {{
    border-right-color: {colors["similar_border"]};
}}

.run_start.espript_conserved {{
    border-left-color: {colors["conserved_bg"]};
}}

.run_end.espript_conserved {{
    border-right-color: {colors["conserved_bg"]};
}}

.run_start.espript_gap {{
    border-left-color: #BBBBBB;
}}

.run_end.espript_gap {{
    border-right-color: #BBBBBB;
}}
"""
