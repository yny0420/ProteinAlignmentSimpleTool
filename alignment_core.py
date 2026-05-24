from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices


AA_GROUPS = {
    "hydrophobic": set("AVLIM"),
    "aromatic": set("FWY"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "polar": set("STNQ"),
    "special": set("GPC"),
    "amide": set("NQ"),
    "small": set("AGS"),
    "tiny": set("AGCST"),
}


AA_GROUP_DISPLAY = {
    "hydrophobic": "Hydrophobic",
    "aromatic": "Aromatic",
    "positive": "Positive",
    "negative": "Negative",
    "polar": "Polar",
    "special": "Special",
    "amide": "Amide",
    "small": "Small",
    "tiny": "Tiny",
}


ANALYSIS_MODE_CONSERVATIVE = "conservative"
ANALYSIS_MODE_MUTATION = "mutation"


def clean_sequence(seq: str) -> str:
    """
    清理蛋白序列。
    """
    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO-")
    seq = seq.upper()
    return "".join([c for c in seq if c in allowed])


def get_aa_primary_class(aa: str) -> str:
    """
    返回氨基酸主要类型，用于 conservative analysis 着色。
    """
    aa = aa.upper()

    if aa == "-":
        return "gap"

    if aa in set("AVLIM"):
        return "hydrophobic"

    if aa in set("FWY"):
        return "aromatic"

    if aa in set("KRH"):
        return "positive"

    if aa in set("DE"):
        return "negative"

    if aa in set("STNQ"):
        return "polar"

    if aa in set("GPC"):
        return "special"

    return "unknown"


def shared_groups(a: str, b: str):
    """
    返回两个氨基酸共享的理化性质分组。
    """
    if a == "-" or b == "-":
        return []

    groups = []

    for group_name, residues in AA_GROUPS.items():
        if a in residues and b in residues:
            groups.append(group_name)

    return groups


def is_physicochemical_similar(a: str, b: str) -> bool:
    if a == b:
        return True

    if a == "-" or b == "-":
        return False

    return len(shared_groups(a, b)) > 0


def get_blosum62_score(a: str, b: str):
    """
    获取两个氨基酸在 BLOSUM62 中的替换分数。
    gap 返回 None。
    """
    if a == "-" or b == "-":
        return None

    try:
        matrix = substitution_matrices.load("BLOSUM62")
        return matrix[a, b]
    except Exception:
        return None


def classify_substitution_conservative(a: str, b: str):
    """
    保守型分析模式：
    |  identical
    :  strong conservative, BLOSUM62 > 0
    .  weak / neutral similarity
       non-conservative or gap
    """
    a = a.upper()
    b = b.upper()

    if a == "-" or b == "-":
        return {
            "symbol": " ",
            "type": "gap",
            "description": "Gap",
            "blosum62": None,
            "shared_groups": [],
        }

    if a == b:
        return {
            "symbol": "|",
            "type": "identical",
            "description": "Identical residue",
            "blosum62": get_blosum62_score(a, b),
            "shared_groups": shared_groups(a, b),
        }

    blosum_score = get_blosum62_score(a, b)
    groups = shared_groups(a, b)

    if blosum_score is not None and blosum_score > 0:
        return {
            "symbol": ":",
            "type": "strong_conservative",
            "description": "Conservative substitution with positive BLOSUM62 score",
            "blosum62": blosum_score,
            "shared_groups": groups,
        }

    if len(groups) > 0:
        return {
            "symbol": ".",
            "type": "weak_similar",
            "description": "Physicochemically similar substitution",
            "blosum62": blosum_score,
            "shared_groups": groups,
        }

    if blosum_score == 0:
        return {
            "symbol": ".",
            "type": "neutral",
            "description": "Neutral substitution in BLOSUM62",
            "blosum62": blosum_score,
            "shared_groups": groups,
        }

    return {
        "symbol": " ",
        "type": "non_conservative",
        "description": "Non-conservative substitution",
        "blosum62": blosum_score,
        "shared_groups": groups,
    }


def classify_substitution_mutation(a: str, b: str):
    """
    突变位点比对模式：
    只区分一样、不一样、gap。
    不做保守性判断，不看 BLOSUM62，不看理化性质。
    """
    a = a.upper()
    b = b.upper()

    if a == "-" or b == "-":
        return {
            "symbol": "-",
            "type": "gap",
            "description": "Gap",
            "blosum62": None,
            "shared_groups": [],
        }

    if a == b:
        return {
            "symbol": "|",
            "type": "identical",
            "description": "Identical residue",
            "blosum62": None,
            "shared_groups": [],
        }

    return {
        "symbol": "*",
        "type": "mutation",
        "description": "Different residue",
        "blosum62": None,
        "shared_groups": [],
    }


def classify_substitution(a: str, b: str, analysis_mode: str = ANALYSIS_MODE_CONSERVATIVE):
    """
    根据分析模式分类 alignment column。
    """
    if analysis_mode == ANALYSIS_MODE_MUTATION:
        return classify_substitution_mutation(a, b)

    return classify_substitution_conservative(a, b)


def is_conservative_substitution(a: str, b: str) -> bool:
    info = classify_substitution_conservative(a, b)
    return info["type"] in ["identical", "strong_conservative", "weak_similar", "neutral"]


def align_pairwise(seq_a: str, seq_b: str, analysis_mode: str = ANALYSIS_MODE_CONSERVATIVE):
    """
    两条蛋白序列 global pairwise alignment。
    alignment 本身仍使用 BLOSUM62。
    但输出解释模式由 analysis_mode 决定。
    """
    seq_a = clean_sequence(seq_a)
    seq_b = clean_sequence(seq_b)

    if not seq_a or not seq_b:
        raise ValueError("Both Sequence A and Sequence B are required.")

    aligner = PairwiseAligner()
    aligner.mode = "global"

    matrix = substitution_matrices.load("BLOSUM62")
    aligner.substitution_matrix = matrix

    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq_a, seq_b)

    if len(alignments) == 0:
        raise RuntimeError("No alignment result was generated.")

    best_alignment = alignments[0]

    aligned_a, aligned_b = format_alignment_strings(best_alignment, seq_a, seq_b)
    match_line = build_match_line(aligned_a, aligned_b, analysis_mode=analysis_mode)
    column_annotations = build_column_annotations(aligned_a, aligned_b, analysis_mode=analysis_mode)

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        stats = calculate_mutation_statistics(aligned_a, aligned_b)
    else:
        stats = calculate_conservative_statistics(aligned_a, aligned_b)

    return {
        "analysis_mode": analysis_mode,
        "aligned_seq_a": aligned_a,
        "match_line": match_line,
        "aligned_seq_b": aligned_b,
        "column_annotations": column_annotations,
        **stats,
        "score": best_alignment.score,
    }


def format_alignment_strings(alignment, seq_a: str, seq_b: str):
    coordinates = alignment.coordinates

    aligned_a = []
    aligned_b = []

    for i in range(coordinates.shape[1] - 1):
        start_a, end_a = coordinates[0, i], coordinates[0, i + 1]
        start_b, end_b = coordinates[1, i], coordinates[1, i + 1]

        len_a = end_a - start_a
        len_b = end_b - start_b

        if len_a > 0 and len_b > 0:
            aligned_a.append(seq_a[start_a:end_a])
            aligned_b.append(seq_b[start_b:end_b])

        elif len_a > 0 and len_b == 0:
            aligned_a.append(seq_a[start_a:end_a])
            aligned_b.append("-" * len_a)

        elif len_a == 0 and len_b > 0:
            aligned_a.append("-" * len_b)
            aligned_b.append(seq_b[start_b:end_b])

    return "".join(aligned_a), "".join(aligned_b)


def build_column_annotations(
    aligned_a: str,
    aligned_b: str,
    analysis_mode: str = ANALYSIS_MODE_CONSERVATIVE,
):
    annotations = []

    for index, (a, b) in enumerate(zip(aligned_a, aligned_b), start=1):
        info = classify_substitution(a, b, analysis_mode=analysis_mode)

        annotations.append({
            "position": index,
            "aa_a": a,
            "aa_b": b,
            **info,
        })

    return annotations


def build_match_line(
    aligned_a: str,
    aligned_b: str,
    analysis_mode: str = ANALYSIS_MODE_CONSERVATIVE,
) -> str:
    line = []

    for a, b in zip(aligned_a, aligned_b):
        info = classify_substitution(a, b, analysis_mode=analysis_mode)
        line.append(info["symbol"])

    return "".join(line)


def calculate_conservative_statistics(aligned_a: str, aligned_b: str):
    total = len(aligned_a)

    counts = {
        "identical": 0,
        "strong_conservative": 0,
        "weak_similar": 0,
        "neutral": 0,
        "non_conservative": 0,
        "gap": 0,
    }

    if total == 0:
        return {
            "identity": 0.0,
            "strong_conservative": 0.0,
            "weak_similar": 0.0,
            "similarity": 0.0,
            "non_conservative": 0.0,
            "gap_percent": 0.0,
            "mutation_percent": 0.0,
            "counts": counts,
        }

    for a, b in zip(aligned_a, aligned_b):
        info = classify_substitution_conservative(a, b)
        sub_type = info["type"]

        if sub_type in counts:
            counts[sub_type] += 1

    identical = counts["identical"]
    strong = counts["strong_conservative"]
    weak = counts["weak_similar"] + counts["neutral"]
    gaps = counts["gap"]
    non_conservative = counts["non_conservative"]

    identity = identical / total * 100
    strong_conservative = strong / total * 100
    weak_similar = weak / total * 100
    similarity = (identical + strong + weak) / total * 100
    non_conservative_percent = non_conservative / total * 100
    gap_percent = gaps / total * 100

    return {
        "identity": identity,
        "strong_conservative": strong_conservative,
        "weak_similar": weak_similar,
        "similarity": similarity,
        "non_conservative": non_conservative_percent,
        "gap_percent": gap_percent,
        "mutation_percent": non_conservative_percent,
        "counts": counts,
    }


def calculate_mutation_statistics(aligned_a: str, aligned_b: str):
    """
    突变位点比对模式统计：
    只统计 identical、mutation、gap。
    """
    total = len(aligned_a)

    counts = {
        "identical": 0,
        "mutation": 0,
        "gap": 0,
    }

    if total == 0:
        return {
            "identity": 0.0,
            "mutation_percent": 0.0,
            "gap_percent": 0.0,
            "strong_conservative": 0.0,
            "weak_similar": 0.0,
            "similarity": 0.0,
            "non_conservative": 0.0,
            "counts": counts,
        }

    for a, b in zip(aligned_a, aligned_b):
        info = classify_substitution_mutation(a, b)
        counts[info["type"]] += 1

    identity = counts["identical"] / total * 100
    mutation_percent = counts["mutation"] / total * 100
    gap_percent = counts["gap"] / total * 100

    return {
        "identity": identity,
        "mutation_percent": mutation_percent,
        "gap_percent": gap_percent,
        "strong_conservative": 0.0,
        "weak_similar": 0.0,
        "similarity": identity,
        "non_conservative": mutation_percent,
        "counts": counts,
    }


def calculate_pairwise_identity_from_aligned(seq_a: str, seq_b: str) -> float:
    if len(seq_a) != len(seq_b):
        raise ValueError("Aligned sequences must have the same length.")

    comparable = 0
    identical = 0

    for a, b in zip(seq_a, seq_b):
        if a == "-" and b == "-":
            continue

        comparable += 1

        if a == b and a != "-":
            identical += 1

    if comparable == 0:
        return 0.0

    return identical / comparable * 100


def calculate_pairwise_identity_matrix(records):
    matrix = []

    for i, rec_i in enumerate(records):
        row = []

        for j, rec_j in enumerate(records):
            if i == j:
                row.append(100.0)
            else:
                result = align_pairwise(
                    rec_i["sequence"],
                    rec_j["sequence"],
                    analysis_mode=ANALYSIS_MODE_MUTATION,
                )
                row.append(result["identity"])

        matrix.append(row)

    return matrix


def format_alignment_output(result: dict, line_width: int = 80) -> str:
    """
    根据 analysis_mode 输出不同风格的纯文本结果。
    """
    analysis_mode = result.get("analysis_mode", ANALYSIS_MODE_CONSERVATIVE)

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        return format_mutation_alignment_output(result, line_width=line_width)

    return format_conservative_alignment_output(result, line_width=line_width)


def format_conservative_alignment_output(result: dict, line_width: int = 80) -> str:
    aligned_a = result["aligned_seq_a"]
    match_line = result["match_line"]
    aligned_b = result["aligned_seq_b"]

    counts = result.get("counts", {})

    output_lines = []

    output_lines.append("Pairwise Protein Sequence Alignment")
    output_lines.append("Analysis mode: Conservative analysis")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(f"Identity:                  {result['identity']:.2f}%")
    output_lines.append(f"Strong conservative:        {result['strong_conservative']:.2f}%")
    output_lines.append(f"Weak / neutral similarity:  {result['weak_similar']:.2f}%")
    output_lines.append(f"Total similarity:           {result['similarity']:.2f}%")
    output_lines.append(f"Non-conservative:           {result['non_conservative']:.2f}%")
    output_lines.append(f"Gaps:                       {result['gap_percent']:.2f}%")
    output_lines.append(f"Alignment score:            {result['score']:.2f}")
    output_lines.append("")
    output_lines.append("Counts:")
    output_lines.append(f"  Identical:                {counts.get('identical', 0)}")
    output_lines.append(f"  Strong conservative:      {counts.get('strong_conservative', 0)}")
    output_lines.append(f"  Weak similar:             {counts.get('weak_similar', 0)}")
    output_lines.append(f"  Neutral:                  {counts.get('neutral', 0)}")
    output_lines.append(f"  Non-conservative:         {counts.get('non_conservative', 0)}")
    output_lines.append(f"  Gap:                      {counts.get('gap', 0)}")
    output_lines.append("")
    output_lines.append("Symbols:")
    output_lines.append("  |  Identical residue")
    output_lines.append("  :  Strong conservative substitution, usually positive BLOSUM62 score")
    output_lines.append("  .  Weak / neutral similarity")
    output_lines.append("     Non-conservative substitution or gap")
    output_lines.append("")
    output_lines.append("=" * 70)
    output_lines.append("")

    for i in range(0, len(aligned_a), line_width):
        block_a = aligned_a[i:i + line_width]
        block_m = match_line[i:i + line_width]
        block_b = aligned_b[i:i + line_width]

        start_pos = i + 1
        end_pos = i + len(block_a)

        output_lines.append(f"Position {start_pos}-{end_pos}")
        output_lines.append(f"SeqA  {block_a}")
        output_lines.append(f"      {block_m}")
        output_lines.append(f"SeqB  {block_b}")
        output_lines.append("")

    return "\n".join(output_lines)


def format_mutation_alignment_output(result: dict, line_width: int = 80) -> str:
    aligned_a = result["aligned_seq_a"]
    match_line = result["match_line"]
    aligned_b = result["aligned_seq_b"]

    counts = result.get("counts", {})

    output_lines = []

    output_lines.append("Pairwise Protein Sequence Alignment")
    output_lines.append("Analysis mode: Mutation site comparison")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(f"Identity:        {result['identity']:.2f}%")
    output_lines.append(f"Mutation sites:  {result['mutation_percent']:.2f}%")
    output_lines.append(f"Gaps:            {result['gap_percent']:.2f}%")
    output_lines.append(f"Alignment score: {result['score']:.2f}")
    output_lines.append("")
    output_lines.append("Counts:")
    output_lines.append(f"  Identical:     {counts.get('identical', 0)}")
    output_lines.append(f"  Mutation:      {counts.get('mutation', 0)}")
    output_lines.append(f"  Gap:           {counts.get('gap', 0)}")
    output_lines.append("")
    output_lines.append("Symbols:")
    output_lines.append("  |  Identical")
    output_lines.append("  *  Different")
    output_lines.append("  -  Gap")
    output_lines.append("")
    output_lines.append("=" * 70)
    output_lines.append("")

    for i in range(0, len(aligned_a), line_width):
        block_a = aligned_a[i:i + line_width]
        block_m = match_line[i:i + line_width]
        block_b = aligned_b[i:i + line_width]

        start_pos = i + 1
        end_pos = i + len(block_a)

        output_lines.append(f"Position {start_pos}-{end_pos}")
        output_lines.append(f"SeqA  {block_a}")
        output_lines.append(f"      {block_m}")
        output_lines.append(f"SeqB  {block_b}")
        output_lines.append("")

    return "\n".join(output_lines)


def format_column_detail_table(result: dict, max_rows: int = 500) -> str:
    annotations = result.get("column_annotations", [])

    lines = []
    lines.append("Column-level substitution annotation")
    lines.append("=" * 90)
    lines.append("Pos\tSeqA\tSeqB\tSymbol\tType\tBLOSUM62\tShared groups")
    lines.append("-" * 90)

    for ann in annotations[:max_rows]:
        groups = ann.get("shared_groups", [])
        group_text = ",".join(groups) if groups else "-"
        blosum = ann.get("blosum62")
        blosum_text = "-" if blosum is None else str(blosum)

        lines.append(
            f"{ann['position']}\t"
            f"{ann['aa_a']}\t"
            f"{ann['aa_b']}\t"
            f"{ann['symbol']}\t"
            f"{ann['type']}\t"
            f"{blosum_text}\t"
            f"{group_text}"
        )

    if len(annotations) > max_rows:
        lines.append("")
        lines.append(f"Only the first {max_rows} columns are shown out of {len(annotations)} total columns.")

    return "\n".join(lines)


def parse_domain_text(domain_text: str):
    domains = []

    lines = domain_text.strip().splitlines()

    for line in lines:
        line = line.strip()

        if not line:
            continue

        if ":" not in line:
            raise ValueError(f"Invalid domain format: {line}")

        name, region = line.split(":", 1)
        name = name.strip()
        region = region.strip()

        if "-" not in region:
            raise ValueError(f"Invalid domain range: {line}")

        start_str, end_str = region.split("-", 1)

        try:
            start = int(start_str.strip())
            end = int(end_str.strip())
        except ValueError:
            raise ValueError(f"Invalid domain numbers: {line}")

        if start < 1 or end < start:
            raise ValueError(f"Invalid domain range: {line}")

        domains.append({
            "name": name,
            "start": start,
            "end": end,
        })

    if len(domains) == 0:
        raise ValueError("No valid domain definition was found.")

    return domains


def extract_domains_from_records(records, domains):
    output = {}

    for domain in domains:
        domain_name = domain["name"]
        start = domain["start"]
        end = domain["end"]

        output[domain_name] = []

        for record in records:
            seq = clean_sequence(record["sequence"])

            if end > len(seq):
                raise ValueError(
                    f"Domain {domain_name} {start}-{end} exceeds sequence length "
                    f"of {record['id']} ({len(seq)} aa)."
                )

            subseq = seq[start - 1:end]

            output[domain_name].append({
                "id": f"{record['id']}_{domain_name}_{start}-{end}",
                "display_name": f"{record.get('display_name', record.get('id', 'sequence'))}_{domain_name}",
                "description": f"{record.get('description', record['id'])} | {domain_name}:{start}-{end}",
                "sequence": subseq,
            })

    return output


def records_to_fasta(records) -> str:
    lines = []

    for record in records:
        description = record.get("display_name", record.get("description", record.get("id", "sequence")))
        sequence = record["sequence"]

        lines.append(f">{description}")

        for i in range(0, len(sequence), 80):
            lines.append(sequence[i:i + 80])

    return "\n".join(lines)
