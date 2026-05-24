import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from fasta_utils import write_fasta_file, read_fasta_file
from alignment_core import ANALYSIS_MODE_CONSERVATIVE, ANALYSIS_MODE_MUTATION


def _candidate_runtime_roots():
    roots = []

    try:
        roots.append(Path(__file__).resolve().parent)
    except Exception:
        pass

    if getattr(sys, "frozen", False):
        if hasattr(sys, "_MEIPASS"):
            roots.append(Path(sys._MEIPASS).resolve())

        try:
            exe_dir = Path(sys.executable).resolve().parent
            roots.append(exe_dir)
            roots.append(exe_dir.parent)
            roots.append(exe_dir.parent / "Resources")
            roots.append(exe_dir.parent / "Frameworks")
        except Exception:
            pass

    seen = set()
    unique_roots = []

    for root in roots:
        key = str(root)
        if key not in seen:
            seen.add(key)
            unique_roots.append(root)

    return unique_roots


def _bundled_mafft_dirs():
    dirs = []

    for root in _candidate_runtime_roots():
        dirs.extend([
            root / "mafft_bundle",
            root / "_internal" / "mafft_bundle",
            root / "Resources" / "mafft_bundle",
            root / "Frameworks" / "mafft_bundle",
        ])

    seen = set()
    unique_dirs = []

    for d in dirs:
        key = str(d)
        if key not in seen:
            seen.add(key)
            unique_dirs.append(d)

    return unique_dirs


def _bundled_mafft_candidates():
    """
    Prefer libexec/bin/mafft over bin/mafft for bundled Homebrew MAFFT.

    The top-level bin/mafft wrapper can be non-relocatable on some old macOS systems.
    The libexec/bin/mafft script plus a correctly configured MAFFT_BINARIES path is
    more reliable inside a PyInstaller .app bundle.
    """
    candidates = []

    for bundle_dir in _bundled_mafft_dirs():
        candidates.extend([
            bundle_dir / "libexec" / "bin" / "mafft",
            bundle_dir / "bin" / "mafft",
            bundle_dir / "mafft",
        ])

    return candidates


def get_bundled_mafft_bundle_dir(mafft_path=None):
    if mafft_path:
        p = Path(mafft_path).resolve()

        for parent in p.parents:
            if parent.name == "mafft_bundle":
                return parent

    for d in _bundled_mafft_dirs():
        if d.exists():
            return d

    return None


def get_bundled_mafft_path():
    for candidate in _bundled_mafft_candidates():
        if candidate.exists() and candidate.is_file():
            try:
                candidate.chmod(candidate.stat().st_mode | 0o111)
            except Exception:
                pass
            return str(candidate)

    return None


def check_mafft_available() -> bool:
    return get_mafft_path() is not None


def get_mafft_path():
    bundled = get_bundled_mafft_path()

    if bundled:
        return bundled

    return shutil.which("mafft")


def sanitize_fasta_name(name: str) -> str:
    name = name.strip()

    if not name:
        name = "sequence"

    safe = []

    for c in name:
        if c.isalnum() or c in ["_", "-", ".", "|"]:
            safe.append(c)
        else:
            safe.append("_")

    return "".join(safe)


def _find_bundled_mafft_binaries_dir(bundle_dir: Path):
    """
    Find the directory containing MAFFT helper binaries.

    MAFFT installations differ slightly by package source. We try several common
    layouts and then scan for helper executables.
    """
    common_dirs = [
        bundle_dir / "libexec" / "core",
        bundle_dir / "libexec" / "mafft",
        bundle_dir / "libexec" / "bin",
        bundle_dir / "lib" / "mafft",
        bundle_dir / "share" / "mafft",
        bundle_dir / "bin",
    ]

    helper_names = [
        "dvtditr",
        "disttbfast",
        "tbfast",
        "mafft-profile",
        "mafft-distance",
        "pairlocalalign",
    ]

    for d in common_dirs:
        if not d.exists() or not d.is_dir():
            continue

        for helper in helper_names:
            if (d / helper).exists():
                return d

    # Deep scan fallback. Limit to mafft_bundle to avoid expensive global search.
    try:
        for d in bundle_dir.rglob("*"):
            if not d.is_dir():
                continue

            for helper in helper_names:
                if (d / helper).exists():
                    return d
    except Exception:
        pass

    return None


def _mafft_environment(mafft_path):
    env = os.environ.copy()

    # Remove user/system MAFFT settings first.
    # Old machines sometimes define MAFFT_BINARIES in shell config files.
    # That breaks bundled MAFFT by forcing it to use incompatible external helpers.
    env.pop("MAFFT_BINARIES", None)
    env.pop("MAFFT_HOME", None)
    env.pop("MAFFT_BINARIES_DIR", None)

    mafft_file = Path(mafft_path).resolve()
    extra_paths = [str(mafft_file.parent)]

    bundle_dir = get_bundled_mafft_bundle_dir(mafft_path)

    if bundle_dir is not None and bundle_dir.exists():
        # Add bundled paths so shell scripts can find helper programs.
        for path in [
            bundle_dir / "bin",
            bundle_dir / "libexec" / "bin",
            bundle_dir / "libexec" / "core",
            bundle_dir / "libexec" / "mafft",
            bundle_dir / "lib" / "mafft",
        ]:
            if path.exists():
                extra_paths.append(str(path))

        binaries_dir = _find_bundled_mafft_binaries_dir(bundle_dir)

        if binaries_dir is not None:
            # This is intentionally set for the bundled app. It prevents the
            # MAFFT wrapper from using stale hard-coded Homebrew paths or user
            # shell values.
            env["MAFFT_BINARIES"] = str(binaries_dir)

    existing_path = env.get("PATH", "")
    env["PATH"] = os.pathsep.join([p for p in extra_paths if p]) + os.pathsep + existing_path
    env.setdefault("TMPDIR", tempfile.gettempdir())

    return env


def run_mafft(records, mafft_mode="auto"):
    if len(records) < 2:
        raise ValueError("At least two sequences are required for MAFFT alignment.")

    mafft_path = get_mafft_path()

    if not mafft_path:
        raise RuntimeError(
            "MAFFT was not found.\n\n"
            "If you are using the packaged app, please use the bundled-MAFFT build.\n"
            "If you are running from source, install MAFFT first:\n"
            "  brew install mafft"
        )

    temp_records = []
    temp_name_to_record = {}

    for idx, record in enumerate(records, start=1):
        user_name = record.get("display_name", record.get("id", f"seq{idx}"))
        temp_id = f"seq_{idx:04d}_{sanitize_fasta_name(user_name)}"

        temp_record = {
            "id": temp_id,
            "display_name": temp_id,
            "description": temp_id,
            "sequence": record["sequence"],
            "_original_record": record,
        }

        temp_records.append(temp_record)
        temp_name_to_record[temp_id] = record

    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = os.path.join(tmpdir, "input.fasta")
        output_fasta = os.path.join(tmpdir, "aligned.fasta")

        write_fasta_file(temp_records, input_fasta)

        cmd = [mafft_path]

        if mafft_mode == "auto":
            cmd.append("--auto")
        elif mafft_mode == "localpair":
            cmd.extend(["--localpair", "--maxiterate", "1000"])
        elif mafft_mode == "globalpair":
            cmd.extend(["--globalpair", "--maxiterate", "1000"])
        elif mafft_mode == "genafpair":
            cmd.extend(["--genafpair", "--maxiterate", "1000"])
        else:
            cmd.append("--auto")

        cmd.append("--inputorder")
        cmd.append(input_fasta)

        env = _mafft_environment(mafft_path)

        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                env=env,
            )
        except Exception as e:
            raise RuntimeError(f"Failed to run MAFFT: {e}")

        if result.returncode != 0:
            debug_lines = [
                "MAFFT failed.",
                "",
                f"MAFFT path:",
                f"{mafft_path}",
                "",
                "Command:",
                " ".join(cmd),
                "",
                "Bundled MAFFT_BINARIES used by app:",
                env.get("MAFFT_BINARIES", "<unset>"),
                "",
                "Error message:",
                result.stderr,
            ]
            raise RuntimeError("\n".join(debug_lines))

        with open(output_fasta, "w", encoding="utf-8") as f:
            f.write(result.stdout)

        aligned_temp_records = read_fasta_file(output_fasta)

    aligned_records = []

    for aligned_temp in aligned_temp_records:
        temp_id = aligned_temp["id"]
        original = temp_name_to_record.get(temp_id)

        if original is None:
            aligned_records.append({
                "id": temp_id,
                "display_name": temp_id,
                "description": temp_id,
                "sequence": aligned_temp["sequence"],
            })
        else:
            aligned_records.append({
                "id": original.get("id", temp_id),
                "display_name": original.get("display_name", original.get("id", temp_id)),
                "description": original.get("description", original.get("display_name", temp_id)),
                "sequence": aligned_temp["sequence"],
            })

    order_map = {
        record.get("display_name", record.get("id", str(i))): i
        for i, record in enumerate(records)
    }

    aligned_records.sort(
        key=lambda r: order_map.get(r.get("display_name", r.get("id", "")), 999999)
    )

    return aligned_records, result.stderr


def build_msa_mutation_line(records, start, end):
    symbols = []

    for col_index in range(start, end):
        column = [record["sequence"][col_index] for record in records]

        if "-" in column:
            symbols.append("-")
        elif len(set(column)) == 1:
            symbols.append("|")
        else:
            symbols.append("*")

    return "".join(symbols)


def format_msa_text(records, line_width=80, analysis_mode: str = ANALYSIS_MODE_CONSERVATIVE):
    if not records:
        return ""

    max_name_len = max(len(record.get("display_name", record.get("id", ""))) for record in records)
    aln_len = len(records[0]["sequence"])

    lines = []
    lines.append("Multiple Sequence Alignment")

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        lines.append("Analysis mode: Mutation site comparison")
    else:
        lines.append("Analysis mode: Conservative analysis")

    lines.append("=" * 60)
    lines.append(f"Number of sequences: {len(records)}")
    lines.append(f"Alignment length:    {aln_len}")
    lines.append("=" * 60)
    lines.append("")

    if analysis_mode == ANALYSIS_MODE_MUTATION:
        lines.append("Symbols:")
        lines.append("  |  Same column")
        lines.append("  *  Different column")
        lines.append("  -  Gap-containing column")
        lines.append("")
        lines.append("=" * 60)
        lines.append("")

    for start in range(0, aln_len, line_width):
        end = min(start + line_width, aln_len)

        lines.append(f"Position {start + 1}-{end}")

        for record in records:
            name = record.get("display_name", record.get("id", "sequence")).ljust(max_name_len)
            block = record["sequence"][start:end]
            lines.append(f"{name}  {block}")

        if analysis_mode == ANALYSIS_MODE_MUTATION:
            marker = build_msa_mutation_line(records, start, end)
            lines.append(f"{'Mutation'.ljust(max_name_len)}  {marker}")

        lines.append("")

    return "\n".join(lines)
