from pathlib import Path

from Bio.PDB import PDBParser, MMCIFParser


STANDARD_THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "MSE": "M",
    "SEP": "S",
    "TPO": "T",
    "PTR": "Y",
    "CSO": "C",
    "CME": "C",
    "HYP": "P",
    "MLY": "K",
    "MLZ": "K",
    "KCX": "K",
    "LLP": "K",
    "PCA": "E",
    "GLX": "Z",
    "ASX": "B",
    "UNK": "X",
}


def residue_to_one_letter(residue):
    resname = residue.get_resname().strip().upper()

    if resname in STANDARD_THREE_TO_ONE:
        return STANDARD_THREE_TO_ONE[resname]

    # Try Biopython's extended dictionary when available.
    try:
        from Bio.Data.PDBData import protein_letters_3to1_extended

        if resname in protein_letters_3to1_extended:
            return protein_letters_3to1_extended[resname]
    except Exception:
        pass

    return None


def read_structure_file_sequences(file_path: str):
    """
    Extract chain sequences from a PDB or mmCIF file.

    Returns records compatible with the app:
    [
        {
            "id": "...",
            "display_name": "...",
            "description": "...",
            "sequence": "...",
            "source_file": "...",
            "source_chain": "A",
        },
        ...
    ]
    """
    path = Path(file_path)
    suffix = path.suffix.lower()

    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif suffix in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure file type: {path.suffix}")

    structure = parser.get_structure(path.stem, str(path))

    models = list(structure.get_models())

    if not models:
        raise ValueError("No model was found in this structure file.")

    model = models[0]
    records = []

    for chain in model:
        chain_id = chain.id.strip() if chain.id.strip() else "blank"

        residues = []
        residue_numbers = []

        for residue in chain:
            aa = residue_to_one_letter(residue)

            if aa is None:
                continue

            residues.append(aa)

            resseq = residue.id[1]
            icode = residue.id[2].strip()
            residue_numbers.append(f"{resseq}{icode}" if icode else str(resseq))

        if not residues:
            continue

        sequence = "".join(residues)
        first_res = residue_numbers[0] if residue_numbers else "?"
        last_res = residue_numbers[-1] if residue_numbers else "?"

        display_name = f"{path.stem}_chain_{chain_id}"
        description = (
            f"{display_name} | extracted from {path.name} | "
            f"chain {chain_id} | residues {first_res}-{last_res} | length {len(sequence)}"
        )

        records.append({
            "id": display_name,
            "display_name": display_name,
            "description": description,
            "sequence": sequence,
            "source_file": path.name,
            "source_chain": chain_id,
            "source_residue_range": f"{first_res}-{last_res}",
        })

    if not records:
        raise ValueError(
            "No protein chain sequence was extracted. "
            "The file may contain only nucleic acids, ligands, waters, or unsupported residue names."
        )

    return records


def write_records_to_fasta(records, file_path: str):
    """
    Save extracted chain sequence records as FASTA.
    """
    with open(file_path, "w", encoding="utf-8") as f:
        for record in records:
            description = record.get("description", record.get("display_name", record.get("id", "sequence")))
            sequence = record["sequence"]

            f.write(f">{description}\n")

            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")
