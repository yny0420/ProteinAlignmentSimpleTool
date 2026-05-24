from Bio import SeqIO


def read_fasta_file(file_path: str):
    """
    读取 FASTA 文件。
    每条序列增加 display_name 字段，方便用户在界面中改名。
    """
    records = []

    try:
        for record in SeqIO.parse(file_path, "fasta"):
            records.append({
                "id": record.id,
                "display_name": record.id,
                "description": record.description,
                "sequence": str(record.seq),
            })
    except Exception as e:
        raise RuntimeError(f"Failed to read FASTA file: {e}")

    if len(records) == 0:
        raise ValueError("No FASTA sequences were found in this file.")

    return records


def format_record_for_display(record: dict) -> str:
    """
    把 FASTA record 格式化为显示文本。
    """
    name = record.get("display_name", record.get("id", "sequence"))
    sequence = record["sequence"]
    return f">{name}\n{sequence}"


def write_fasta_file(records, file_path: str):
    """
    写出 FASTA 文件。
    优先使用 display_name。
    """
    with open(file_path, "w", encoding="utf-8") as f:
        for record in records:
            name = record.get("display_name", record.get("id", "sequence"))
            sequence = record["sequence"]

            f.write(f">{name}\n")

            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")
