import os
import sys
from pathlib import Path

# Clear potentially conflicting Qt environment variables.
# This is usually safer than forcing a fixed Qt plugin path on macOS.
os.environ.pop("QT_PLUGIN_PATH", None)
os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH", None)

from PyQt6.QtWidgets import (
    QApplication,
    QWidget,
    QLabel,
    QTextEdit,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QFileDialog,
    QMessageBox,
    QGroupBox,
    QGridLayout,
    QComboBox,
    QListWidget,
    QListWidgetItem,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QAbstractItemView,
    QInputDialog,
    QColorDialog,
)
from PyQt6.QtGui import QFont, QColor
from PyQt6.QtCore import Qt

from alignment_core import (
    ANALYSIS_MODE_CONSERVATIVE,
    ANALYSIS_MODE_MUTATION,
    align_pairwise,
    format_alignment_output,
    parse_domain_text,
    extract_domains_from_records,
    calculate_pairwise_identity_matrix,
)
from fasta_utils import read_fasta_file, format_record_for_display, write_fasta_file
from structure_utils import read_structure_file_sequences, write_records_to_fasta
from mafft_utils import run_mafft, format_msa_text, check_mafft_available, get_mafft_path
from html_exporter import (
    pairwise_result_to_html,
    msa_to_html,
    DISPLAY_STYLE_AA_CLASS,
    DISPLAY_STYLE_ESPRIPT,
    DEFAULT_ESPRIPT_COLORS,
)
from pdf_exporter import export_pairwise_pdf, export_msa_pdf


class ProteinAlignmentApp(QWidget):
    def __init__(self):
        super().__init__()

        self.fasta_records = []

        self.current_pairwise_result = None
        self.current_pairwise_text = ""
        self.current_pairwise_seq_a_name = "SeqA"
        self.current_pairwise_seq_b_name = "SeqB"
        self.current_pairwise_analysis_mode = ANALYSIS_MODE_CONSERVATIVE

        self.current_msa_records = []
        self.current_msa_text = ""
        self.current_msa_analysis_mode = ANALYSIS_MODE_CONSERVATIVE
        self.current_msa_display_style = DISPLAY_STYLE_AA_CLASS

        self.espript_colors = DEFAULT_ESPRIPT_COLORS.copy()

        self.current_domain_outputs = {}
        self.current_custom_domain_records = []

        self.current_matrix_records = []
        self.current_identity_matrix = []

        self.setWindowTitle("Protein Sequence Alignment Tool")
        self.resize(1350, 940)
        self.setAcceptDrops(True)

        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout()

        title_label = QLabel("Protein Sequence Alignment Tool")
        title_label.setFont(QFont("Arial", 18))
        main_layout.addWidget(title_label)

        self.status_label = QLabel(self.build_status_text())
        self.status_label.setFont(QFont("Arial", 10))
        main_layout.addWidget(self.status_label)

        self.drop_hint_label = QLabel(
            "Drag and drop FASTA / PDB / CIF files anywhere in this window. "
            "Multiple files are supported. Structure files will be converted to chain FASTA records."
        )
        self.drop_hint_label.setFont(QFont("Arial", 10))
        main_layout.addWidget(self.drop_hint_label)

        self.tabs = QTabWidget()

        self.pairwise_tab = QWidget()
        self.msa_tab = QWidget()
        self.domain_tab = QWidget()
        self.matrix_tab = QWidget()

        self.init_pairwise_tab()
        self.init_msa_tab()
        self.init_domain_tab()
        self.init_matrix_tab()

        self.tabs.addTab(self.pairwise_tab, "Pairwise Alignment")
        self.tabs.addTab(self.msa_tab, "Multiple Alignment")
        self.tabs.addTab(self.domain_tab, "Domain Extraction")
        self.tabs.addTab(self.matrix_tab, "Identity Matrix")

        main_layout.addWidget(self.tabs)
        self.setLayout(main_layout)

    def build_status_text(self):
        if check_mafft_available():
            return f"MAFFT detected: {get_mafft_path()}"

        return (
            "MAFFT not detected. Pairwise alignment works, "
            "but multiple sequence alignment requires MAFFT. "
            "Install with: brew install mafft"
        )

    # ============================================================
    # Shared FASTA functions
    # ============================================================

    def load_fasta(self):
        """
        Load one or more FASTA files and append their records to the current sequence list.
        """
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Open FASTA File(s)",
            "",
            "Sequence Files (*.fasta *.fa *.faa *.fas *.txt);;All Files (*)",
        )

        if not file_paths:
            return

        self.process_input_files(file_paths)

    def load_structure_files(self):
        """
        Load one or more PDB/mmCIF files, extract each chain as an independent sequence,
        append to the current sequence list, and optionally export extracted chains as FASTA.
        """
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Open PDB / mmCIF File(s)",
            "",
            "Structure Files (*.pdb *.ent *.cif *.mmcif);;All Files (*)",
        )

        if not file_paths:
            return

        self.process_input_files(file_paths)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                path = Path(url.toLocalFile())
                if self.is_supported_input_file(path):
                    event.acceptProposedAction()
                    return

        event.ignore()

    def dropEvent(self, event):
        paths = []

        for url in event.mimeData().urls():
            local_path = url.toLocalFile()

            if local_path:
                path = Path(local_path)

                if self.is_supported_input_file(path):
                    paths.append(str(path))

        if paths:
            self.process_input_files(paths)
            event.acceptProposedAction()
        else:
            event.ignore()

    def is_supported_input_file(self, path: Path):
        suffix = path.suffix.lower()
        return suffix in {
            ".fasta",
            ".fa",
            ".faa",
            ".fas",
            ".txt",
            ".pdb",
            ".ent",
            ".cif",
            ".mmcif",
        }

    def process_input_files(self, file_paths):
        """
        Import multiple FASTA and/or structure files.
        FASTA records are appended.
        PDB/mmCIF files are parsed chain-by-chain and appended as sequence records.
        """
        fasta_suffixes = {".fasta", ".fa", ".faa", ".fas", ".txt"}
        structure_suffixes = {".pdb", ".ent", ".cif", ".mmcif"}

        added_records = []
        extracted_structure_records = []
        errors = []

        for file_path in file_paths:
            path = Path(file_path)
            suffix = path.suffix.lower()

            try:
                if suffix in fasta_suffixes:
                    records = read_fasta_file(str(path))
                    records = self.add_source_to_records(records, path)
                    added_records.extend(records)

                elif suffix in structure_suffixes:
                    records = read_structure_file_sequences(str(path))
                    extracted_structure_records.extend(records)
                    added_records.extend(records)

                else:
                    errors.append(f"{path.name}: unsupported file type")

            except Exception as e:
                errors.append(f"{path.name}: {e}")

        if added_records:
            self.fasta_records.extend(added_records)
            self.refresh_all_fasta_controls()
            self.fill_pairwise_inputs_after_import()

        message_lines = []

        if added_records:
            message_lines.append(f"Added {len(added_records)} sequence record(s).")

        if extracted_structure_records:
            message_lines.append(
                f"Extracted {len(extracted_structure_records)} chain sequence(s) from structure file(s)."
            )

        if errors:
            message_lines.append("")
            message_lines.append("Some files were not imported:")
            message_lines.extend(errors)

        if message_lines:
            if errors and not added_records:
                QMessageBox.critical(self, "Import Error", "\n".join(message_lines))
            elif errors:
                QMessageBox.warning(self, "Import Completed with Warnings", "\n".join(message_lines))
            else:
                QMessageBox.information(self, "Import Completed", "\n".join(message_lines))

        if extracted_structure_records:
            self.offer_export_extracted_structure_fasta(extracted_structure_records)

    def add_source_to_records(self, records, path: Path):
        """
        Add source file information to FASTA records.
        This keeps user-facing names compact but traceable.
        """
        updated = []

        for record in records:
            rec = dict(record)
            rec.setdefault("source_file", path.name)

            if not rec.get("display_name"):
                rec["display_name"] = rec.get("id", path.stem)

            rec.setdefault("description", rec.get("display_name", rec.get("id", "sequence")))
            updated.append(rec)

        return updated

    def fill_pairwise_inputs_after_import(self):
        """
        Fill pairwise input boxes with the first two loaded records.
        """
        if len(self.fasta_records) >= 1:
            self.seq_a_text.setPlainText(format_record_for_display(self.fasta_records[0]))
            self.current_pairwise_seq_a_name = self.fasta_records[0].get(
                "display_name",
                self.fasta_records[0].get("id", "SeqA"),
            )

        if len(self.fasta_records) >= 2:
            self.seq_b_text.setPlainText(format_record_for_display(self.fasta_records[1]))
            self.record_b_combo.setCurrentIndex(1)
            self.current_pairwise_seq_b_name = self.fasta_records[1].get(
                "display_name",
                self.fasta_records[1].get("id", "SeqB"),
            )

    def offer_export_extracted_structure_fasta(self, extracted_records):
        """
        After structure import, ask whether to immediately save extracted chain sequences.
        """
        reply = QMessageBox.question(
            self,
            "Export extracted chain FASTA?",
            "Structure chain sequences have been extracted.\n\n"
            "Do you want to save these extracted chains as a separate FASTA file now?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )

        if reply != QMessageBox.StandardButton.Yes:
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Extracted Chain FASTA",
            "extracted_chains.fasta",
            "FASTA Files (*.fasta *.fa);;All Files (*)",
        )

        if not file_path:
            return

        try:
            write_records_to_fasta(extracted_records, file_path)
            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def record_label(self, record):
        name = record.get("display_name", record.get("id", "sequence"))
        return f"{name} | length: {len(record['sequence'])}"

    def refresh_all_fasta_controls(self):
        self.refresh_pairwise_combos()
        self.refresh_sequence_list(self.msa_sequence_list)
        self.refresh_sequence_list(self.matrix_sequence_list)

        if hasattr(self, "domain_segment_table"):
            self.refresh_domain_segment_table()

    def refresh_pairwise_combos(self):
        self.record_a_combo.clear()
        self.record_b_combo.clear()

        for record in self.fasta_records:
            display_name = self.record_label(record)
            self.record_a_combo.addItem(display_name)
            self.record_b_combo.addItem(display_name)

    def refresh_sequence_list(self, list_widget):
        list_widget.clear()

        for index, record in enumerate(self.fasta_records):
            item = QListWidgetItem(self.record_label(record))
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsEnabled
            )
            item.setCheckState(Qt.CheckState.Checked)
            item.setData(Qt.ItemDataRole.UserRole, index)
            list_widget.addItem(item)

    def get_checked_records_from_list(self, list_widget):
        records = []

        for i in range(list_widget.count()):
            item = list_widget.item(i)

            if item.checkState() == Qt.CheckState.Checked:
                index = item.data(Qt.ItemDataRole.UserRole)

                if 0 <= index < len(self.fasta_records):
                    records.append(self.fasta_records[index])

        return records

    def get_current_item_index(self, list_widget):
        current_item = list_widget.currentItem()

        if current_item is None:
            QMessageBox.warning(self, "No Selection", "Please select a sequence first.")
            return None

        index = current_item.data(Qt.ItemDataRole.UserRole)

        if not (0 <= index < len(self.fasta_records)):
            QMessageBox.warning(self, "Invalid Selection", "Selected sequence is invalid.")
            return None

        return index

    def rename_record_from_list(self, list_widget):
        index = self.get_current_item_index(list_widget)

        if index is None:
            return

        record = self.fasta_records[index]
        old_name = record.get("display_name", record.get("id", ""))

        new_name, ok = QInputDialog.getText(
            self,
            "Rename Sequence",
            "New sequence name:",
            text=old_name,
        )

        if not ok:
            return

        new_name = new_name.strip()

        if not new_name:
            QMessageBox.warning(self, "Invalid Name", "Sequence name cannot be empty.")
            return

        record["display_name"] = new_name
        self.refresh_all_fasta_controls()

    def move_record_from_list(self, list_widget, direction):
        current_row = list_widget.currentRow()

        if current_row < 0:
            QMessageBox.warning(self, "No Selection", "Please select a sequence first.")
            return

        new_row = current_row + direction

        if new_row < 0 or new_row >= list_widget.count():
            return

        current_item = list_widget.item(current_row)
        target_item = list_widget.item(new_row)

        current_index = current_item.data(Qt.ItemDataRole.UserRole)
        target_index = target_item.data(Qt.ItemDataRole.UserRole)

        if not (0 <= current_index < len(self.fasta_records)):
            return

        if not (0 <= target_index < len(self.fasta_records)):
            return

        self.fasta_records[current_index], self.fasta_records[target_index] = (
            self.fasta_records[target_index],
            self.fasta_records[current_index],
        )

        self.refresh_all_fasta_controls()
        list_widget.setCurrentRow(new_row)

    # ============================================================
    # ESPript-like color customization
    # ============================================================

    def set_espript_color(self, key, title):
        current_hex = self.espript_colors.get(key, DEFAULT_ESPRIPT_COLORS[key])
        color = QColorDialog.getColor(QColor(current_hex), self, title)

        if not color.isValid():
            return

        self.espript_colors[key] = color.name().upper()
        self.update_espript_color_button_labels()

    def reset_espript_colors(self):
        self.espript_colors = DEFAULT_ESPRIPT_COLORS.copy()
        self.update_espript_color_button_labels()

    def update_espript_color_button_labels(self):
        if not hasattr(self, "espript_conserved_bg_button"):
            return

        self.espript_conserved_bg_button.setText(
            f"Conserved bg: {self.espript_colors['conserved_bg']}"
        )
        self.espript_conserved_text_button.setText(
            f"Conserved text: {self.espript_colors['conserved_text']}"
        )
        self.espript_similar_text_button.setText(
            f"Similar text: {self.espript_colors['similar_text']}"
        )
        self.espript_similar_border_button.setText(
            f"Similar border: {self.espript_colors['similar_border']}"
        )
        self.espript_gap_bg_button.setText(
            f"Gap bg: {self.espript_colors['gap_bg']}"
        )


    def get_checked_indices_from_list(self, list_widget):
        """
        Return global record indices for checked boxes in the current list.
        The grey highlighted row is only for Rename / Move Up / Move Down.
        Checked boxes are used for Compare / Delete.
        """
        indices = []

        for i in range(list_widget.count()):
            item = list_widget.item(i)

            if item.checkState() == Qt.CheckState.Checked:
                index = item.data(Qt.ItemDataRole.UserRole)

                if 0 <= index < len(self.fasta_records):
                    indices.append(index)

        return indices

    def delete_checked_records_from_list(self, list_widget):
        """
        Delete all checked sequences from the global sequence list.
        This matches the comparison logic: checked boxes define the active sequence set.
        """
        indices = self.get_checked_indices_from_list(list_widget)

        if not indices:
            QMessageBox.warning(
                self,
                "No Checked Sequences",
                "Please tick the checkbox in front of each sequence you want to delete.",
            )
            return

        names = [
            self.fasta_records[index].get("display_name", self.fasta_records[index].get("id", "sequence"))
            for index in indices
        ]

        preview_names = "\n".join(names[:12])

        if len(names) > 12:
            preview_names += f"\n... and {len(names) - 12} more"

        reply = QMessageBox.question(
            self,
            "Delete checked sequences?",
            f"Delete {len(indices)} checked sequence(s) from the current list?\n\n{preview_names}",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )

        if reply != QMessageBox.StandardButton.Yes:
            return

        for index in sorted(indices, reverse=True):
            del self.fasta_records[index]

        self.refresh_all_fasta_controls()

        if self.fasta_records:
            self.fill_pairwise_inputs_after_import()
        else:
            self.seq_a_text.clear()
            self.seq_b_text.clear()
            self.result_text.clear()
            self.current_pairwise_result = None
            self.current_pairwise_text = ""


    def analysis_mode_slug(self, analysis_mode):
        if analysis_mode == ANALYSIS_MODE_MUTATION:
            return "mutation_site_comparison"
        return "conservative_analysis"

    def display_style_slug(self, display_style):
        if display_style == DISPLAY_STYLE_ESPRIPT:
            return "espript_like"
        return "amino_acid_class"

    def default_pairwise_filename(self, extension):
        mode = self.analysis_mode_slug(self.current_pairwise_analysis_mode)

        seq_a = self.safe_filename_text(self.current_pairwise_seq_a_name)
        seq_b = self.safe_filename_text(self.current_pairwise_seq_b_name)

        return f"pairwise_{mode}_{seq_a}_vs_{seq_b}.{extension}"

    def default_msa_filename(self, extension):
        mode = self.analysis_mode_slug(self.current_msa_analysis_mode)
        style = self.display_style_slug(self.current_msa_display_style)

        seq_count = len(self.current_msa_records) if self.current_msa_records else len(self.get_checked_records_from_list(self.msa_sequence_list))
        return f"msa_{mode}_{style}_{seq_count}seq.{extension}"

    def default_matrix_filename(self):
        seq_count = len(self.current_matrix_records) if self.current_matrix_records else len(self.get_checked_records_from_list(self.matrix_sequence_list))
        return f"identity_matrix_{seq_count}seq.csv"

    def safe_filename_text(self, text):
        text = str(text).strip()

        if not text:
            return "sequence"

        safe = []

        for char in text:
            if char.isalnum() or char in ["-", "_", "."]:
                safe.append(char)
            else:
                safe.append("_")

        result = "".join(safe).strip("_")
        return result[:50] if result else "sequence"

    # ============================================================
    # Pairwise Tab
    # ============================================================

    def init_pairwise_tab(self):
        layout = QVBoxLayout()

        input_group = QGroupBox("Input Protein Sequences")
        input_layout = QGridLayout()

        self.seq_a_label = QLabel("Sequence A")
        self.seq_b_label = QLabel("Sequence B")

        self.seq_a_text = QTextEdit()
        self.seq_b_text = QTextEdit()

        self.seq_a_text.setPlaceholderText("Paste protein sequence A here, or load from FASTA...")
        self.seq_b_text.setPlaceholderText("Paste protein sequence B here, or load from FASTA...")

        mono_font = QFont("Courier New", 11)
        self.seq_a_text.setFont(mono_font)
        self.seq_b_text.setFont(mono_font)

        input_layout.addWidget(self.seq_a_label, 0, 0)
        input_layout.addWidget(self.seq_b_label, 0, 1)
        input_layout.addWidget(self.seq_a_text, 1, 0)
        input_layout.addWidget(self.seq_b_text, 1, 1)

        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        mode_layout = QHBoxLayout()

        self.pairwise_analysis_mode_combo = QComboBox()
        self.pairwise_analysis_mode_combo.addItem("Conservative analysis", ANALYSIS_MODE_CONSERVATIVE)
        self.pairwise_analysis_mode_combo.addItem("Mutation site comparison", ANALYSIS_MODE_MUTATION)

        mode_layout.addWidget(QLabel("Analysis mode:"))
        mode_layout.addWidget(self.pairwise_analysis_mode_combo)
        mode_layout.addStretch()

        layout.addLayout(mode_layout)

        button_layout = QHBoxLayout()

        self.load_fasta_button = QPushButton("Load FASTA(s)")
        self.load_structure_button = QPushButton("Load PDB/CIF")
        self.align_button = QPushButton("Align")
        self.clear_button = QPushButton("Clear")
        self.export_txt_button = QPushButton("Export TXT")
        self.export_html_button = QPushButton("Export HTML")
        self.export_pdf_button = QPushButton("Export PDF")

        self.load_fasta_button.clicked.connect(self.load_fasta)
        self.load_structure_button.clicked.connect(self.load_structure_files)
        self.align_button.clicked.connect(self.run_pairwise_alignment)
        self.clear_button.clicked.connect(self.clear_pairwise)
        self.export_txt_button.clicked.connect(self.export_pairwise_txt)
        self.export_html_button.clicked.connect(self.export_pairwise_html)
        self.export_pdf_button.clicked.connect(self.export_pairwise_pdf)

        button_layout.addWidget(self.load_fasta_button)
        button_layout.addWidget(self.load_structure_button)
        button_layout.addWidget(self.align_button)
        button_layout.addWidget(self.clear_button)
        button_layout.addWidget(self.export_txt_button)
        button_layout.addWidget(self.export_html_button)
        button_layout.addWidget(self.export_pdf_button)

        layout.addLayout(button_layout)

        fasta_group = QGroupBox("FASTA Records")
        fasta_layout = QHBoxLayout()

        self.record_a_combo = QComboBox()
        self.record_b_combo = QComboBox()

        self.record_a_combo.currentIndexChanged.connect(self.select_record_a)
        self.record_b_combo.currentIndexChanged.connect(self.select_record_b)

        fasta_layout.addWidget(QLabel("Use as Sequence A:"))
        fasta_layout.addWidget(self.record_a_combo)
        fasta_layout.addWidget(QLabel("Use as Sequence B:"))
        fasta_layout.addWidget(self.record_b_combo)

        fasta_group.setLayout(fasta_layout)
        layout.addWidget(fasta_group)

        output_group = QGroupBox("Pairwise Alignment Result")
        output_layout = QVBoxLayout()

        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        self.result_text.setFont(QFont("Courier New", 10))

        output_layout.addWidget(self.result_text)
        output_group.setLayout(output_layout)

        layout.addWidget(output_group)

        self.pairwise_tab.setLayout(layout)

    def select_record_a(self):
        index = self.record_a_combo.currentIndex()

        if self.fasta_records and 0 <= index < len(self.fasta_records):
            record = self.fasta_records[index]
            self.current_pairwise_seq_a_name = record.get("display_name", record.get("id", "SeqA"))
            self.seq_a_text.setPlainText(format_record_for_display(record))

    def select_record_b(self):
        index = self.record_b_combo.currentIndex()

        if self.fasta_records and 0 <= index < len(self.fasta_records):
            record = self.fasta_records[index]
            self.current_pairwise_seq_b_name = record.get("display_name", record.get("id", "SeqB"))
            self.seq_b_text.setPlainText(format_record_for_display(record))

    def extract_sequence_from_input(self, text: str) -> str:
        lines = text.strip().splitlines()
        seq_lines = []

        for line in lines:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                continue

            seq_lines.append(line)

        return "".join(seq_lines)

    def extract_name_from_input(self, text: str, default_name: str):
        lines = text.strip().splitlines()

        for line in lines:
            line = line.strip()

            if line.startswith(">"):
                name = line[1:].strip()
                return name if name else default_name

        return default_name

    def run_pairwise_alignment(self):
        seq_a_raw = self.seq_a_text.toPlainText()
        seq_b_raw = self.seq_b_text.toPlainText()

        seq_a = self.extract_sequence_from_input(seq_a_raw)
        seq_b = self.extract_sequence_from_input(seq_b_raw)

        self.current_pairwise_seq_a_name = self.extract_name_from_input(seq_a_raw, "SeqA")
        self.current_pairwise_seq_b_name = self.extract_name_from_input(seq_b_raw, "SeqB")

        analysis_mode = self.pairwise_analysis_mode_combo.currentData()
        self.current_pairwise_analysis_mode = analysis_mode

        try:
            result = align_pairwise(seq_a, seq_b, analysis_mode=analysis_mode)
            result_text = format_alignment_output(result, line_width=80)

            self.current_pairwise_result = result
            self.current_pairwise_text = result_text
            self.result_text.setPlainText(result_text)

        except Exception as e:
            QMessageBox.critical(self, "Alignment Error", str(e))

    def export_pairwise_txt(self):
        if not self.current_pairwise_text.strip():
            QMessageBox.warning(self, "No Result", "Please run pairwise alignment first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Pairwise Alignment TXT",
            self.default_pairwise_filename("txt"),
            "Text Files (*.txt);;All Files (*)",
        )

        if not file_path:
            return

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(self.current_pairwise_text)

            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def export_pairwise_html(self):
        if self.current_pairwise_result is None:
            QMessageBox.warning(self, "No Result", "Please run pairwise alignment first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Pairwise Alignment HTML",
            self.default_pairwise_filename("html"),
            "HTML Files (*.html);;All Files (*)",
        )

        if not file_path:
            return

        try:
            html_text = pairwise_result_to_html(
                self.current_pairwise_result,
                seq_a_name=self.current_pairwise_seq_a_name,
                seq_b_name=self.current_pairwise_seq_b_name,
            )

            with open(file_path, "w", encoding="utf-8") as f:
                f.write(html_text)

            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def export_pairwise_pdf(self):
        if self.current_pairwise_result is None:
            QMessageBox.warning(self, "No Result", "Please run pairwise alignment first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Pairwise Alignment PDF",
            self.default_pairwise_filename("pdf"),
            "PDF Files (*.pdf);;All Files (*)",
        )

        if not file_path:
            return

        try:
            export_pairwise_pdf(
                self.current_pairwise_result,
                file_path,
                seq_a_name=self.current_pairwise_seq_a_name,
                seq_b_name=self.current_pairwise_seq_b_name,
            )

            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "PDF Export Error", str(e))

    def clear_pairwise(self):
        self.seq_a_text.clear()
        self.seq_b_text.clear()
        self.result_text.clear()
        self.current_pairwise_result = None
        self.current_pairwise_text = ""

    # ============================================================
    # MSA Tab
    # ============================================================

    def init_msa_tab(self):
        layout = QVBoxLayout()

        top_group = QGroupBox("Multiple Sequence Alignment")
        top_layout = QVBoxLayout()

        instruction_label = QLabel(
            "Load FASTA/PDB/CIF files. Tick checkboxes to choose sequences for comparison or deletion. "
            "The grey highlighted row is only for Rename / Move Up / Move Down. Double-click a row to rename it."
        )
        top_layout.addWidget(instruction_label)

        self.msa_sequence_list = QListWidget()
        self.msa_sequence_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.msa_sequence_list.itemDoubleClicked.connect(
            lambda item: self.rename_record_from_list(self.msa_sequence_list)
        )
        top_layout.addWidget(self.msa_sequence_list)

        edit_layout = QHBoxLayout()

        self.load_msa_fasta_button = QPushButton("Load FASTA(s)")
        self.load_msa_structure_button = QPushButton("Load PDB/CIF")
        self.rename_msa_button = QPushButton("Rename Selected")
        self.delete_msa_button = QPushButton("Delete Checked")
        self.move_msa_up_button = QPushButton("Move Up")
        self.move_msa_down_button = QPushButton("Move Down")

        self.load_msa_fasta_button.clicked.connect(self.load_fasta)
        self.load_msa_structure_button.clicked.connect(self.load_structure_files)
        self.rename_msa_button.clicked.connect(lambda: self.rename_record_from_list(self.msa_sequence_list))
        self.delete_msa_button.clicked.connect(lambda: self.delete_checked_records_from_list(self.msa_sequence_list))
        self.move_msa_up_button.clicked.connect(lambda: self.move_record_from_list(self.msa_sequence_list, -1))
        self.move_msa_down_button.clicked.connect(lambda: self.move_record_from_list(self.msa_sequence_list, 1))

        edit_layout.addWidget(self.load_msa_fasta_button)
        edit_layout.addWidget(self.load_msa_structure_button)
        edit_layout.addWidget(self.rename_msa_button)
        edit_layout.addWidget(self.delete_msa_button)
        edit_layout.addWidget(self.move_msa_up_button)
        edit_layout.addWidget(self.move_msa_down_button)

        top_layout.addLayout(edit_layout)

        msa_analysis_mode_layout = QHBoxLayout()

        self.msa_analysis_mode_combo = QComboBox()
        self.msa_analysis_mode_combo.addItem("Conservative analysis", ANALYSIS_MODE_CONSERVATIVE)
        self.msa_analysis_mode_combo.addItem("Mutation site comparison", ANALYSIS_MODE_MUTATION)

        msa_analysis_mode_layout.addWidget(QLabel("Analysis mode:"))
        msa_analysis_mode_layout.addWidget(self.msa_analysis_mode_combo)
        msa_analysis_mode_layout.addStretch()

        top_layout.addLayout(msa_analysis_mode_layout)

        msa_style_layout = QHBoxLayout()

        self.msa_display_style_combo = QComboBox()
        self.msa_display_style_combo.addItem("Amino acid class style", DISPLAY_STYLE_AA_CLASS)
        self.msa_display_style_combo.addItem("ESPript-like conservation style", DISPLAY_STYLE_ESPRIPT)

        msa_style_layout.addWidget(QLabel("MSA display style:"))
        msa_style_layout.addWidget(self.msa_display_style_combo)
        msa_style_layout.addStretch()

        top_layout.addLayout(msa_style_layout)

        color_group = QGroupBox("ESPript-like color customization")
        color_layout = QHBoxLayout()

        self.espript_conserved_bg_button = QPushButton()
        self.espript_conserved_text_button = QPushButton()
        self.espript_similar_text_button = QPushButton()
        self.espript_similar_border_button = QPushButton()
        self.espript_gap_bg_button = QPushButton()
        self.espript_reset_button = QPushButton("Reset ESPript colors")

        self.espript_conserved_bg_button.clicked.connect(
            lambda: self.set_espript_color("conserved_bg", "Choose conserved background color")
        )
        self.espript_conserved_text_button.clicked.connect(
            lambda: self.set_espript_color("conserved_text", "Choose conserved text color")
        )
        self.espript_similar_text_button.clicked.connect(
            lambda: self.set_espript_color("similar_text", "Choose similar residue text color")
        )
        self.espript_similar_border_button.clicked.connect(
            lambda: self.set_espript_color("similar_border", "Choose similar column border color")
        )
        self.espript_gap_bg_button.clicked.connect(
            lambda: self.set_espript_color("gap_bg", "Choose gap background color")
        )
        self.espript_reset_button.clicked.connect(self.reset_espript_colors)

        color_layout.addWidget(self.espript_conserved_bg_button)
        color_layout.addWidget(self.espript_conserved_text_button)
        color_layout.addWidget(self.espript_similar_text_button)
        color_layout.addWidget(self.espript_similar_border_button)
        color_layout.addWidget(self.espript_gap_bg_button)
        color_layout.addWidget(self.espript_reset_button)

        color_group.setLayout(color_layout)
        top_layout.addWidget(color_group)

        self.update_espript_color_button_labels()

        controls_layout = QHBoxLayout()

        self.mafft_mode_combo = QComboBox()
        self.mafft_mode_combo.addItems(["auto", "localpair", "globalpair", "genafpair"])

        self.run_mafft_button = QPushButton("Run MAFFT")
        self.export_msa_txt_button = QPushButton("Export MSA TXT")
        self.export_msa_fasta_button = QPushButton("Export MSA FASTA")
        self.export_msa_html_button = QPushButton("Export MSA HTML")
        self.export_msa_pdf_button = QPushButton("Export MSA PDF")

        self.run_mafft_button.clicked.connect(self.run_msa_mafft)
        self.export_msa_txt_button.clicked.connect(self.export_msa_txt)
        self.export_msa_fasta_button.clicked.connect(self.export_msa_fasta)
        self.export_msa_html_button.clicked.connect(self.export_msa_html)
        self.export_msa_pdf_button.clicked.connect(self.export_msa_pdf)

        controls_layout.addWidget(QLabel("MAFFT mode:"))
        controls_layout.addWidget(self.mafft_mode_combo)
        controls_layout.addWidget(self.run_mafft_button)
        controls_layout.addWidget(self.export_msa_txt_button)
        controls_layout.addWidget(self.export_msa_fasta_button)
        controls_layout.addWidget(self.export_msa_html_button)
        controls_layout.addWidget(self.export_msa_pdf_button)

        top_layout.addLayout(controls_layout)

        top_group.setLayout(top_layout)
        layout.addWidget(top_group)

        output_group = QGroupBox("MSA Result")
        output_layout = QVBoxLayout()

        self.msa_result_text = QTextEdit()
        self.msa_result_text.setReadOnly(True)
        self.msa_result_text.setFont(QFont("Courier New", 10))

        output_layout.addWidget(self.msa_result_text)
        output_group.setLayout(output_layout)

        layout.addWidget(output_group)

        self.msa_tab.setLayout(layout)

    def run_msa_mafft(self):
        records = self.get_checked_records_from_list(self.msa_sequence_list)

        if len(records) < 2:
            QMessageBox.warning(
                self,
                "Not Enough Sequences",
                "Please load a multi-sequence FASTA/PDB/CIF file and tick at least two checkboxes.",
            )
            return

        if not check_mafft_available():
            QMessageBox.critical(
                self,
                "MAFFT Not Detected",
                "MAFFT was not found.\n\n"
                "Please install it first:\n\n"
                "brew install mafft\n\n"
                "Then restart this app.",
            )
            return

        mafft_mode = self.mafft_mode_combo.currentText()
        analysis_mode = self.msa_analysis_mode_combo.currentData()
        display_style = self.msa_display_style_combo.currentData()

        self.current_msa_analysis_mode = analysis_mode
        self.current_msa_display_style = display_style

        try:
            aligned_records, stderr_text = run_mafft(records, mafft_mode=mafft_mode)
            msa_text = format_msa_text(
                aligned_records,
                analysis_mode=analysis_mode,
            )

            self.current_msa_records = aligned_records
            self.current_msa_text = msa_text
            self.msa_result_text.setPlainText(msa_text)

        except Exception as e:
            QMessageBox.critical(self, "MAFFT Error", str(e))

    def export_msa_txt(self):
        if not self.current_msa_text.strip():
            QMessageBox.warning(self, "No MSA Result", "Please run MAFFT first.")
            return

        self.current_msa_display_style = self.msa_display_style_combo.currentData()

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save MSA TXT",
            self.default_msa_filename("txt"),
            "Text Files (*.txt);;All Files (*)",
        )

        if not file_path:
            return

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(self.current_msa_text)

            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def export_msa_fasta(self):
        if not self.current_msa_records:
            QMessageBox.warning(self, "No MSA Result", "Please run MAFFT first.")
            return

        self.current_msa_display_style = self.msa_display_style_combo.currentData()

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save MSA FASTA",
            self.default_msa_filename("fasta"),
            "FASTA Files (*.fasta *.fa);;All Files (*)",
        )

        if not file_path:
            return

        try:
            write_fasta_file(self.current_msa_records, file_path)
            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def export_msa_html(self):
        if not self.current_msa_records:
            QMessageBox.warning(self, "No MSA Result", "Please run MAFFT first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save MSA HTML",
            self.default_msa_filename("html"),
            "HTML Files (*.html);;All Files (*)",
        )

        if not file_path:
            return

        try:
            self.current_msa_display_style = self.msa_display_style_combo.currentData()
            html_text = msa_to_html(
                self.current_msa_records,
                analysis_mode=self.current_msa_analysis_mode,
                display_style=self.current_msa_display_style,
                color_config=self.espript_colors,
            )

            with open(file_path, "w", encoding="utf-8") as f:
                f.write(html_text)

            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def export_msa_pdf(self):
        if not self.current_msa_records:
            QMessageBox.warning(self, "No MSA Result", "Please run MAFFT first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save MSA PDF",
            self.default_msa_filename("pdf"),
            "PDF Files (*.pdf);;All Files (*)",
        )

        if not file_path:
            return

        try:
            self.current_msa_display_style = self.msa_display_style_combo.currentData()
            export_msa_pdf(
                self.current_msa_records,
                file_path,
                analysis_mode=self.current_msa_analysis_mode,
                display_style=self.current_msa_display_style,
                color_config=self.espript_colors,
            )
            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "PDF Export Error", str(e))

    # ============================================================
    # Domain Tab
    # ============================================================

    def init_domain_tab(self):
        layout = QVBoxLayout()

        input_group = QGroupBox("Domain / Segment Extraction")
        input_layout = QVBoxLayout()

        info_label = QLabel(
            "Two extraction modes are available:\n"
            "1. Shared domain definition: one residue range is applied to all loaded sequences.\n"
            "2. Per-sequence segment table: set a different start/end range for each sequence, then extract those segments for alignment."
        )
        input_layout.addWidget(info_label)

        self.domain_text = QTextEdit()
        self.domain_text.setFont(QFont("Courier New", 10))
        self.domain_text.setPlaceholderText(
            "GTD: 1-540\nCPD: 541-758\nDRBD: 830-1756\nCST: 1756-1776"
        )
        self.domain_text.setPlainText(
            "GTD: 1-540\nCPD: 541-758\nDRBD: 830-1756\nCST: 1756-1776"
        )

        input_layout.addWidget(self.domain_text)

        button_layout = QHBoxLayout()

        self.load_domain_fasta_button = QPushButton("Load FASTA(s)")
        self.load_domain_structure_button = QPushButton("Load PDB/CIF")
        self.extract_domains_button = QPushButton("Extract Shared Domains")
        self.export_all_domains_button = QPushButton("Export Shared Domain FASTA")
        self.export_selected_domain_button = QPushButton("Export Selected Shared Domain")
        self.domain_select_combo = QComboBox()

        self.load_domain_fasta_button.clicked.connect(self.load_fasta)
        self.load_domain_structure_button.clicked.connect(self.load_structure_files)
        self.extract_domains_button.clicked.connect(self.extract_domains)
        self.export_all_domains_button.clicked.connect(self.export_all_domains)
        self.export_selected_domain_button.clicked.connect(self.export_selected_domain)

        button_layout.addWidget(self.load_domain_fasta_button)
        button_layout.addWidget(self.load_domain_structure_button)
        button_layout.addWidget(self.extract_domains_button)
        button_layout.addWidget(QLabel("Selected shared domain:"))
        button_layout.addWidget(self.domain_select_combo)
        button_layout.addWidget(self.export_selected_domain_button)
        button_layout.addWidget(self.export_all_domains_button)

        input_layout.addLayout(button_layout)

        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        segment_group = QGroupBox("Per-sequence segment extraction for variable domain positions")
        segment_layout = QVBoxLayout()

        segment_info = QLabel(
            "Use this table when the same domain is conserved but located at different residue positions in different sequences. "
            "Tick rows, set Domain name / Start / End for each sequence, then extract selected segments."
        )
        segment_layout.addWidget(segment_info)

        self.domain_segment_table = QTableWidget()
        self.domain_segment_table.setColumnCount(6)
        self.domain_segment_table.setHorizontalHeaderLabels(
            ["Use", "Sequence name", "Length", "Domain name", "Start", "End"]
        )
        self.domain_segment_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.domain_segment_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        segment_layout.addWidget(self.domain_segment_table)

        segment_button_layout = QHBoxLayout()

        self.refresh_segment_table_button = QPushButton("Refresh Table from Loaded Sequences")
        self.extract_segments_button = QPushButton("Extract Checked Segments")
        self.add_segments_to_msa_button = QPushButton("Add Extracted Segments to Sequence List")
        self.export_segments_fasta_button = QPushButton("Export Extracted Segments FASTA")
        self.clear_segments_button = QPushButton("Clear Extracted Segments")

        self.refresh_segment_table_button.clicked.connect(self.refresh_domain_segment_table)
        self.extract_segments_button.clicked.connect(self.extract_custom_segments)
        self.add_segments_to_msa_button.clicked.connect(self.add_extracted_segments_to_sequence_list)
        self.export_segments_fasta_button.clicked.connect(self.export_custom_segments_fasta)
        self.clear_segments_button.clicked.connect(self.clear_custom_segments)

        segment_button_layout.addWidget(self.refresh_segment_table_button)
        segment_button_layout.addWidget(self.extract_segments_button)
        segment_button_layout.addWidget(self.add_segments_to_msa_button)
        segment_button_layout.addWidget(self.export_segments_fasta_button)
        segment_button_layout.addWidget(self.clear_segments_button)

        segment_layout.addLayout(segment_button_layout)

        segment_group.setLayout(segment_layout)
        layout.addWidget(segment_group)

        output_group = QGroupBox("Domain / Segment Extraction Result")
        output_layout = QVBoxLayout()

        self.domain_result_text = QTextEdit()
        self.domain_result_text.setReadOnly(True)
        self.domain_result_text.setFont(QFont("Courier New", 10))

        output_layout.addWidget(self.domain_result_text)
        output_group.setLayout(output_layout)

        layout.addWidget(output_group)

        self.domain_tab.setLayout(layout)

    def extract_domains(self):
        if not self.fasta_records:
            QMessageBox.warning(self, "No FASTA", "Please load a FASTA file first.")
            return

        try:
            domains = parse_domain_text(self.domain_text.toPlainText())
            domain_outputs = extract_domains_from_records(self.fasta_records, domains)

            self.current_domain_outputs = domain_outputs

            self.domain_select_combo.clear()

            for domain_name in domain_outputs.keys():
                self.domain_select_combo.addItem(domain_name)

            lines = []
            lines.append("Domain Extraction Result")
            lines.append("=" * 60)
            lines.append("")

            for domain_name, records in domain_outputs.items():
                lines.append(f"[{domain_name}]")

                for record in records:
                    name = record.get("display_name", record.get("id", "sequence"))
                    lines.append(f"{name}    length: {len(record['sequence'])}")

                lines.append("")

            self.domain_result_text.setPlainText("\n".join(lines))

        except Exception as e:
            QMessageBox.critical(self, "Domain Extraction Error", str(e))

    def export_selected_domain(self):
        if not self.current_domain_outputs:
            QMessageBox.warning(self, "No Domain Output", "Please extract domains first.")
            return

        domain_name = self.domain_select_combo.currentText()

        if not domain_name:
            QMessageBox.warning(self, "No Domain Selected", "Please select a domain.")
            return

        records = self.current_domain_outputs.get(domain_name, [])

        if not records:
            QMessageBox.warning(self, "No Records", "No records found for selected domain.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Selected Domain FASTA",
            f"{domain_name}.fasta",
            "FASTA Files (*.fasta *.fa);;All Files (*)",
        )

        if not file_path:
            return

        try:
            write_fasta_file(records, file_path)
            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def export_all_domains(self):
        if not self.current_domain_outputs:
            QMessageBox.warning(self, "No Domain Output", "Please extract domains first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save All Domain FASTA",
            "all_domains.fasta",
            "FASTA Files (*.fasta *.fa);;All Files (*)",
        )

        if not file_path:
            return

        try:
            all_records = []

            for domain_name, records in self.current_domain_outputs.items():
                all_records.extend(records)

            write_fasta_file(all_records, file_path)
            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))


    # ============================================================
    # Per-sequence custom segment extraction
    # ============================================================

    def refresh_domain_segment_table(self):
        """
        Populate the per-sequence segment table from the currently loaded sequence list.
        """
        if not hasattr(self, "domain_segment_table"):
            return

        self.domain_segment_table.setRowCount(len(self.fasta_records))

        for row, record in enumerate(self.fasta_records):
            name = record.get("display_name", record.get("id", f"sequence_{row + 1}"))
            sequence = record.get("sequence", "")
            length = len(sequence)

            use_item = QTableWidgetItem()
            use_item.setFlags(
                Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsSelectable
            )
            use_item.setCheckState(Qt.CheckState.Checked)
            self.domain_segment_table.setItem(row, 0, use_item)

            name_item = QTableWidgetItem(name)
            name_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            self.domain_segment_table.setItem(row, 1, name_item)

            length_item = QTableWidgetItem(str(length))
            length_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            length_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.domain_segment_table.setItem(row, 2, length_item)

            domain_item = QTableWidgetItem("domain")
            self.domain_segment_table.setItem(row, 3, domain_item)

            start_item = QTableWidgetItem("1")
            self.domain_segment_table.setItem(row, 4, start_item)

            end_item = QTableWidgetItem(str(length))
            self.domain_segment_table.setItem(row, 5, end_item)

        self.domain_segment_table.resizeColumnsToContents()

    def extract_custom_segments(self):
        """
        Extract independently defined regions from each checked sequence row.
        Each row can have its own domain name, start, and end.
        """
        if not self.fasta_records:
            QMessageBox.warning(self, "No Sequences", "Please load FASTA/PDB/CIF sequences first.")
            return

        if self.domain_segment_table.rowCount() == 0:
            self.refresh_domain_segment_table()

        extracted = []
        errors = []

        for row in range(self.domain_segment_table.rowCount()):
            use_item = self.domain_segment_table.item(row, 0)

            if use_item is None or use_item.checkState() != Qt.CheckState.Checked:
                continue

            if row >= len(self.fasta_records):
                continue

            record = self.fasta_records[row]
            source_name = record.get("display_name", record.get("id", f"sequence_{row + 1}"))
            sequence = record.get("sequence", "")
            length = len(sequence)

            domain_item = self.domain_segment_table.item(row, 3)
            start_item = self.domain_segment_table.item(row, 4)
            end_item = self.domain_segment_table.item(row, 5)

            domain_name = domain_item.text().strip() if domain_item is not None else "domain"

            if not domain_name:
                domain_name = "domain"

            try:
                start = int(start_item.text().strip()) if start_item is not None else 1
                end = int(end_item.text().strip()) if end_item is not None else length
            except ValueError:
                errors.append(f"{source_name}: start/end must be integers.")
                continue

            if start < 1 or end < start or end > length:
                errors.append(
                    f"{source_name}: invalid range {start}-{end}; valid range is 1-{length}."
                )
                continue

            segment = sequence[start - 1:end]
            segment_name = f"{source_name}_{domain_name}_{start}-{end}"

            extracted.append({
                "id": segment_name,
                "display_name": segment_name,
                "description": (
                    f"{segment_name} | extracted from {source_name} | "
                    f"{domain_name}:{start}-{end} | length {len(segment)}"
                ),
                "sequence": segment,
                "source_record": source_name,
                "domain_name": domain_name,
                "source_range": f"{start}-{end}",
            })

        if errors:
            QMessageBox.warning(
                self,
                "Some segments were not extracted",
                "\n".join(errors[:20]) + (f"\n... and {len(errors) - 20} more" if len(errors) > 20 else ""),
            )

        if not extracted:
            QMessageBox.warning(
                self,
                "No Segments Extracted",
                "No valid checked segment was extracted. Please check the Use boxes and start/end values.",
            )
            return

        self.current_custom_domain_records = extracted
        self.show_custom_segment_result(extracted)

        QMessageBox.information(
            self,
            "Segments Extracted",
            f"Extracted {len(extracted)} segment(s).\n\n"
            "You can now add them to the sequence list for alignment or export them as FASTA.",
        )

    def show_custom_segment_result(self, records):
        lines = []
        lines.append("Per-sequence Segment Extraction Result")
        lines.append("=" * 70)
        lines.append("")

        for record in records:
            name = record.get("display_name", record.get("id", "segment"))
            source = record.get("source_record", "-")
            domain_name = record.get("domain_name", "-")
            source_range = record.get("source_range", "-")
            length = len(record.get("sequence", ""))

            lines.append(f"{name}")
            lines.append(f"  Source: {source}")
            lines.append(f"  Domain: {domain_name}")
            lines.append(f"  Range:  {source_range}")
            lines.append(f"  Length: {length}")
            lines.append("")

        self.domain_result_text.setPlainText("\n".join(lines))

    def add_extracted_segments_to_sequence_list(self):
        """
        Append extracted custom segments to the global sequence list so the user can
        select them for pairwise alignment, MSA, or identity matrix.
        """
        if not self.current_custom_domain_records:
            QMessageBox.warning(
                self,
                "No Extracted Segments",
                "Please click Extract Checked Segments first.",
            )
            return

        self.fasta_records.extend([dict(record) for record in self.current_custom_domain_records])
        self.refresh_all_fasta_controls()
        self.fill_pairwise_inputs_after_import()

        QMessageBox.information(
            self,
            "Added to Sequence List",
            f"Added {len(self.current_custom_domain_records)} extracted segment(s) to the sequence list.\n\n"
            "Go to Multiple Alignment and tick these segment records for MAFFT alignment.",
        )

    def export_custom_segments_fasta(self):
        if not self.current_custom_domain_records:
            QMessageBox.warning(
                self,
                "No Extracted Segments",
                "Please click Extract Checked Segments first.",
            )
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Extracted Segment FASTA",
            "custom_extracted_segments.fasta",
            "FASTA Files (*.fasta *.fa);;All Files (*)",
        )

        if not file_path:
            return

        try:
            write_fasta_file(self.current_custom_domain_records, file_path)
            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def clear_custom_segments(self):
        self.current_custom_domain_records = []
        self.domain_result_text.clear()
        QMessageBox.information(self, "Cleared", "Extracted segment records were cleared.")

    # ============================================================
    # Matrix Tab
    # ============================================================

    def init_matrix_tab(self):
        layout = QVBoxLayout()

        top_group = QGroupBox("Pairwise Identity Matrix")
        top_layout = QVBoxLayout()

        info_label = QLabel(
            "Load FASTA/PDB/CIF files. Tick checkboxes to choose sequences for matrix calculation or deletion. "
            "The grey highlighted row is only for Rename / Move Up / Move Down. Double-click a row to rename it."
        )
        top_layout.addWidget(info_label)

        self.matrix_sequence_list = QListWidget()
        self.matrix_sequence_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.matrix_sequence_list.itemDoubleClicked.connect(
            lambda item: self.rename_record_from_list(self.matrix_sequence_list)
        )
        top_layout.addWidget(self.matrix_sequence_list)

        edit_layout = QHBoxLayout()

        self.load_matrix_fasta_button = QPushButton("Load FASTA(s)")
        self.load_matrix_structure_button = QPushButton("Load PDB/CIF")
        self.rename_matrix_button = QPushButton("Rename Selected")
        self.delete_matrix_button = QPushButton("Delete Checked")
        self.move_matrix_up_button = QPushButton("Move Up")
        self.move_matrix_down_button = QPushButton("Move Down")

        self.load_matrix_fasta_button.clicked.connect(self.load_fasta)
        self.load_matrix_structure_button.clicked.connect(self.load_structure_files)
        self.rename_matrix_button.clicked.connect(lambda: self.rename_record_from_list(self.matrix_sequence_list))
        self.delete_matrix_button.clicked.connect(lambda: self.delete_checked_records_from_list(self.matrix_sequence_list))
        self.move_matrix_up_button.clicked.connect(lambda: self.move_record_from_list(self.matrix_sequence_list, -1))
        self.move_matrix_down_button.clicked.connect(lambda: self.move_record_from_list(self.matrix_sequence_list, 1))

        edit_layout.addWidget(self.load_matrix_fasta_button)
        edit_layout.addWidget(self.load_matrix_structure_button)
        edit_layout.addWidget(self.rename_matrix_button)
        edit_layout.addWidget(self.delete_matrix_button)
        edit_layout.addWidget(self.move_matrix_up_button)
        edit_layout.addWidget(self.move_matrix_down_button)

        top_layout.addLayout(edit_layout)

        button_layout = QHBoxLayout()

        self.calculate_matrix_button = QPushButton("Calculate Identity Matrix")
        self.export_matrix_csv_button = QPushButton("Export Matrix CSV")

        self.calculate_matrix_button.clicked.connect(self.calculate_identity_matrix)
        self.export_matrix_csv_button.clicked.connect(self.export_identity_matrix_csv)

        button_layout.addWidget(self.calculate_matrix_button)
        button_layout.addWidget(self.export_matrix_csv_button)

        top_layout.addLayout(button_layout)

        top_group.setLayout(top_layout)
        layout.addWidget(top_group)

        self.identity_matrix_table = QTableWidget()
        layout.addWidget(self.identity_matrix_table)

        self.matrix_tab.setLayout(layout)

    def calculate_identity_matrix(self):
        records = self.get_checked_records_from_list(self.matrix_sequence_list)

        if len(records) < 2:
            QMessageBox.warning(
                self,
                "Not Enough Sequences",
                "Please load a FASTA/PDB/CIF file and tick at least two checkboxes.",
            )
            return

        try:
            matrix = calculate_pairwise_identity_matrix(records)

            self.current_matrix_records = records
            self.current_identity_matrix = matrix

            n = len(records)
            self.identity_matrix_table.clear()
            self.identity_matrix_table.setRowCount(n)
            self.identity_matrix_table.setColumnCount(n)

            labels = [record.get("display_name", record.get("id", "")) for record in records]

            self.identity_matrix_table.setHorizontalHeaderLabels(labels)
            self.identity_matrix_table.setVerticalHeaderLabels(labels)

            for i in range(n):
                for j in range(n):
                    item = QTableWidgetItem(f"{matrix[i][j]:.2f}")
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.identity_matrix_table.setItem(i, j, item)

            self.identity_matrix_table.resizeColumnsToContents()

        except Exception as e:
            QMessageBox.critical(self, "Matrix Error", str(e))

    def export_identity_matrix_csv(self):
        if not self.current_identity_matrix or not self.current_matrix_records:
            QMessageBox.warning(self, "No Matrix", "Please calculate identity matrix first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Identity Matrix CSV",
            self.default_matrix_filename(),
            "CSV Files (*.csv);;All Files (*)",
        )

        if not file_path:
            return

        try:
            labels = [
                record.get("display_name", record.get("id", ""))
                for record in self.current_matrix_records
            ]

            with open(file_path, "w", encoding="utf-8") as f:
                f.write("," + ",".join(labels) + "\n")

                for label, row in zip(labels, self.current_identity_matrix):
                    row_text = ",".join([f"{value:.2f}" for value in row])
                    f.write(f"{label},{row_text}\n")

            QMessageBox.information(self, "Export Successful", f"Saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))


def main():
    app = QApplication(sys.argv)
    window = ProteinAlignmentApp()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
