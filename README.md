# Protein Alignment Tool v1.0.0

Protein Alignment Tool is a local graphical application for protein sequence alignment, sequence comparison, domain and segment extraction, mutation-site comparison, and publication-ready HTML/PDF export.

This release provides both macOS and Windows usage options. macOS users can use the packaged DMG installer. Windows users can use the Windows source release package together with the official Windows MAFFT all-in-one package.

## Download

### macOS

Download: `Protein_Alignment_Tool_macOS.dmg`

Installation:

1. Open the DMG file.
2. Drag `Protein Alignment Tool.app` to the `Applications` folder.
3. Launch the app from `Applications`.
4. If macOS shows a security warning, right-click the app and choose **Open** for the first launch.

The macOS package includes the required runtime environment and MAFFT.

### Windows

Download: `Protein_Alignment_Tool_Windows_source_release.zip`

Windows users also need to install the official Windows MAFFT all-in-one package for multiple sequence alignment.

Recommended MAFFT location: `C:\Tools\mafft-win\mafft.bat`

After setup, start the application by double-clicking `run_app_windows.bat`.

For detailed Windows setup instructions, see `README_Windows.md`.

## Main Features

### Sequence input

- Load one or multiple FASTA files.
- Drag and drop FASTA files into the application window.
- Load PDB, CIF, and mmCIF structure files.
- Extract protein sequences from individual chains in PDB, CIF, and mmCIF files.
- Convert extracted structure-chain sequences into FASTA records for downstream comparison.

### Pairwise alignment

- Perform pairwise protein sequence alignment.
- Display aligned sequences with residue-level comparison.
- Support conservative analysis mode.
- Support mutation-site comparison mode.
- Export pairwise alignment results as HTML or PDF.

### Multiple sequence alignment

- Perform multiple sequence alignment through MAFFT.
- Preserve user-defined sequence names.
- Allow users to manually adjust sequence order before alignment.
- Support multiple display modes for aligned residues.
- Export MSA results as HTML, PDF, TXT, or FASTA.

### Analysis modes

This release includes two major comparison modes.

Conservative analysis highlights conserved and chemically similar residues. It is useful for comparing related protein families, conserved motifs, or homologous domains.

Mutation site comparison is designed for point-mutant or variant comparison. Identical residues remain unmarked, different residues are highlighted with red background and white text, and gaps are displayed as unmarked gap characters.

### ESPript-like display style

- Includes an ESPript-like conservation visualization mode.
- Conserved residues and similar residues are highlighted in a publication-oriented style.
- Colors are softer than high-saturation red schemes.
- Users can adjust display colors for ESPript-like outputs.

### Domain and segment extraction

- Define shared domain boundaries using residue numbers.
- Extract domains from all loaded sequences.
- Support per-sequence segment extraction when conserved domains occur at different sequence positions.
- Assign custom domain names to extracted segments.
- Add extracted segments back to the sequence list for downstream alignment.
- Export selected or all extracted domain and segment FASTA files.

### Sequence management

- Rename each sequence manually.
- Delete selected sequences.
- Reorder sequences manually.
- Select sequences using checkboxes for alignment, deletion, and matrix calculation.
- Keep row highlighting separate from checkbox-based selection.

### Identity matrix

- Calculate pairwise sequence identity matrix.
- Select which sequences to include.
- Export identity matrix as CSV.

### Export

Supported output formats include HTML, PDF, FASTA, TXT, and CSV. The HTML and PDF outputs are designed for convenient inspection, sharing, and figure preparation.

## Platform Notes

### macOS

The macOS release is distributed as a DMG package. The macOS application bundle includes Python runtime, PyQt, Biopython, ReportLab, MAFFT, and application source code/resources.

No separate Python or MAFFT installation is required for macOS users.

If macOS reports that the application is from an unidentified developer, right-click the app and select **Open** for the first launch.

If macOS reports that the app is damaged, remove the quarantine attribute manually:

    xattr -cr "/Applications/Protein Alignment Tool.app"

Then right-click the app and choose **Open**.

### Windows

The Windows release is distributed as a source package.

Windows users need Miniconda or Anaconda, a conda environment named `protein_align`, PyQt5, Biopython, ReportLab, and the official Windows MAFFT all-in-one package.

Recommended Windows MAFFT path: `C:\Tools\mafft-win\mafft.bat`

Important notes for Windows users:

- Do not mix files from different software versions.
- Use the full Windows source release package as provided.
- Native Windows conda often cannot install MAFFT from Bioconda.
- Use the official Windows MAFFT all-in-one package instead.
- If MAFFT is not detected, Multiple Alignment will be unavailable, but other functions can still be used.

## Recommended Windows Setup Summary

1. Install Miniconda or Anaconda.
2. Download the official Windows MAFFT all-in-one package.
3. Extract MAFFT to `C:\Tools\mafft-win`.
4. Confirm that `C:\Tools\mafft-win\mafft.bat` exists.
5. Extract `Protein_Alignment_Tool_Windows_source_release.zip`.
6. Double-click `run_app_windows.bat`.
7. Confirm that the app shows `MAFFT detected`.

Then the Multiple Alignment function is ready to use.

## Known Limitations

- Windows users currently use the source release rather than a standalone `.exe` installer.
- Windows multiple sequence alignment requires separately installing the official Windows MAFFT all-in-one package.
- macOS builds may be architecture-dependent. A DMG built on Apple Silicon is recommended for Apple Silicon Macs.
- The app is not notarized with an Apple Developer ID in this release, so macOS users may need to right-click and choose **Open** on first launch.
- Very large sequence sets or very long proteins may require additional time for MAFFT alignment and PDF/HTML export.

## Files Included in This Release

Recommended release assets:

- `Protein_Alignment_Tool_macOS.dmg`
- `Protein_Alignment_Tool_Windows_source_release.zip`
- `README_macOS_CN.md`
- `README_macOS_EN.md`
- `README_Windows.md`

Windows source package contents include:

- `main.py`
- `alignment_core.py`
- `fasta_utils.py`
- `mafft_utils.py`
- `html_exporter.py`
- `pdf_exporter.py`
- `structure_utils.py`
- `requirements_windows.txt`
- `run_app_windows.bat`
- `run_app_windows_conda.bat`
- `README_Windows_CN.md`
- `COPYRIGHT_NOTICE.txt`
- `VERSION.txt`
- `example_data/example_sequences.fasta`

## Citation and Third-party Tools

Protein Alignment Tool uses or depends on several third-party tools and libraries. If results generated by this software are used in research, reports, or publications, please cite the relevant tools according to the functions used.

### MAFFT

Used for multiple sequence alignment.

Katoh, K., Misawa, K., Kuma, K. and Miyata, T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. *Nucleic Acids Research* 30, 3059–3066, 2002.

Katoh, K. and Standley, D.M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30, 772–780, 2013.

### Biopython

Used for biological sequence and structure file parsing.

Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25, 1422–1423, 2009.

### PyQt / Qt

Used for the graphical user interface.

### ReportLab

Used for PDF export.

### Python

Used as the main software runtime.

## Copyright

Copyright © 2026 yangyu. All rights reserved.

Protein Alignment Tool, including the software interface, workflow integration, documentation, and release materials, is owned by **yangyu**.

The software may be used for research, teaching, and non-commercial purposes.

Commercial redistribution, repackaging, or release under another name requires permission from the copyright holder.

## Version

Release version: `v1.0.0`

Initial public release of Protein Alignment Tool for macOS and Windows users.
