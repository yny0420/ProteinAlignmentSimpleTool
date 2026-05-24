#!/usr/bin/env python3
"""
Build macOS .app with bundled MAFFT.

Run from the project folder:

    cd /Users/yangyu/Desktop/Projects/protein_alignment_app
    source .venv/bin/activate
    python build_macos_bundled_mafft.py

Output:

    dist/Protein Alignment Tool.app

Notes:
- This script is designed for macOS.
- It bundles the Homebrew MAFFT installation into the app.
- Build on Apple Silicon for Apple Silicon users.
- Build on Intel Mac for Intel users.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


APP_NAME = "Protein Alignment Tool"


def run(cmd, *, check=True):
    print("\n>", " ".join(str(x) for x in cmd))
    return subprocess.run(cmd, check=check)


def command_output(cmd):
    return subprocess.check_output(cmd, text=True).strip()


def find_mafft_prefix():
    """
    Locate Homebrew MAFFT prefix.

    Usually:
        /opt/homebrew/opt/mafft
    or resolved:
        /opt/homebrew/Cellar/mafft/<version>
    """
    # Best: brew --prefix mafft
    try:
        prefix = command_output(["brew", "--prefix", "mafft"])
        p = Path(prefix).resolve()
        if (p / "bin" / "mafft").exists() or (p / "libexec" / "bin" / "mafft").exists():
            return p
    except Exception:
        pass

    # Fallback: mafft executable path.
    mafft = shutil.which("mafft")

    if not mafft:
        raise RuntimeError(
            "MAFFT was not found.\n\n"
            "Install it first:\n"
            "  brew install mafft"
        )

    real = Path(mafft).resolve()

    for parent in real.parents:
        if (parent / "bin" / "mafft").exists() or (parent / "libexec" / "bin" / "mafft").exists():
            if parent.name.lower() != "bin":
                return parent

    raise RuntimeError(f"Could not determine MAFFT prefix from: {real}")


def ensure_pyinstaller():
    try:
        import PyInstaller  # noqa: F401
        return
    except Exception:
        run([sys.executable, "-m", "pip", "install", "pyinstaller"])


def build():
    project_dir = Path(__file__).resolve().parent
    os.chdir(project_dir)

    mafft_prefix = find_mafft_prefix()
    print(f"\nMAFFT prefix detected:\n{mafft_prefix}")

    ensure_pyinstaller()

    # Clean previous build.
    for path in ["build", "dist"]:
        p = project_dir / path
        if p.exists():
            shutil.rmtree(p)

    spec_file = project_dir / f"{APP_NAME}.spec"
    if spec_file.exists():
        spec_file.unlink()

    add_data_arg = f"{mafft_prefix}:mafft_bundle"

    cmd = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--name",
        APP_NAME,
        "--windowed",
        "--onedir",
        "--clean",
        "--collect-all",
        "PyQt6",
        "--collect-all",
        "Bio",
        "--collect-all",
        "reportlab",
        "--add-data",
        add_data_arg,
        "main.py",
    ]

    run(cmd)

    app_path = project_dir / "dist" / f"{APP_NAME}.app"

    if not app_path.exists():
        raise RuntimeError(f"Build finished but app was not found: {app_path}")

    # Remove quarantine on the local machine.
    run(["xattr", "-cr", str(app_path)], check=False)

    print("\nBuild complete:")
    print(app_path)
    print("\nTest it with:")
    print(f'  open "{app_path}"')
    print("\nTo zip for sharing:")
    print("  cd dist")
    print(f'  zip -r "{APP_NAME} macOS bundled MAFFT.zip" "{APP_NAME}.app"')


if __name__ == "__main__":
    build()
