from __future__ import annotations

from pathlib import Path
from typing import Optional, Iterable
import subprocess
import shutil
import pandas as pd


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_dataframe(df: pd.DataFrame, path: Path, sep: str = "\t") -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep=sep, index=False)


def load_dataframe(path: Path, sep: str = "\t") -> pd.DataFrame:
    return pd.read_csv(path, sep=sep)


def which(executable: str) -> Optional[str]:
    return shutil.which(executable)


def run_command(
    cmd: list[str],
    cwd: Optional[Path] = None,
    check: bool = True,
    capture_output: bool = False,
) -> subprocess.CompletedProcess:
    """
    Thin subprocess wrapper. Raises a readable error if command fails.
    """
    try:
        return subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            check=check,
            text=True,
            capture_output=capture_output,
        )
    except FileNotFoundError as e:
        raise RuntimeError(f"Executable not found: {cmd[0]}") from e
    except subprocess.CalledProcessError as e:
        msg = f"Command failed: {' '.join(cmd)}"
        if e.stdout:
            msg += f"\nSTDOUT:\n{e.stdout}"
        if e.stderr:
            msg += f"\nSTDERR:\n{e.stderr}"
        raise RuntimeError(msg) from e
