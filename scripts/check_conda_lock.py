import filecmp
import subprocess
import tempfile
from pathlib import Path


def check_conda_lock():
    # 現在のconda-lock.ymlの一時コピーを作成
    current_lock = Path("conda-lock.yml")
    if not current_lock.exists():
        print("No existing conda-lock.yml found")
        return 1

    with tempfile.NamedTemporaryFile(suffix=".yml") as temp_file:
        # 新しいlockファイルを一時ファイルに生成
        try:
            subprocess.run(
                [
                    "conda-lock",
                    "-f",
                    "environment.yml",
                    "--platform",
                    "linux-64",
                    "-f",
                    temp_file.name,
                ],
                capture_output=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"Error generating new lock file: {e.stderr}")
            return 1

        # ファイルを比較
        if not filecmp.cmp(current_lock, temp_file.name):
            print(
                "conda-lock.yml is out of date. Please run 'conda-lock -f environment.yml' to update it."
            )
            return 1

    print("conda-lock.yml is up to date")
    return 0
