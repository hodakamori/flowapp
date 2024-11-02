import subprocess
import sys


def update_conda_lock():
    # conda-lockファイルの更新
    try:
        result = subprocess.run(
            [
                "conda-lock",
                "-f",
                "environment.yml",
                "--platform",
                "linux-64",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        print("Successfully updated conda-lock files")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Error updating conda-lock files: {e.stderr}")
        return 1


if __name__ == "__main__":
    sys.exit(update_conda_lock())
