import subprocess
import pytest

DATA_DIR = "tests/_data"
FASTQ_FILE = DATA_DIR + "/example_fastq.fastq"


@pytest.mark.parametrize(
    "args, should_succeed",
    [
        (["-i", str(FASTQ_FILE)], True),
        (["-i", str(FASTQ_FILE), "--gc", "30", "60"], True),
        (["-i", str(FASTQ_FILE), "--gc", "70"], True),
        (["-i", str(FASTQ_FILE), "-l", "100", "500"], True),
        (["-i", str(FASTQ_FILE), "-l", "300"], True),
        (["-i", str(FASTQ_FILE), "-q", "30"], True),
        (["--input", str(FASTQ_FILE), "--output", "filtered.fastq"], True),
        (["-i", str(FASTQ_FILE), "--gc", "abc"], False),
        (["-i", str(FASTQ_FILE), "-l", "abc"], False),
        (["-i", str(FASTQ_FILE), "--gc"], False),
        (["-i", str(FASTQ_FILE), "-l"], False),
        (["-i", str(FASTQ_FILE), "-q"], False),
    ],
)
def test_cli_run(args, should_succeed):
    result = subprocess.run(
        ["python", "filter_fastq_CLI.py"] + args, capture_output=True, text=True
    )

    if should_succeed:
        assert (
            result.returncode == 0
        ), f"Failed with stdout: {result.stdout}, stderr: {result.stderr}"
    else:
        assert result.returncode != 0, "Should have failed but didn't"
