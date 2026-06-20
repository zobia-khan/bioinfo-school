# AGENTS.md — bioinfo-school

## Environment & Package Management
- **Python ≥ 3.8** is required; target 3.11+ for local dev.
- Use **`uv`** as the package and environment manager (`uv venv`, `uv pip install`).
- The main QC project lives in `exercises/week2/project/` and has its own
  [`pyproject.toml`](exercises/week2/project/pyproject.toml). Install it editably with
  `uv pip install -e exercises/week2/project/`.

## Data & File Conventions
- Raw sequencing data (`.fastq`, `.fastq.gz`) lives alongside the script that
  processes it (e.g. `exercises/week2/project/test.fastq`). Do **not** commit
  large FASTQ or BAM files; they are listed in `.gitignore`.
- VCF files follow standard VCF 4.x format. Chromosome names use UCSC style
  (`chr1`, `chrX`) unless otherwise noted in the exercise instructions.
- **Genome build: GRCh38 / hg38** for all human data.

## Quality Scores & Bioinformatics Conventions
- Quality scores are **Phred-33 encoded** (Illumina 1.8+). Assume this unless
  a file explicitly states otherwise.
- A read/base is considered **low-quality** at Q < 20 (FAIL) or Q < 25 (WARN),
  matching the thresholds in `fastq_qc.py`.

## Validation Checks Before Claiming Success
1. Run the QC tool on the included mock data and confirm an HTML report is produced:
   `python exercises/week2/project/fastq_qc.py exercises/week2/project/test.fastq`
2. Confirm the script exits 0 and the output HTML opens without JavaScript errors.
3. For VCF work, verify `parse_vcf.py` prints a per-chromosome summary table
   without Python tracebacks.
4. Check that no new files > 5 MB are staged (`git diff --stat --cached`).

## Files the Agent Must Not Edit or Commit
- `.venv/` — virtual environment directory, never commit.
- `*.fastq`, `*.fastq.gz`, `*.bam`, `*.bai` — raw sequencing data.
- `my_report.html` and `test_qc_report.html` — generated artefacts, not source.
- `LICENSE` — do not modify.
