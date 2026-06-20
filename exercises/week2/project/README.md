# FASTQ Quality Control Tool

A standalone Python tool that parses a FASTQ file (plain or gzip-compressed) and produces a self-contained, interactive HTML quality-control report — similar to FastQC but requiring no external bioinformatics dependencies.

---

## Required inputs

| Argument | Description |
|---|---|
| `<path_to_fastq>` | Path to a `.fastq` or `.fastq.gz` file (Phred-33 encoded, standard 4-line FASTQ format) |
| `[output_report.html]` | *(optional)* Path for the output HTML report. Defaults to `<stem>_qc_report.html` next to the input file. |

---

## Outputs

| File | Description |
|---|---|
| `<stem>_qc_report.html` | Self-contained HTML report with interactive charts (no server needed — open in any browser) |

The report includes:

- **Overview** — total reads, total bases, average read length, mean GC %, QC status table (PASS / WARN / FAIL), overrepresented sequences
- **Per-Base Quality** — Phred score distribution (median, IQR, P10/P90) at every cycle position
- **Per-Sequence Quality** — histogram of mean quality scores across all reads
- **GC Content Distribution** — per-read GC % histogram
- **N-Content & Length** — ambiguous base (N) rate per position and read-length histogram

---

## How to rerun

### 1. Run QC on an existing FASTQ file

```bash
# Minimal (output defaults to test_qc_report.html)
python fastq_qc.py test.fastq

# Explicit output path
python fastq_qc.py test.fastq my_report.html

# Gzip input
python fastq_qc.py reads.fastq.gz reads_qc_report.html
```

### 2. Generate mock test data, then run QC

```bash
# Creates test.fastq with 1 000 synthetic reads (seed = 42, reproducible)
python generate_mock_fastq.py

# Then run QC on the generated file
python fastq_qc.py test.fastq
```

---

## Dependencies

### Runtime (stdlib only — no `pip install` required)

| Module | Purpose |
|---|---|
| `sys`, `os` | argument parsing, path manipulation |
| `gzip` | transparent decompression of `.gz` inputs |
| `json` | serialising QC data for embedding in the HTML report |
| `pathlib.Path` | cross-platform file path handling |
| `collections.Counter`, `defaultdict` | per-position quality / base frequency tracking |
| `random` *(generate_mock_fastq.py only)* | reproducible synthetic read generation |

> **All dependencies ship with Python ≥ 3.8 — no third-party packages are needed.**

### Browser-side (CDN, no install)

| Library | Purpose |
|---|---|
| [Chart.js](https://www.chartjs.org/) (CDN) | interactive charts in the HTML report |
| [Google Fonts — Outfit + Plus Jakarta Sans](https://fonts.google.com/) (CDN) | typography |

> The HTML report requires an internet connection the **first time** it is opened (to load Chart.js and fonts from CDN). Subsequent opens work offline once the browser has cached them.

---

## Fresh-environment setup with `uv`

```bash
# 1. Create and activate a virtual environment
uv venv .venv
# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate

# 2. Install the package (editable, no extra deps needed)
uv pip install -e .

# 3. Use the installed console script
fastq-qc test.fastq
```

---

## QC thresholds

| Module | WARN | FAIL |
|---|---|---|
| Per-base quality (lowest median) | < Q25 | < Q20 |
| Per-sequence quality (mean) | < Q27 | < Q20 |
| N-content (max per position) | > 2 % | > 5 % |
| Overrepresented sequences | ≥ 0.5 % of reads (reported, no PASS/WARN/FAIL) | — |
