#!/usr/bin/env python3
import sys
import os
import gzip
import json
from pathlib import Path
from collections import Counter, defaultdict

def get_fastq_iterator(filepath):
    """
    Returns an iterator over the lines of a FASTQ file.
    Supports gzipped files automatically if name ends in .gz.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"FASTQ file '{filepath}' not found.")
        
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt")
    else:
        return open(filepath, "r")

def compute_percentiles(freq_dict, percentiles=[10, 25, 50, 75, 90]):
    """
    Computes percentiles from a frequency dict of integer quality scores.
    """
    total = sum(freq_dict.values())
    if total == 0:
        return [0] * len(percentiles)
        
    sorted_scores = sorted(freq_dict.keys())
    results = {}
    
    current_cumulative = 0
    for score in sorted_scores:
        current_cumulative += freq_dict[score]
        for p in percentiles:
            if p not in results and (current_cumulative / total) >= (p / 100.0):
                results[p] = score
                
    # Fill in any missing (in case of rounding edge cases)
    for p in percentiles:
        if p not in results:
            results[p] = sorted_scores[-1] if sorted_scores else 0
            
    return [results[p] for p in percentiles]

def run_qc(filepath):
    """
    Parses a FASTQ file and accumulates QC metrics.
    """
    total_reads = 0
    total_bases = 0
    
    # Distributions
    gc_counts = Counter()  # GC content per read rounded to nearest integer percentage
    read_lengths = Counter()
    mean_quals = Counter()  # Mean Phred score per read rounded to nearest integer
    
    # Per-position tracking
    # Dict of position -> Counter of Phred score occurrences
    pos_qual_freqs = defaultdict(Counter)
    pos_n_count = defaultdict(int)
    pos_base_count = defaultdict(int)
    
    # Overrepresented sequences (we track the first 30bp of reads)
    seq_prefixes = Counter()
    
    # Read FASTQ in 4-line blocks
    try:
        f = get_fastq_iterator(filepath)
    except Exception as e:
        print(f"Error opening file: {e}")
        sys.exit(1)
        
    line_idx = 0
    seq = None
    header = None
    
    for line in f:
        line = line.strip()
        mod = line_idx % 4
        
        if mod == 0:
            header = line
        elif mod == 1:
            seq = line
        elif mod == 2:
            pass  # + separator
        elif mod == 3:
            # Quality scores line
            qual = line
            
            # Sequence validation
            if len(seq) != len(qual):
                print(f"Warning: Read length ({len(seq)}) does not match quality length ({len(qual)}) at read {total_reads + 1}.")
                line_idx += 1
                continue
                
            total_reads += 1
            read_len = len(seq)
            total_bases += read_len
            read_lengths[read_len] += 1
            
            # GC Content calculation
            g_count = seq.count("G") + seq.count("g")
            c_count = seq.count("C") + seq.count("c")
            gc_pct = int(round((g_count + c_count) / read_len * 100)) if read_len > 0 else 0
            gc_counts[gc_pct] += 1
            
            # Track prefix for overrepresented sequences (limit to first 30bp)
            if read_len >= 30:
                seq_prefixes[seq[:30]] += 1
            else:
                seq_prefixes[seq] += 1
                
            # Quality score calculations (Phred-33)
            read_qual_sum = 0
            for pos, (char, base) in enumerate(zip(qual, seq)):
                q_score = ord(char) - 33
                read_qual_sum += q_score
                
                pos_qual_freqs[pos][q_score] += 1
                pos_base_count[pos] += 1
                if base.upper() == "N":
                    pos_n_count[pos] += 1
                    
            mean_q = int(round(read_qual_sum / read_len)) if read_len > 0 else 0
            mean_quals[mean_q] += 1
            
        line_idx += 1
        
    f.close()
    
    if total_reads == 0:
        print("Error: No valid reads found in FASTQ file.")
        sys.exit(1)
        
    # Process per-position statistics
    max_position = max(pos_qual_freqs.keys()) if pos_qual_freqs else 0
    positions = list(range(max_position + 1))
    
    per_base_stats = []
    for pos in positions:
        freqs = pos_qual_freqs[pos]
        total_at_pos = pos_base_count[pos]
        
        # Calculate mean quality score
        qual_sum = sum(score * count for score, count in freqs.items())
        mean_qual = qual_sum / total_at_pos if total_at_pos > 0 else 0
        
        # Compute percentiles
        p10, p25, p50, p75, p90 = compute_percentiles(freqs)
        
        # N percent
        n_pct = (pos_n_count[pos] / total_at_pos * 100) if total_at_pos > 0 else 0.0
        
        per_base_stats.append({
            "position": pos + 1,
            "mean": round(mean_qual, 2),
            "p10": p10,
            "p25": p25,
            "median": p50,
            "p75": p75,
            "p90": p90,
            "n_pct": round(n_pct, 4)
        })
        
    # GC Stats
    overall_gc = sum(pct * count for pct, count in gc_counts.items()) / total_reads if total_reads > 0 else 0.0
    
    # Overrepresented sequences (keep top 5 that are > 0.5% of total reads)
    overrepresented = []
    for sequence, count in seq_prefixes.most_common(5):
        pct = (count / total_reads * 100)
        if pct >= 0.5:
            overrepresented.append({
                "sequence": sequence,
                "count": count,
                "percentage": round(pct, 2)
            })
            
    # Sequence length summary
    lengths = sorted(read_lengths.keys())
    min_len = lengths[0] if lengths else 0
    max_len = lengths[-1] if lengths else 0
    avg_len = total_bases / total_reads if total_reads > 0 else 0
    
    # Determine Statuses (Pass/Warn/Fail)
    # 1. Per-base quality status: Fail if median quality at any position < 20, Warn if < 25
    lowest_median = min(stat["median"] for stat in per_base_stats) if per_base_stats else 0
    if lowest_median < 20:
        base_qual_status = "FAIL"
    elif lowest_median < 25:
        base_qual_status = "WARN"
    else:
        base_qual_status = "PASS"
        
    # 2. Per-sequence quality status: Fail if mean read quality < 20, Warn if < 27
    avg_read_qual = sum(score * count for score, count in mean_quals.items()) / total_reads if total_reads > 0 else 0
    if avg_read_qual < 20:
        seq_qual_status = "FAIL"
    elif avg_read_qual < 27:
        seq_qual_status = "WARN"
    else:
        seq_qual_status = "PASS"
        
    # 3. N-Content status: Fail if max N-content > 5%, Warn if > 2%
    max_n_pct = max(stat["n_pct"] for stat in per_base_stats) if per_base_stats else 0.0
    if max_n_pct > 5.0:
        n_content_status = "FAIL"
    elif max_n_pct > 2.0:
        n_content_status = "WARN"
    else:
        n_content_status = "PASS"
        
    qc_data = {
        "filename": os.path.basename(filepath),
        "total_reads": total_reads,
        "total_bases": total_bases,
        "min_len": min_len,
        "max_len": max_len,
        "avg_len": round(avg_len, 2),
        "overall_gc": round(overall_gc, 2),
        "avg_read_qual": round(avg_read_qual, 2),
        "lowest_median": lowest_median,
        "max_n_pct": round(max_n_pct, 4),
        "statuses": {
            "base_quality": base_qual_status,
            "seq_quality": seq_qual_status,
            "n_content": n_content_status
        },
        "distributions": {
            "gc": [{"gc": pct, "count": count} for pct, count in sorted(gc_counts.items())],
            "read_lengths": [{"length": l, "count": count} for l, count in sorted(read_lengths.items())],
            "mean_quals": [{"qual": q, "count": count} for q, count in sorted(mean_quals.items())]
        },
        "per_base_stats": per_base_stats,
        "overrepresented": overrepresented
    }
    
    return qc_data

def generate_html_report(qc_data, output_path):
    """
    Generates a premium, beautiful HTML report containing the QC results.
    """
    json_data_str = json.dumps(qc_data, indent=2)
    
    # Render overrepresented sequences table rows
    overrepresented_rows = ""
    if qc_data['overrepresented']:
        for item in qc_data['overrepresented']:
            overrepresented_rows += f"""
            <tr>
                <td><span class="overrepresented-seq">{item['sequence']}</span></td>
                <td>{item['count']:,}</td>
                <td>{item['percentage']}%</td>
            </tr>"""
    else:
        overrepresented_rows = """
        <tr>
            <td colspan="3" style="text-align: center; color: #94a3b8;">No overrepresented sequences found (above 0.5% threshold)</td>
        </tr>"""
        
    # Standard HTML template using simple placeholders instead of python f-string formatting
    # to avoid single curly braces clash with Javascript/CSS
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FASTQ Quality Control Report - __FILENAME__</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=Plus+Jakarta+Sans:wght@300;400;500;600;700&display=swap" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --primary: #8b5cf6;
            --primary-light: #a78bfa;
            --success: #10b981;
            --success-bg: rgba(16, 185, 129, 0.1);
            --warning: #f59e0b;
            --warning-bg: rgba(245, 158, 11, 0.1);
            --danger: #ef4444;
            --danger-bg: rgba(239, 68, 68, 0.1);
            --accent-teal: #14b8a6;
            --accent-pink: #ec4899;
        }

        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }

        body {
            font-family: 'Plus Jakarta Sans', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem;
            min-height: 100vh;
        }

        .container {
            max-width: 1300px;
            margin: 0 auto;
        }

        /* Header block */
        header {
            background: linear-gradient(135deg, #1e1b4b 0%, #1e293b 100%);
            border: 1px solid var(--card-border);
            border-radius: 1.25rem;
            padding: 2.5rem;
            margin-bottom: 2rem;
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.25);
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 1.5rem;
        }

        .header-title h1 {
            font-family: 'Outfit', sans-serif;
            font-size: 2.25rem;
            font-weight: 700;
            background: linear-gradient(to right, #a78bfa, #f472b6);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.5rem;
        }

        .header-title p {
            color: var(--text-secondary);
            font-size: 1rem;
            font-weight: 500;
        }

        .header-badge {
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            padding: 0.75rem 1.25rem;
            border-radius: 50px;
            font-weight: 600;
            font-size: 0.9rem;
        }

        .status-dot {
            width: 10px;
            height: 10px;
            border-radius: 50%;
            display: inline-block;
        }

        .status-PASS { background-color: var(--success); box-shadow: 0 0 10px var(--success); }
        .status-WARN { background-color: var(--warning); box-shadow: 0 0 10px var(--warning); }
        .status-FAIL { background-color: var(--danger); box-shadow: 0 0 10px var(--danger); }

        /* Tabs and Navigation */
        .tabs {
            display: flex;
            gap: 0.5rem;
            margin-bottom: 2rem;
            background: rgba(30, 41, 59, 0.5);
            padding: 0.35rem;
            border-radius: 0.75rem;
            border: 1px solid var(--card-border);
            overflow-x: auto;
        }

        .tab-btn {
            background: transparent;
            border: none;
            color: var(--text-secondary);
            padding: 0.75rem 1.5rem;
            font-size: 0.95rem;
            font-weight: 600;
            border-radius: 0.5rem;
            cursor: pointer;
            transition: all 0.2s ease;
            white-space: nowrap;
        }

        .tab-btn:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.05);
        }

        .tab-btn.active {
            background: var(--primary);
            color: var(--text-primary);
            box-shadow: 0 4px 12px rgba(139, 92, 246, 0.3);
        }

        .tab-content {
            display: none;
        }

        .tab-content.active {
            display: block;
            animation: fadeIn 0.4s ease;
        }

        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }

        /* Dashboard Overview Grid */
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }

        .stat-card {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 1rem;
            padding: 1.75rem;
            position: relative;
            overflow: hidden;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }

        .stat-card:hover {
            transform: translateY(-2px);
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.3), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
        }

        .stat-card::after {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 4px;
            background: linear-gradient(to right, var(--primary), var(--accent-teal));
        }

        .stat-title {
            color: var(--text-secondary);
            font-size: 0.85rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-bottom: 0.5rem;
        }

        .stat-value {
            font-family: 'Outfit', sans-serif;
            font-size: 2rem;
            font-weight: 700;
            color: var(--text-primary);
        }

        /* Main panels and Layouts */
        .panel {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 1rem;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        }

        .panel-header {
            margin-bottom: 1.5rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }

        .panel-title h2 {
            font-family: 'Outfit', sans-serif;
            font-size: 1.35rem;
            font-weight: 600;
            color: var(--text-primary);
        }

        .panel-title p {
            color: var(--text-secondary);
            font-size: 0.85rem;
            margin-top: 0.25rem;
        }

        /* Chart container */
        .chart-container {
            position: relative;
            width: 100%;
            height: 450px;
        }

        /* Tables */
        .table-wrapper {
            width: 100%;
            overflow-x: auto;
            border-radius: 0.75rem;
            border: 1px solid var(--card-border);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            text-align: left;
            font-size: 0.9rem;
        }

        th {
            background-color: rgba(15, 23, 42, 0.6);
            color: var(--text-primary);
            padding: 1rem 1.25rem;
            font-weight: 600;
            border-bottom: 1px solid var(--card-border);
        }

        td {
            padding: 1rem 1.25rem;
            color: var(--text-secondary);
            border-bottom: 1px solid var(--card-border);
        }

        tr:last-child td {
            border-bottom: none;
        }

        tr:hover td {
            background-color: rgba(255, 255, 255, 0.02);
            color: var(--text-primary);
        }

        .status-pill {
            display: inline-flex;
            align-items: center;
            gap: 0.35rem;
            padding: 0.35rem 0.75rem;
            border-radius: 50px;
            font-size: 0.75rem;
            font-weight: 700;
            letter-spacing: 0.05em;
        }

        .status-pill-PASS { background-color: var(--success-bg); color: var(--success); }
        .status-pill-WARN { background-color: var(--warning-bg); color: var(--warning); }
        .status-pill-FAIL { background-color: var(--danger-bg); color: var(--danger); }

        /* Custom UI cards for specific values */
        .overrepresented-seq {
            font-family: monospace;
            background: #0f172a;
            padding: 0.4rem 0.6rem;
            border-radius: 0.25rem;
            color: var(--accent-pink);
            font-size: 0.85rem;
            word-break: break-all;
        }

        /* Responsive styling */
        @media(max-width: 768px) {
            body { padding: 1rem; }
            header { padding: 1.5rem; }
            .chart-container { height: 350px; }
        }
    </style>
</head>
<body>
    <div class="container">
        <!-- Dashboard Header -->
        <header>
            <div class="header-title">
                <h1>FASTQ Quality Control</h1>
                <p>File: __FILENAME__</p>
            </div>
            <div>
                <div class="header-badge">
                    <span>Overall Quality:</span>
                    <span class="status-dot status-__OVERALL_STATUS__"></span>
                    <span>__OVERALL_STATUS__</span>
                </div>
            </div>
        </header>

        <!-- Navigation Tabs -->
        <div class="tabs">
            <button class="tab-btn active" onclick="switchTab('overview')">Overview</button>
            <button class="tab-btn" onclick="switchTab('base-quality')">Per-Base Quality</button>
            <button class="tab-btn" onclick="switchTab('seq-quality')">Per-Seq Quality</button>
            <button class="tab-btn" onclick="switchTab('gc-content')">GC Distribution</button>
            <button class="tab-btn" onclick="switchTab('n-content')">N-Content & Length</button>
        </div>

        <!-- OVERVIEW TAB -->
        <div id="overview" class="tab-content active">
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-title">Total Reads</div>
                    <div class="stat-value">__TOTAL_READS__</div>
                </div>
                <div class="stat-card">
                    <div class="stat-title">Total Bases</div>
                    <div class="stat-value">__TOTAL_BASES__</div>
                </div>
                <div class="stat-card">
                    <div class="stat-title">Avg Read Length</div>
                    <div class="stat-value">__AVG_LEN__ bp</div>
                </div>
                <div class="stat-card">
                    <div class="stat-title">Overall GC</div>
                    <div class="stat-value">__OVERALL_GC__%</div>
                </div>
            </div>

            <div class="panel">
                <div class="panel-header">
                    <div class="panel-title">
                        <h2>QC Status Summary</h2>
                        <p>Module outcomes based on standard thresholds</p>
                    </div>
                </div>
                <div class="table-wrapper">
                    <table>
                        <thead>
                            <tr>
                                <th>QC Assessment Module</th>
                                <th>Status</th>
                                <th>Details</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>Per-Base Quality Scores</td>
                                <td>
                                    <span class="status-pill status-pill-__BASE_QUAL_STATUS__">
                                        __BASE_QUAL_STATUS__
                                    </span>
                                </td>
                                <td>Lowest median Phred score along read is __LOWEST_MEDIAN__</td>
                            </tr>
                            <tr>
                                <td>Per-Sequence Quality Scores</td>
                                <td>
                                    <span class="status-pill status-pill-__SEQ_QUAL_STATUS__">
                                        __SEQ_QUAL_STATUS__
                                    </span>
                                </td>
                                <td>Average quality score per read is __AVG_READ_QUAL__</td>
                            </tr>
                            <tr>
                                <td>N-Content per Base Position</td>
                                <td>
                                    <span class="status-pill status-pill-__N_CONTENT_STATUS__">
                                        __N_CONTENT_STATUS__
                                    </span>
                                </td>
                                <td>Maximum N-content at any cycle is __MAX_N_PCT__%</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Overrepresented Sequences panel -->
            <div class="panel">
                <div class="panel-header">
                    <div class="panel-title">
                        <h2>Overrepresented Sequences</h2>
                        <p>Sequences (first 30bp) matching more than 0.5% of total reads (often adapters or low diversity)</p>
                    </div>
                </div>
                <div class="table-wrapper">
                    <table>
                        <thead>
                            <tr>
                                <th>Sequence Prefix (First 30bp)</th>
                                <th>Read Count</th>
                                <th>Percentage</th>
                            </tr>
                        </thead>
                        <tbody>
                            __OVERREPRESENTED_ROWS__
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <!-- PER-BASE QUALITY TAB -->
        <div id="base-quality" class="tab-content">
            <div class="panel">
                <div class="panel-header">
                    <div class="panel-title">
                        <h2>Per-Base Quality Scores</h2>
                        <p>Distribution of Phred quality scores ($Q$) at each position/cycle. Higher is better (Q &gt; 30 represents 99.9% accuracy).</p>
                    </div>
                </div>
                <div class="chart-container">
                    <canvas id="baseQualChart"></canvas>
                </div>
            </div>
        </div>

        <!-- PER-SEQUENCE QUALITY TAB -->
        <div id="seq-quality" class="tab-content">
            <div class="panel">
                <div class="panel-header">
                    <div class="panel-title">
                        <h2>Per-Sequence Mean Quality Scores</h2>
                        <p>Frequency of reads with specific average quality scores.</p>
                    </div>
                </div>
                <div class="chart-container">
                    <canvas id="seqQualChart"></canvas>
                </div>
            </div>
        </div>

        <!-- GC DISTRIBUTION TAB -->
        <div id="gc-content" class="tab-content">
            <div class="panel">
                <div class="panel-header">
                    <div class="panel-title">
                        <h2>GC Content Distribution</h2>
                        <p>Percentage GC content of reads compared to a normal distribution.</p>
                    </div>
                </div>
                <div class="chart-container">
                    <canvas id="gcChart"></canvas>
                </div>
            </div>
        </div>

        <!-- N-CONTENT & LENGTH TAB -->
        <div id="n-content" class="tab-content">
            <div class="panel">
                <div class="panel-header">
                    <div class="panel-title">
                        <h2>Ambiguous Base (N) Content per Position</h2>
                        <p>Percentage of N bases called at each cycle position.</p>
                    </div>
                </div>
                <div class="chart-container" style="height: 300px; margin-bottom: 2rem;">
                    <canvas id="nContentChart"></canvas>
                </div>
            </div>

            <div class="panel">
                <div class="panel-header">
                    <div class="panel-title">
                        <h2>Sequence Length Distribution</h2>
                        <p>Length of reads in the library.</p>
                    </div>
                </div>
                <div class="chart-container" style="height: 300px;">
                    <canvas id="lengthChart"></canvas>
                </div>
            </div>
        </div>
    </div>

    <!-- Inject data for charts -->
    <script>
        const qcData = __JSON_DATA__;

        function switchTab(tabId) {
            document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(content => content.classList.remove('active'));
            
            const eventBtn = Array.from(document.querySelectorAll('.tab-btn')).find(btn => btn.getAttribute('onclick').includes(tabId));
            if (eventBtn) eventBtn.classList.add('active');
            
            const activeContent = document.getElementById(tabId);
            if (activeContent) activeContent.classList.add('active');
        }

        document.addEventListener('DOMContentLoaded', () => {
            // Configure global Chart defaults for premium aesthetics
            Chart.defaults.color = '#94a3b8';
            Chart.defaults.font.family = "'Plus Jakarta Sans', sans-serif";
            Chart.defaults.font.size = 12;
            Chart.defaults.borderColor = 'rgba(51, 65, 85, 0.4)';
            
            const xLabels = qcData.per_base_stats.map(s => s.position);
            
            // 1. Per-Base Quality Box Chart (Line-based simulation for simplicity and max compatibility)
            const ctxQual = document.getElementById('baseQualChart').getContext('2d');
            new Chart(ctxQual, {
                type: 'line',
                data: {
                    labels: xLabels,
                    datasets: [
                        {
                            label: 'Median Quality',
                            data: qcData.per_base_stats.map(s => s.median),
                            borderColor: '#8b5cf6',
                            backgroundColor: 'rgba(139, 92, 246, 0.2)',
                            borderWidth: 3,
                            fill: false,
                            tension: 0.1,
                            yAxisID: 'y'
                        },
                        {
                            label: 'Mean Quality',
                            data: qcData.per_base_stats.map(s => s.mean),
                            borderColor: '#14b8a6',
                            borderWidth: 2,
                            borderDash: [5, 5],
                            fill: false,
                            tension: 0.1
                        },
                        {
                            label: '25th - 75th Percentile Range',
                            data: qcData.per_base_stats.map(s => s.p75),
                            borderColor: 'transparent',
                            backgroundColor: 'rgba(139, 92, 246, 0.08)',
                            fill: '+1',
                            tension: 0.1
                        },
                        {
                            label: '25th Percentile boundary',
                            data: qcData.per_base_stats.map(s => s.p25),
                            borderColor: 'transparent',
                            backgroundColor: 'transparent',
                            fill: false,
                            tension: 0.1
                        }
                    ]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        legend: {
                            labels: {
                                filter: function(item, chart) {
                                    // Filter out helper boundaries for cleaner legend
                                    return !item.text.includes('boundary');
                                }
                            }
                        }
                    },
                    scales: {
                        y: {
                            min: 0,
                            max: 42,
                            title: {
                                display: true,
                                text: 'Quality Score (Phred)'
                            },
                            grid: {
                                color: function(context) {
                                    if (context.tick.value < 20) return 'rgba(239, 68, 68, 0.15)';
                                    if (context.tick.value < 28) return 'rgba(245, 158, 11, 0.15)';
                                    return 'rgba(16, 185, 129, 0.15)';
                                }
                            }
                        },
                        x: {
                            title: {
                                display: true,
                                text: 'Position in Read (bp)'
                            }
                        }
                    }
                }
            });

            // 2. Per-Sequence Quality Chart
            const ctxSeqQual = document.getElementById('seqQualChart').getContext('2d');
            const meanQualData = qcData.distributions.mean_quals;
            new Chart(ctxSeqQual, {
                type: 'bar',
                data: {
                    labels: meanQualData.map(d => d.qual),
                    datasets: [{
                        label: 'Read Count',
                        data: meanQualData.map(d => d.count),
                        backgroundColor: 'rgba(20, 184, 166, 0.6)',
                        borderColor: '#14b8a6',
                        borderWidth: 1.5,
                        borderRadius: 4
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            beginAtZero: true,
                            title: { display: true, text: 'Number of Reads' }
                        },
                        x: {
                            title: { display: true, text: 'Mean Quality Score' }
                        }
                    }
                }
            });

            // 3. GC Content Distribution
            const ctxGC = document.getElementById('gcChart').getContext('2d');
            const gcData = qcData.distributions.gc;
            new Chart(ctxGC, {
                type: 'line',
                data: {
                    labels: gcData.map(d => d.gc),
                    datasets: [{
                        label: 'GC Content (%)',
                        data: gcData.map(d => d.count),
                        borderColor: '#ec4899',
                        backgroundColor: 'rgba(236, 72, 153, 0.1)',
                        fill: true,
                        tension: 0.3,
                        borderWidth: 2.5
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            beginAtZero: true,
                            title: { display: true, text: 'Number of Reads' }
                        },
                        x: {
                            title: { display: true, text: 'GC % per Read' }
                        }
                    }
                }
            });

            // 4. N Content Chart
            const ctxN = document.getElementById('nContentChart').getContext('2d');
            new Chart(ctxN, {
                type: 'line',
                data: {
                    labels: xLabels,
                    datasets: [{
                        label: 'Percentage N',
                        data: qcData.per_base_stats.map(s => s.n_pct),
                        borderColor: '#ef4444',
                        backgroundColor: 'rgba(239, 68, 68, 0.05)',
                        fill: true,
                        tension: 0.1,
                        borderWidth: 2
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            min: 0,
                            title: { display: true, text: '% N Called' }
                        },
                        x: {
                            title: { display: true, text: 'Position (bp)' }
                        }
                    }
                }
            });

            // 5. Sequence Length Distribution
            const ctxLen = document.getElementById('lengthChart').getContext('2d');
            const lenData = qcData.distributions.read_lengths;
            new Chart(ctxLen, {
                type: 'bar',
                data: {
                    labels: lenData.map(d => d.length),
                    datasets: [{
                        label: 'Read Count',
                        data: lenData.map(d => d.count),
                        backgroundColor: 'rgba(245, 158, 11, 0.6)',
                        borderColor: '#f59e0b',
                        borderWidth: 1.5,
                        borderRadius: 4
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            beginAtZero: true,
                            title: { display: true, text: 'Number of Reads' }
                        },
                        x: {
                            title: { display: true, text: 'Read Length (bp)' }
                        }
                    }
                }
            });
        });
    </script>
</body>
</html>
"""
    
    # Perform string replacements to construct final page
    html_content = html_template
    html_content = html_content.replace("__FILENAME__", str(qc_data['filename']))
    html_content = html_content.replace("__OVERALL_STATUS__", str(qc_data['statuses']['base_quality']))
    html_content = html_content.replace("__TOTAL_READS__", f"{qc_data['total_reads']:,}")
    html_content = html_content.replace("__TOTAL_BASES__", f"{qc_data['total_bases']:,}")
    html_content = html_content.replace("__AVG_LEN__", str(qc_data['avg_len']))
    html_content = html_content.replace("__OVERALL_GC__", str(qc_data['overall_gc']))
    html_content = html_content.replace("__BASE_QUAL_STATUS__", str(qc_data['statuses']['base_quality']))
    html_content = html_content.replace("__SEQ_QUAL_STATUS__", str(qc_data['statuses']['seq_quality']))
    html_content = html_content.replace("__N_CONTENT_STATUS__", str(qc_data['statuses']['n_content']))
    html_content = html_content.replace("__LOWEST_MEDIAN__", str(qc_data['lowest_median']))
    html_content = html_content.replace("__AVG_READ_QUAL__", str(qc_data['avg_read_qual']))
    html_content = html_content.replace("__MAX_N_PCT__", f"{qc_data['max_n_pct']:.2f}")
    html_content = html_content.replace("__OVERREPRESENTED_ROWS__", overrepresented_rows)
    html_content = html_content.replace("__JSON_DATA__", json_data_str)
    
    with open(output_path, "w") as f:
        f.write(html_content)
        
    print(f"HTML report successfully generated at '{output_path}'.")

def main():
    if len(sys.argv) < 2:
        print("Usage: python fastq_qc.py <path_to_fastq_or_fastq.gz> [output_report.html]")
        sys.exit(1)
        
    fastq_path = sys.argv[1]
    if not Path(fastq_path).exists():
        print(f"Error: FASTQ file '{fastq_path}' not found.")
        sys.exit(1)
        
    # Set output HTML path
    if len(sys.argv) >= 3:
        output_html = sys.argv[2]
    else:
        # Defaults to filename_qc_report.html in same directory
        base = os.path.basename(fastq_path)
        if base.endswith(".gz"):
            base = base[:-3]
        if base.endswith(".fastq"):
            base = base[:-6]
        elif base.endswith(".fq"):
            base = base[:-3]
        output_html = os.path.join(os.path.dirname(fastq_path) or ".", f"{base}_qc_report.html")
        
    print(f"Analyzing FASTQ file: {fastq_path}...")
    qc_data = run_qc(fastq_path)
    
    print("Generating report...")
    generate_html_report(qc_data, output_html)

if __name__ == "__main__":
    main()
