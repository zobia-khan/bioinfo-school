## Week 1

### From the materials

Karpathy Deep Dive: Pretraining is where 99% of the compute goes, turning web text into a base model that is just a highly advanced internet-mimic.

What I want to test: I want to see how a base model (if available) vs. a post-trained model handles completing a highly specific genomic sequence string when given no conversational context.

Karpathy Deep Dive: Tokenization is why LLMs struggle with discrete character-level tasks because chunks of characters are compressed into single integer tokens.

What I want to test: Test how different tokenizers (e.g., Tiktoken vs. Llama) chunk an identical string of amino acid residues to see if they split functional domains awkwardly.

Karpathy Deep Dive: Post-training (SFT and RLHF) changes the behavior and tone of the model from a document-completer to an assistant, but it doesn't fundamentally give it a "database" of perfect facts.

What I want to test: Prompt a model to write a Python script for reverse-complementing a sequence, first using standard chat, and then forcing it to output character by character to see if its logic changes.

Domain Expertise Lecture Note: * Example from my domain: In bioinformatics pipelines (such as extracting a coding sequence from genome coordinates) , a vanilla AI agent will confidently output running code like seq[start:end]. However, an expert knows a hidden trap: GFF/GTF coordinates are commonly 1-based and inclusive, whereas Python slicing is 0-based and end-exclusive. Without domain expertise to enforce biological invariants (e.g., verifying that the reading frame is preserved, start and stop codons are correct, and the sequence length is strictly divisible by 3) , the workflow falls into Trap 1: working code that produces a plausible output from a fundamentally wrong analysis. Execution is not validation ; domain judgment is the only apparatus that keeps the AI's accelerated generation connected to physical reality.

Reflection Exercise Part 1 (Miscounted residues in a sequence): * Why it's hard/easy: This is notoriously hard for an LLM because of tokenization. The model does not see "A-T-G-C" as individual characters; it might see "ATG" as token 3492 and "C" as part of another. When asked to count residues, it is guessing based on statistical associations of those token frequencies, not physically counting items in an array.

Domain check needed: Never trust an LLM's string length or indexing. Run a local len(sequence) script or use a regex pattern tool to get a precise character count.

Reflection Exercise Part 2 (Confused coordinate systems - 0-based vs 1-based): * Why it's hard/easy: The model has read thousands of papers and code repos using both BED format (0-based) and GTF/VCF format (1-based). Because both conventions are highly prevalent in its pretraining corpus, the next-token distribution for a coordinate calculation can easily bleed from one convention into the other, creating a devastating off-by-one error.

Domain check needed: Explicitly assert the file type and coordinates in a small test suite. For example, grab a known feature (like a specific Exon 1 boundary) and check if the generated code correctly hits the exact nucleotide index in a controlled local environment.

Reflection Exercise Part 3 (Classifier learning hospital/site/batch instead of biology): * Why it's hard/easy: Neural networks are fundamentally lazy; they optimize for the path of least resistance to minimize loss during training. If a synthetic biomedical dataset has a batch effect (e.g., all cancer samples came from Lab A on an Illumina sequencer and all controls came from Lab B on a Pacific Biosciences sequencer), the network will classify the metadata/instrument noise rather than the biological signal, yielding a deceptively high metric.

Domain check needed: Perform a principal component analysis (PCA) on the embeddings or features to visually check for clustering by batch/site. Run the model on an independent, strictly external validation dataset from a completely separate institution to see if the high performance collapses.

---

### Surprises

2026-06-08 · Vanilla Chatbot (No Code Execution)

What I asked: "Write a Python snippet to extract the transcript ID and gene name from this GTF line: chr1 HAVANA transcript 11869 14409 . + . gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";"

What happened: The model generated a regex string re.search(r'transcript_id "([^"]+)"', gtf_line), re.search(r'gene_name "([^"]+)"', gtf_line). While it successfully found the fields, it assumed every GTF file strictly uses double quotes and trailing semicolons in that exact spatial order, failing to account for rows that might have extra spaces or missing optional attributes. It did this purely by predicting the most common text pattern it had memorized.

Takeaway: Vanilla text generation defaults to the most generic text structure and fails to validate edge cases.

2026-06-08 · Chatbot with Code Execution (Gemini with Code Execution)

What I asked: "Write and run a Python snippet to extract the transcript ID and gene name from this GTF line: chr1 HAVANA transcript 11869 14409 . + . gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";"

What happened: The model did not just guess the regex. It silently wrote a Python script, spun up a local sandbox environment, ran the example line through its code, saw that it correctly extracted ENST00000456328 and DDX11L1, and then presented the working code to me. The gap was immense; the code execution tool served as a real-time sanity check against its own token-prediction hallucinations.

Takeaway: Code execution forces the LLM to pass its own syntax/runtime validation filter

Scientific validation check: To make sure it's correct for the right reasons, I need to test a negative case: what happens if transcript_id is missing from the attributes column (e.g., a structural RNA line)? The chatbot's script would throw an AttributeError on .group(1). A robust scientific pipeline requires explicit exception handling for missing attributes, which the agent skipped because it only optimized for my single positive example.


## Week 2
### Surprises

The agent did not make the expected mistake, it wrote the correct slicing logic in the script by looking at annotations.gff3 and genome.fa files.

start_idx = cds["start"] - 1
end_idx = cds["end"]


### From the materials
1. Other "Looks Right But Isn't" Failures in AI-Generated Code
Beyond strand, coordinate, and format confusion, agents frequently stumble over the subtle, domain-specific logic of biology:

Chromosome Prefix Mismatches: Forgetting to harmonize chr1 versus 1 across different annotation files, which silently drops massive amounts of data during table joins without throwing an error.

Ignoring Reverse Complements: Writing motif-searching logic that only scans the forward strand, completely missing valid biological targets in unstranded data.

Incorrect Background Sets: Calculating enrichment metrics against a universe of all annotated genes rather than restricting the background to only the genes successfully detected or expressed in that specific assay.

2. Biological Invariants for Routine Validation (Genomics/Transcriptomics)
When validating an agent's logic, anchor your tests to immutable biological rules. For transcriptomics, three reliable invariants are:

The Modulo-3 Rule: If an agent extracts intact, translatable Coding Sequences (CDS), the length of every extracted sequence must be perfectly divisible by 3.

Coordinate Logic: A sequence's start coordinate must mathematically always be less than or equal to its end coordinate, regardless of whether the gene lies on the positive or the negative strand.

Alphabet Integrity and Arithmetic: Sequence features must strictly consist of valid characters (A, C, G, T, N). Furthermore, if the agent calculates fractional compositions (like GC content), those fractions must mathematically sum exactly to 1.0.

3. Scaling Validation for 10,000 CDS Features
You cannot visually inspect 10,000 rows of data. To scale validation safely, you have to transition from manual reading to programmatic checking:

Automated Assertion Pipelines: Write strict unit tests for the invariants mentioned above (e.g., assert no negative lengths, assert expected data types, assert values bounded between 0 and 1) so the computer automatically flags violations across the entire dataset.

Distribution Sanity Checks: Plot global distributions (histograms, box plots) for every engineered feature. You don't need to read every row to visually spot an impossible outlier (like a GC content of 110%) or a bizarrely bimodal distribution.

"Golden" Control Subsets: Run the agent's code on a small, hand-curated dataset of 50-100 sequences where you already know the exact biological ground truth, and mandate a 100% match before scaling up to the full 10,000.

Stress-Testing Edge Cases: Sample a small fraction of the data, but heavily over-index on known biological edge cases (e.g., genes with overlapping exons, exceptionally short transcripts, or sequences peppered with 'N' bases) to see how the agent's logic handles exceptions.


## Week 3

### From the materials

**Jumper Nobel Lecture (AlphaFold):** The key insight I want to test: pLDDT is not a generic confidence score — it is a predicted per-residue accuracy relative to an idealised structure. A high pLDDT region doesn't mean that region is biologically ordered in the cell; it means the model is sure about its geometry. Disordered regions (signal peptides, IDRs) reliably score below 50. I want to check this on a known IDP (intrinsically disordered protein) to see whether the score correctly tracks experimental disorder data.

**CARBON tech report:** CARBON presents task-specific benchmarks but all within its own evaluation suite. The claim I'd most want to verify on my own data: do the CARBON embeddings outperform ESM-2 embeddings for family classification when the test proteins are genuinely out-of-distribution (different species or synthetic sequences)? The in-distribution benchmarks look strong, but they can mask brittleness.

**Interpretability lecture (Stefan):** A clean cluster in a UMAP plot is not proof of biological mechanism — it proves the model learned something that correlates with the family label. The same cluster could reflect sequence length bias, amino acid composition, or training data coverage rather than true functional similarity. One thing I'd test first: replace protein embeddings with random vectors of the same dimensionality and see whether UMAP can still produce "clusters" by chance with only 45 points. That sets the baseline for what structured separation actually means.


---

### Surprises

**2026-06-20 · Antigravity agent + ESM-2 notebook (B_protein_embeddings_esm2_output.ipynb)**


What happened: (45 sequences, 320-dim CLS embeddings, UMAP shape 45×2), identified the five protein families from the TSV, and produced a factual write-up that correctly described which families clustered and why kinases were spread — all without running any code itself. It also embedded the saved UMAP plot. Impressive for structured document extraction, but it could not check whether the UMAP axes were meaningful or whether a different random seed would have changed the topology.

---

## Week 4

### From the materials

**Tim Berglund — Agent Skills vs MCP (~10 min):** The key distinction: a *skill* is a bundle of instructions and scripts that lives inside the agent's own context and gets executed in its environment; an *MCP tool* is a typed, callable function exposed over a protocol, where the agent does not own the implementation. Skills are better for stable conventions you keep repeating (e.g., "always use `uv` and Phred-33 quality thresholds"). MCP tools are better for live data sources, external services, or anything where the implementation should be maintained independently (e.g., a database query against Ensembl or a UniProt lookup). The overlap is real — anything simple enough can be either — but the question to ask is: *"does the logic belong to the agent, or to a service?"*

---

### Surprises

**2026-06-20 · Antigravity agent — AGENTS.md and repo summary**

Asked the agent: *"Summarise the repo."* After `AGENTS.md` was committed it correctly led with the Python/uv constraint, flagged the Phred-33 convention, listed the files it must not commit, and explicitly mentioned the validation checks. Before `AGENTS.md` existed, the same prompt produced a generic directory listing with no domain constraints at all. The difference was immediate and required no change to the prompt itself.

Takeaway: `AGENTS.md` is genuinely load-bearing. A 20-line conventions file changes agent behaviour more reliably than a long system prompt, because it lives in version control and applies to every new conversation automatically.
