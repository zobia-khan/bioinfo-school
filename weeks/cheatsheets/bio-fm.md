# Bio Foundation Model Cheat Sheet

Use this as a map of common model families, not as a leaderboard. Availability, APIs, and recommended checkpoints change quickly.

## Protein Sequence Models

| Model family | Typical use | Access |
| --- | --- | --- |
| ESM2 | Protein embeddings, variant scoring, quick baselines | HuggingFace, `fair-esm` |
| ESM-C | Protein embeddings with newer EvolutionaryScale models | EvolutionaryScale API / model docs |
| ESM3 | Protein generation and multimodal sequence-structure-function work | EvolutionaryScale |
| ProtT5 / ProtBert | Older but still useful protein embedding baselines | HuggingFace |

Watch for pooling choices. For sequence-level embeddings, mean-pooling over non-padding residues is usually a better first baseline than accidentally pooling padding tokens.

## Structure Prediction

| Tool | Typical use | Practical note |
| --- | --- | --- |
| AlphaFold2 / ColabFold | Single-chain and some multimer structure prediction | Best free-tier starting point |
| AlphaFold3 | Biomolecular complexes | Web server access and usage limits vary |
| Boltz-1 | Open structure prediction for complexes | GPU recommended |
| Chai-1 | Protein, ligand, and complex prediction | GPU recommended |

Confidence is part of the output. Always inspect pLDDT and, for complexes or multimers, PAE.

## Protein Design

| Model or tool | Typical use | Practical note |
| --- | --- | --- |
| RFdiffusion | Backbone generation and motif scaffolding | More engineering-heavy than Week 3 requires |
| ProteinMPNN | Sequence design for a fixed backbone | Often paired with generated or known structures |
| ESM inverse folding | Sequence design conditioned on backbone coordinates | Good notebook-based entry point |

Generated designs are hypotheses. They need structural, biophysical, and experimental validation.

## Genomic Language Models

| Model family | Typical use | Access |
| --- | --- | --- |
| Nucleotide Transformer | DNA embeddings and classification baselines | HuggingFace |
| DNABERT-2 | DNA sequence classification and embeddings | HuggingFace |
| HyenaDNA | Long-context DNA modeling | HuggingFace / GitHub |
| Evo | Long biological sequence generation and analysis | Arc Institute / model docs |

Check tokenization, context length, padding masks, and train/test splits before interpreting scores.

## Single-Cell Models

| Model family | Typical use | Access |
| --- | --- | --- |
| scGPT | Single-cell embeddings, perturbation tasks | GitHub / model docs |
| Geneformer | Gene and cell embeddings from expression profiles | HuggingFace / GitHub |

These models often depend on gene identifiers, normalization choices, and organism-specific preprocessing. Most failures come from mismatched input conventions rather than the model itself.

## Quick Validation Habits

- Compare against a simple baseline before trusting a foundation model.
- Inspect confidence outputs, not just predictions.
- Plot embeddings and check whether known labels behave sensibly.
- Read model cards for input limits and intended use.
- Record exact model names, versions, seeds, and hardware.
