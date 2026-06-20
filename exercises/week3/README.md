# Week 3 Exercise Starters

This folder contains small starter assets for the Week 3 foundation-model exercises.

## Exercise A: Structure Prediction

Use this public notebook:

- [ColabFold AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

Keep inputs small. The Week 3 page suggests 1-3 proteins or domains under roughly 200 amino acids.

## Exercise B: Protein Embeddings

Build the starter FASTA from UniProt:

```bash
python3 fetch_proteins.py
```

This reads `protein_accessions.tsv` and writes `proteins.fasta`. The list is grouped by coarse protein family so you can check whether embedding-space clusters recover the labels.

In Colab, use [`B_protein_embeddings_esm2.ipynb`](./B_protein_embeddings_esm2.ipynb).

Before running the notebook, upload these two files into the Colab session:

- `~/projects/other/bioinfo-school/exercises/week3/protein_accessions.tsv`
- `~/projects/other/bioinfo-school/exercises/week3/proteins.fasta`

Colab runs on a Google server, not on your machine, so the notebook cannot see local files until you upload them.

## Exercise C: Optional Genomic Benchmarks

Use a fresh notebook or script. These links are useful starting points:

- [Genomic Benchmarks repository](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks)
- [Nucleotide Transformer v2 50M multi-species model card](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-50m-multi-species)
- [HuggingFace Transformers notebooks](https://huggingface.co/docs/transformers/notebooks)

Write your final numbers and short interpretation in `results.md`.
