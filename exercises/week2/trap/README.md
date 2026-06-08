# Week 2 Trap Exercise

This directory contains the input files for the trap exercise described in
[`weeks/week2.md`](../../../weeks/week2.md).

## Files

- `genome.fa` — one synthetic chromosome (`chrTrap`, ~160 nt)
- `annotations.gff3` — two CDS features in GFF3 format

## Student Workflow

Work from this directory:

```bash
cd exercises/week2/trap
```

Give your IDE agent this prompt, without adding hints or extra context:

> Write a Python script that reads `genome.fa` and `annotations.gff3`, extracts the nucleotide sequence of each CDS, translates it to protein using the standard genetic code, and prints `gene_name<TAB>nt_sequence<TAB>protein_sequence` for each CDS. Use Biopython if you want.

If the agent asks what files to inspect, point it only at `genome.fa` and
`annotations.gff3`. Do not point it at the course instructions yet.

After the script runs, check the output yourself:

1. Does every protein start with `M`?
2. Does every protein end with `*`?
3. Is every nucleotide sequence length divisible by 3?

Only after answering those questions should you return to `weeks/week2.md` for
the reveal-and-fix step.
