# LocusRuler

LocusRuler determines whether a reference gene cluster is present in each genome of a target set, and characterizes how each genome differs from that reference. It is intended for comparative surveys in which the cluster of interest varies across the genomes examined. In these surveys, telling genuine biological difference apart from technical artifact carries real interpretive weight.

Three properties distinguish it from a presence-or-absence homology search.

**Assembly of clusters from dispersed pieces.** LocusRuler reconstructs a cluster even when a contig boundary interrupts it, or when it is genuinely split across the chromosome. It groups matches by contig, strand and physical proximity in the target genome, then reconciles the separated pieces into a single call. Recurring arrangements are cataloged as distinct structures.

**Weighting of genes by functional role.** LocusRuler derives a cluster's status from whether each of its genes is intact, degraded, or missing, rather than from one overall similarity score. You may designate a gene as essential to biosynthesis, or as merely part of the cluster without being required, and LocusRuler classifies the cluster's status accordingly. This distinction is also what separates a cluster that is complete but diverged from one that is genuinely decaying, an evolutionary distinction a single similarity cutoff would conflate.

**Classification by synteny context.** LocusRuler records which reference flanking genes are adjacent to each accepted piece in the target genome, producing a per-genome fingerprint of the surrounding region. This distinguishes a cluster in its canonical context from one that has been translocated or rearranged, and separates a flanking gene missing from the assembly from one genuinely absent from the genome. Where no reference anchor is adjacent, LocusRuler records the nearest annotated product instead, so genomes sharing an unfamiliar neighbor can still be grouped together.

## Which genomes carry a given gene cluster?

1. Install. Requires Python 3.11 or later and NCBI BLAST+ (`makeblastdb`, `blastn`, `tblastn`, `blastp`) on PATH. HMMER (`hmmsearch`, `hmmfetch`) is needed only where `[domain_recovery]` is enabled.

   ```bash
   python -m pip install -e .
   ```

2. Start the wizard with a word from a gene of interest, or with a locus_tag.

   ```bash
   locus-ruler-wizard transporter
   ```

3. Follow the prompts. Select a reference genome, select the two genes marking the ends of the cluster, then run. Nothing needs to exist beforehand, as the wizard offers to create `settings.toml` and to download a genome dataset from NCBI where neither is present.

INPUT

* a target genome set, which the wizard can download, or see [How is a genome dataset built?](#how-is-a-genome-dataset-built)
* the two genes marking the ends of the cluster, within one reference genome

OUTPUT

* `output/<target>/<locus_id>/locus_report.xlsx`, containing one row per genome. The status is `COMPLETE` where the cluster is intact, `CONDITIONAL` where it is intact but missing a replaceable component, `DIVERGENT` where it is present but diverged, `DECAYED` where it is partially lost, `ABSENT` where it is not detected, and `UNKNOWN` otherwise
* `OUTPUTS.md` in the same directory, regenerated with every run, documenting every column and cell label produced by that run
* `cluster_heatmap.png` and `cassette_structure.png`, the cluster and cassette figures

## How is the pipeline run with a known target and flank genes?

The same sequence, without prompts.

```bash
python -m pip install -e .
cp settings.example.toml settings.toml

locus-ruler-setup --taxon "Escherichia coli" \
    --outdir input/example --db-name example --add-target settings.toml

locus-ruler-genes --settings settings.toml --target example \
    --accession GCF_000000000.1 --out example_genes.csv

locus-ruler-make-config --locus-id my_locus --accession GCF_000000000.1 \
    --target example --flank-l LEFT_TAG --flank-r RIGHT_TAG \
    --table example_genes.csv --out my_locus.json

locus-ruler --config my_locus.json --settings settings.toml \
    --target example --cpu 8
```

`locus-ruler-genes` exports every gene in the reference genome, which allows the two genes bracketing the cluster to be identified. `locus-ruler-make-config` converts those two locus_tags into the configuration file that `locus-ruler` reads.

## How is a genome dataset built?

```bash
locus-ruler-setup \
  --taxon "Genus species" \
  --outdir input/example \
  --db-name example \
  --seq-count 100 \
  --assembly-level complete chromosome scaffold contig \
  --add-target settings.toml
```

INPUT

* `--taxon`, an NCBI taxon name
* `--seq-count`, the greatest number of contigs a genome may contain before it is excluded, where `0` imposes no limit
* `--assembly-level`, the NCBI assembly levels to accept
* `--dry-run`, which reports how many genomes a given combination returns before anything is downloaded

OUTPUT

* `input/<name>/`, containing a SQLite database together with `gff/`, `genomes/` and `combined_proteins.faa`
* a `[[targets]]` block appended to `settings.toml`, where `--add-target` is supplied

## How is the anchor table curated?

On the first run LocusRuler writes `<locus_id>_anchors.csv` alongside the configuration file. Its `status_role` column sets how much each gene's state counts toward the cluster-level call. `CORE` genes define the cluster and drive the principal calls. `SHARED` genes represent functions that may legitimately be supplied from elsewhere in the genome, such as a precursor produced by an unrelated pathway, so their absence demotes an otherwise intact cluster to `CONDITIONAL` rather than declaring it degraded. `ASSOCIATED` and `CONTEXT` genes are reported without affecting the status, and `IGNORE` genes are excluded entirely. The remaining editable columns are `family`, which groups paralogous genes so that duplicates are counted as one system, `exception`, which keeps a mobile element visible without allowing it to influence cluster assembly, and `lenient`, which relaxes the matching thresholds for one manually verified divergent gene.

Two further columns cover a gene that has diverged past the point where BLAST will find it. `pfam` names the Pfam domains its protein must carry, written as one accession such as `PF00005`, or as several joined by `+` where all are required. Where BLAST returns nothing for that gene, LocusRuler searches the genome for a protein carrying those domains and reports what it finds. Setting `pfam_split` to `TRUE` additionally accepts a protein carrying only some of the required domains, provided the gene beside it carries the rest, which is how a coding sequence broken in two by a frameshift is recovered as one gene rather than lost. This layer needs HMMER and a local copy of Pfam-A, which `locus-ruler-pfam` downloads. The search is cut at each family's own gathering threshold, the same curated cutoff Pfam applies when deciding whether a protein carries a domain, so the call does not depend on a chosen E-value or on how many genomes were searched.

The wizard presents all of these choices one gene at a time, with an explanation of each. The file may also be edited directly.

## Where are results written?

```text
output/<target>/<locus_id>/
  locus_report.xlsx
  cluster_heatmap.png
  cassette_structure.png
  OUTPUTS.md
  tables/
  diagnostics/
```

`OUTPUTS.md` is regenerated with every run and documents every file, column and cell label that run produced. It, rather than this README, is the appropriate reference for interpreting a specific result.

Cached BLAST results are held under `output/_work/<target>/<locus_id>/`. Both `output/` and `input/` are excluded from version control.

## How is the cassette defined?

A gene joins the cassette when its piece was accepted into `pieces.csv` and its diagnostic zone is `cluster`. The rule is the same for single-contig and multi-contig loci, and accepted sibling pieces are included automatically in reference-query order. A homolog outside the accepted pieces is left out and reported separately.

Each arrangement is identified by `structure_signature`, which lists every cluster gene in reference-query order, and by `structure_id`, a deterministic hash of that signature. `assembly_signature` marks accepted piece boundaries with `||`, so fragmentation can be audited without changing `structure_id`.

An optional `cassette_display` block in the locus config sets figure colors, labels, scale, and connector styling through `family_styles` and `product_styles`. It affects the figure only, never membership or identifiers.

## Repository Conventions

The following apply to maintenance of this repository rather than to the running of an analysis.

* Downloaded data belongs under `input/`, and generated results under `output/`.
* Commit `settings.example.toml`. A personal `settings.toml` should never be committed.
* Commit a locus configuration before its first run rather than after, so that it remains free of machine-specific paths.
* Large Pfam and BLAST databases should not be committed.

## Running The Tests

```bash
python -m pip install -e ".[test]"
pytest tests
```

Neither a BLAST installation nor downloaded genomes are required. The suite runs against small built-in fixtures.

## Citation

No formal citation has been assigned.

