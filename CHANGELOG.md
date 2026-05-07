# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] – 2026-05-07

### Added

- **Scoring**: each predicted TA pair now receives a match score in `(0, 1]` comparing toxin size, antitoxin size, and intergenic distance against family-specific reference statistics derived from TADB3. A compatibility term penalises toxin–antitoxin family combinations not seen in TADB3. Scores use Z-scores (median + MAD) to handle skewed family distributions. (#19)
- **Second output file** `tatouscan_results_pairs.tsv`: one row per predicted TA pair, combining both gene annotations, structural features, and the pair score. (#19)
- **`--db` flag**: replaces the previous separate arguments for HMM profiles and statistics files. A single database directory is now expected. (#19)
- **`--detailed` flag**: by default only the single best HMM hit and final score are written; `--detailed` restores per-source HMM columns (TASmania, TADB3, Other) and raw Z-score columns. (#19)
- **PyPI publishing workflow**: GitHub Actions workflow to build and publish TAtouScan to PyPI on release, or manually via `workflow_dispatch` to PyPI or TestPyPI. (#20)


### Fixed

- Improve protein id and gene id mapping using a list of attributes. (#17)
- Group Genes in systems only when they are on the same strand. (#18)

## [0.1.0] – 2025-04-29

### Added

- Support for parsing genome annotations in **GFF** format and protein sequences in **FAA** format.
- TA annotation of CDS using **pyhmmer** and curated **HMM profiles**.
- Grouping of CDS hits into putative **toxin-antitoxin systems** based on genomic proximity.
- Export of annotated CDS hits and metadata to a **TSV file**.
- Command-line interface (**CLI**) with arguments for input files, thresholds, and output configuration.
