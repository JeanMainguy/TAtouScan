# TAtouScan

**TAtouScan** is a command-line tool designed to identify **toxin-antitoxin (TA) systems** in genomes and metagenomes. 


## Installation

### Option 1: Install with `pip` 

1. **Clone the repository:**

```bash
git clone https://github.com/JeanMainguy/TAtouScan.git
cd TAtouScan
```

2. **Create and activate a virtual environment:**

```bash
# Create a virtual environment
python -m venv venv

# Activate the virtual environment
# On Linux/macOS:
source venv/bin/activate

```

1. **Install TAtouScan:**

```bash
pip install .
```

---

### Option 2: Install using `conda`

If you prefer using `conda`, you can create a dedicated environment as follows:

```bash
# Create a new conda environment with Python
conda create -n tatouscan python=3.12

# Activate the environment
conda activate tatouscan

# Clone the repository
git clone https://github.com/JeanMainguy/TAtouScan.git
cd TAtouScan

# Install TAtouScan
pip install -e .
```
> [!NOTE]
> TAtouScan is not yet available via `bioconda`. The above method combines `conda` for environment management and `pip` for installation.


### Download the TAtouScan Database

TAtouScan requires a database directory containing HMM profiles and reference statistics.


Download the database and extract it with:
```bash
wget https://zenodo.org/records/20059258/files/tatouscan_db.tar.gz
tar -xzf tatouscan_db.tar.gz
```

The database directory must contain the following four files:

```
tatouscan_db/
  ta.hmm                 # HMM profiles (HMMER3 format)
  hmm_info.tsv           # profile metadata (name, type, source)
  family_statistics.tsv  # per-family reference statistics for scoring
  known_pairs.tsv        # known toxin–antitoxin family co-occurrences
```


## Usage

After installation and downloading the database, run TAtouScan with:

- a **GFF** file with gene annotations
- a **FAA** file with the corresponding protein sequences
- the **database directory** downloaded above

```bash
tatouscan --gff <genes.gff> --faa <proteins.faa> --db tatouscan_db/
```

By default, results are written to a directory called `tatouscan_results/`. Use `--outdir` to specify a different location:

```bash
tatouscan --gff <genes.gff> --faa <proteins.faa> --db tatouscan_db/ --outdir my_results/
```

Two TSV files are produced inside the output directory:

| File | Description |
|------|-------------|
| `tatouscan_results.tsv` | One row per predicted toxin or antitoxin gene |
| `tatouscan_results_pairs.tsv` | One row per predicted TA pair (two-gene systems only) |



## HMM Database Composition

The HMM database used by TAtouScan is composed of profiles collected from multiple sources, including curated databases and literature. The file `tatouscan_hmm_description.tsv` provides metadata for each profile, indicating its origin and whether it corresponds to a **toxin** or an **antitoxin**.

### Breakdown of the database:

- **682 profiles** were obtained from the [TASmania project](https://doi.org/10.1371/journal.pcbi.1006946):  
  > Akarsu H, Bordes P, Mansour M, Bigot D-J, Genevaux P, Falquet L (2019). *TASmania: A bacterial Toxin-Antitoxin Systems database*. PLoS Comput Biol 15(4): e1006946.  
  > https://doi.org/10.1371/journal.pcbi.1006946

- **3,168 profiles** were generated from sequences in the [TADB 3.0 database](https://bioinfo-mml.sjtu.edu.cn/TADB3/):  
  These sequences were first clustered, and each cluster was then aligned using multiple sequence alignment. HMM profiles were built from the resulting alignments.

  > Guan J, Chen Y, Goh YX, Wang M, Tai C, Deng Z, Song J, Ou HY (2024).  
  > *TADB 3.0: an updated database of bacterial toxin-antitoxin loci and associated mobile genetic elements.*  
  > Nucleic Acids Research, 52(D1): D784–D790.  
  > https://doi.org/10.1093/nar/gkad962

- Additional HMM profiles were manually collected from other sources in the literature.


## Output

TAtouScan writes two TSV files into the output directory.

By default, only the most informative columns are written. Add `--detailed` to include per-source HMM breakdowns and raw Z-score columns.

### `tatouscan_results.tsv` — per-gene results

One row per predicted toxin or antitoxin gene.

| Column | Description |
|--------|-------------|
| `contig_name` | Contig where the gene is located |
| `gene_id` | Gene identifier (from the input GFF) |
| `start` / `end` | Genomic coordinates |
| `strand` | `+` or `-` |
| `length_aa` | Protein length in amino acids |
| `product` | Predicted gene product (if available) |
| `ta_system_id` | ID shared by both genes of a pair (`None` for single-gene predictions) |
| `is_single_gene` | `True` if no paired partner was found |
| `gene_type` | `Toxin` or `Antitoxin` |
| `hmm_name` / `hmm_score` / `hmm_evalue` | Best HMM hit across all database sources |
| `hmm_source` | Database the best hit comes from (`TADB3`, `TASmania`, or other) |
| `hmm_description` | Profile description |
| `pair_is_known` | `1` if this toxin–antitoxin family combination is known in TADB3, `0` if not, `None` if family could not be identified |
| `score` | Unified match score in `(0, 1]` (see [Scoring](#scoring)) |

Scoring columns are `None` for single-gene predictions.

### `tatouscan_results_pairs.tsv` — per-pair results

One row per predicted toxin–antitoxin pair. For systems with more than one toxin or antitoxin, all valid combinations are written as separate rows.

| Column | Description |
|--------|-------------|
| `ta_system_id` | Shared system ID (matches the per-gene file) |
| `contig_name` | Contig where the pair is located |
| `toxin_gene_id` | Toxin gene identifier |
| `toxin_strand` | `+` or `-` |
| `toxin_product` | Predicted gene product |
| `toxin_length_aa` | Toxin protein length in amino acids |
| `toxin_hmm_name` / `_score` / `_evalue` / `_source` / `_description` | Best HMM hit for the toxin |
| `antitoxin_gene_id` | Antitoxin gene identifier |
| `antitoxin_strand` | `+` or `-` |
| `antitoxin_product` | Predicted gene product |
| `antitoxin_length_aa` | Antitoxin protein length in amino acids |
| `antitoxin_hmm_name` / `_score` / `_evalue` / `_source` / `_description` | Best HMM hit for the antitoxin |
| `intergenic_distance` | Distance in nucleotides between the two genes (negative = overlap) |
| `pair_is_known` | `1` / `0` / `None` (see above) |
| `score` | Unified match score in `(0, 1]` |

### Detailed output

With `--detailed`, the following additional columns are written to both files:

- **Per-source HMM hits**: `TASmania_hmm_name/score/evalue/description`, `TADB3_hmm_name/score/evalue/description`, `Other_hmm_name/score/evalue/description` (prefixed with `toxin_` / `antitoxin_` in the pairs file)
- **Raw Z-scores**: `toxin_size_z`, `at_size_z`, `intergenic_distance_z`, `matched_family`, `n_reference_pairs`

The pairs file also adds `toxin_start/end` and `antitoxin_start/end` in detailed mode.

---

## Scoring

Every predicted TA **pair** is compared against reference statistics derived from known TADB3 type-II systems. The score measures how closely the predicted pair resembles a genuine TA system of its family.

### What is compared

Three structural features are measured for each predicted pair and compared against the reference distribution for the matched family:

| Feature | Definition |
|---------|------------|
| `toxin_size` | Toxin protein length (amino acids) |
| `at_size` | Antitoxin protein length (amino acids) |
| `intergenic_distance` | Distance in nucleotides between the two genes (negative = overlap) |

The **toxin family** is determined from its best TADB3 HMM hit. If no TADB3 hit exists or the family has fewer than 20 reference pairs, global statistics computed across all families are used as a fallback.

### Robust Z-scores

For each feature, a Z-score measures how far the predicted value deviates from the family reference:

$$z = \frac{x - \text{median}}{\text{MAD} / 0.6745}$$

Median and MAD (median absolute deviation) are used instead of mean and standard deviation because size distributions in TA families are often skewed. This makes the scores robust to outliers.

### Unified score

All Z-scores are combined into a single **score** in the range $(0, 1]$:

$$\text{score} = \exp\!\left(-\frac{1}{n}\sum_i |z_i|\right)$$

The mean is taken over all available terms: the three structural Z-scores plus a **compatibility term** ($z_{\text{compat}}$) based on whether this toxin–antitoxin family combination has been observed in TADB3:

- `pair_is_known = 1` → $z_{\text{compat}} = 0$ (no penalty)
- `pair_is_known = 0` → $z_{\text{compat}} = 2$ (unknown combination lowers the score)
- `pair_is_known = None` → compatibility term excluded from the mean

**Score interpretation:**

| Score | Meaning |
|-------|---------|
| ~1.0 | Features match the family reference almost exactly, known combination |
| ~0.7 | Moderate structural match, known combination |
| ~0.4 | Moderate structural match, but family combination not seen in TADB3 |
| < 0.2 | Large structural deviations or unknown combination — treat with caution |

A high score supports a genuine TA pair; a low score does not exclude it, but suggests the prediction should be reviewed.



## License

This project is licensed under the [MIT License](LICENSE).
