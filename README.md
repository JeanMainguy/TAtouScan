# TAtouScan

**TAtouScan** is a command-line tool designed to identify **toxin-antitoxin (TA) systems** in genomes and metagenomes. 


## Installation

### Option 1: Install with `pip` 

1. **Clone the repository:**

```bash
git clone https://github.com/JeanMainguy/tatouscan.git
cd tatouscan
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
git clone https://github.com/JeanMainguy/tatouscan.git
cd tatouscan

# Install TAtouScan
pip install -e .
```
> [!NOTE]
> TAtouScan is not yet available via `bioconda`. The above method combines `conda` for environment management and `pip` for installation.


### Download TAtouScan HMM Database

TAtouScan requires a database of HMM profiles to run. You can download the latest version from Zenodo:  
🔗 https://zenodo.org/records/15305313

Download the required files using:

```bash
wget https://zenodo.org/records/15305313/files/tatouscan_hmm_description.tsv
wget https://zenodo.org/records/15305313/files/tatouscan_hmm_profiles.hmm.gz
```


## Usage

After installation and downloading the required HMM database, you can run TAtouScan as follows:

TAtouScan currently requires:
- a **GFF** file with gene annotations
- a **FAA** file with the corresponding protein sequences

To identify toxin-antitoxin systems in a genome, run:

```bash
tatouscan --gff <genes.gff> --faa <proteins.faa> \
  --hmm_db tatouscan_hmm_profiles.hmm.gz \
  --hmm_info tatouscan_hmm_description.tsv
```

This command will produce an output file named:  
📄 `tatouscan_results.tsv` — listing all predicted toxins and antitoxins found in the input genome.



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

TatouScan produces a TSV file (`tatouscan_results.tsv`) summarizing the predicted toxin-antitoxin (TA) genes. The file includes the following columns:

| Column Name                | Description                                                            |
| -------------------------- | ---------------------------------------------------------------------- |
| `contig_name`              | Name of the contig where the gene is located                           |
| `gene_id`                  | Unique identifier of the gene (from the input GFF file)                |
| `start`                    | Start position of the gene on the contig                               |
| `end`                      | End position of the gene on the contig                                 |
| `strand`                   | Strand of the gene (`+` or `-`)                                        |
| `product`                  | Predicted function or product of the gene (if available)               |
| `is_single_gene`           | Whether the gene is a single hit or part of a TA pair (`True`/`False`) |
| `ta_system_id`             | ID of the TA system this gene belongs to (shared between paired genes) |
| `gene_type`                | Type of gene: either `toxin` or `antitoxin`                            |
| `TASmania_hmm_name`        | Name of the matched HMM profile from the TASmania database (if any)    |
| `TASmania_hmm_score`       | Bit score of the TASmania HMM hit                                      |
| `TASmania_hmm_evalue`      | E-value of the TASmania HMM hit                                        |
| `TASmania_hmm_description` | Description of the TASmania HMM profile                                |
| `Other_hmm_name`           | Name of a matched HMM profile from other sources (if any)              |
| `Other_hmm_score`          | Bit score of the "Other" HMM hit                                       |
| `Other_hmm_evalue`         | E-value of the "Other" HMM hit                                         |
| `Other_hmm_description`    | Description of the HMM profile from other sources                      |
| `TADB3_hmm_name`           | Name of the matched HMM profile from the TADB3 database (if any)       |
| `TADB3_hmm_score`          | Bit score of the TADB3 HMM hit                                         |
| `TADB3_hmm_evalue`         | E-value of the TADB3 HMM hit                                           |
| `TADB3_hmm_description`    | Description of the TADB3 HMM profile                                   |



## License

This project is licensed under the [MIT License](LICENSE).
