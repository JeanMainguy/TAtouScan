
# TatouScan

TAtouScan is a command-line tool for identifying toxin-antitoxin (TA) systems in genomes and metagenomes.

## Installation

To install TAtouScan, first clone the repository:

```bash
git clone https://github.com/JeanMainguy/tatouscan.git
cd tatouscan
```
Create a virtual environment:
```bash
# Create a virtual environment
python -m venv venv

# Activate the virtual environment
# On Linux/macOS:
source venv/bin/activate

# On Windows (PowerShell):
# venv\Scripts\Activate

```
And then install TAtouScan using `pip`:
```bash
# Install the tool in editable mode
pip install -e .
```

## Usage

After installation, you can run the tool using:

```bash
tatouscan --help
```

## License

This project is licensed under the [MIT License](LICENSE).
