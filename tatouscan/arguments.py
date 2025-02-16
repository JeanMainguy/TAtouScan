import os
import logging

from typer import Argument, Option, echo, Exit

def version_callback(value: bool):
    """Prints the version and exits if --version is passed."""
    if value:
        echo(f"TAtouScan {__version__}")
        raise Exit()

args = {
    "faa_file" : Argument(help="Path to the input FAA file."),
    "gff_file" : Argument(help="Path to the input GFF file."),
    "hmm_db" : Argument(help="Path to the HMM profile database."),
    "output_directory" : Argument(help="Path to the output directory where results will be saved."),
    "e_value" : Argument(help="The E-value threshold for HMMER hits"),
    "max_sequence_length" : Argument(help="Maximum sequence length (between 30 to 500 amino acids) to consider a gene a putative Toxin or antitoxin."),
    "max_distance" : Argument(help="The maximum distance between paired toxin and antitoxin genes (between 0 to 300 nucleotides)."),
    "cpu" : Argument(help="Number of CPU threads to use for parallel processing."),
    "verbose_mode" : Argument(help="Control verbosity (e.g., INFO, DEBUG, WARNING) for troubleshooting."),
    "version" : Option("--version", callback=version_callback, help="Show the version and exit.")
}


logging_levels = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "ERROR": logging.ERROR
}


def validate_file_path(file_path: str, errors: list, argument: str):
    """Validate file path exists"""
    if not os.path.isfile(file_path):
        errors.append(f"{argument} : {file_path} is not a file path")


def validate_directory_path(dir_path: str, errors: list, argument: str):
    """Validate directory path exists"""
    if not os.path.isdir(dir_path):
        errors.append(f"{argument} : {dir_path} is not a directory path")


def validate_threshold(threshold, min_value: int, max_value: int, errors: list, argument: str):
    """Validate threshold argument is between his min and max tolerated value"""
    threshold_value = float(threshold)
    if (threshold_value < min_value) or (threshold_value > max_value):
        errors.append(f"{argument} : Get {threshold}, expected value between {min_value} and {max_value}")


def validate_cpu_number(cpu: int, errors: list):
    """Validate CPU number is between 1 and the total cpu count"""
    total_cpu = os.cpu_count()
    if (cpu < 1) or (cpu > total_cpu):
        errors.append(f"cpu : Get {cpu}, expected value between 1 and {total_cpu}")

def validate_verbose_mode(verbose_mode: str, errors: list):
    """Validate logging level"""
    if verbose_mode not in logging_levels.keys():
        errors.append(f"verbose_mode : {verbose_mode} is not a valid logging level")



def get_errors_from(faa_file: str, gff_file: str, hmm_db: str, output_directory: str, e_value: float, max_sequence_length: int, max_distance: int, cpu: int, verbose_mode: str):
    """Validate each argument independently and register errors"""
    errors = []
    validate_file_path(faa_file, errors, "faa_file")
    validate_file_path(gff_file, errors, "gff_file")
    validate_file_path(hmm_db, errors, "hmm_db")
    validate_directory_path(output_directory, errors, "output_directory")
    validate_threshold(e_value, 0, 1, errors, "e_value")
    validate_threshold(max_sequence_length, 30, 500, errors, "max_sequence_length")
    validate_threshold(max_distance, 0, 300, errors, "max_distance")
    validate_cpu_number(cpu, errors)
    validate_verbose_mode(verbose_mode, errors)
    return errors


def validate_args(faa_file: str, gff_file: str, hmm_db: str, output_directory: str, e_value: float, max_sequence_length: int, max_distance: int, cpu: int, verbose_mode: str):
    """Validate arguments and raise an exception with errors if there is any"""
    errors = get_errors_from(faa_file, gff_file, hmm_db, output_directory, e_value, max_sequence_length, max_distance, cpu, verbose_mode)
    if errors:
        raise Exception("Invalid arguments :\n\t- " + "\n\t- ".join(errors))