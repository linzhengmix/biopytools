import os
import pandas as pd
import subprocess
import gzip
import shutil  # Add this import


def read_bacteria_info(tsv_file):
    """Read bacterial genome information from a TSV file.

    Args:
        tsv_file: Path to the TSV file.

    Returns:
        DataFrame containing genome information.
    """
    return pd.read_csv(tsv_file, sep='\t')


def download_genome_wget(ftp_url, output_dir, filename, log_file):
    """Download bacterial genome file using wget.

    Args:
        ftp_url: FTP URL.
        output_dir: Directory to save the downloaded file.
        filename: Name of the downloaded file.
        log_file: Path to the wget log file.

    Returns:
        Path to the downloaded file.
    """
    local_filepath = os.path.join(output_dir, filename)
    cmd = ['wget', '-O', local_filepath, ftp_url, '-a', log_file]
    subprocess.run(cmd, check=True)
    return local_filepath


def check_file_integrity(filepath):
    """Check if the downloaded file is complete.

    Args:
        filepath: Path to the file.

    Returns:
        Boolean indicating whether the file is complete.
    """
    try:
        with gzip.open(filepath, 'rb') as f:
            while f.read(1024):
                pass
        return True
    except Exception as e:
        print(f"File integrity check failed for {filepath}: {e}")
        return False


def update_download_lists(complete_file, incomplete_file, complete_downloads, incomplete_downloads):
    """Update the lists of completed and incomplete downloads.

    Args:
        complete_file: Path to the completed downloads list file.
        incomplete_file: Path to the incomplete downloads list file.
        complete_downloads: Set of completed downloads.
        incomplete_downloads: Set of incomplete downloads.
    """
    with open(complete_file, 'w') as f:
        f.write('\n'.join(complete_downloads))

    with open(incomplete_file, 'w') as f:
        f.write('\n'.join(incomplete_downloads))


def update_faa_list(faalist_file, faa_list):
    """Update the list of Assembly Accessions with faa.gz files.

    Args:
        faalist_file: Path to the faa.gz files list file.
        faa_list: Set of Assembly Accessions with faa.gz files.
    """
    with open(faalist_file, 'w') as f:
        f.write('\n'.join(faa_list))


def update_failed_list(failed_file, failed_downloads):
    """Update the list of failed downloads.

    Args:
        failed_file: Path to the failed downloads list file.
        failed_downloads: Set of failed downloads.
    """
    with open(failed_file, 'w') as f:
        f.write('\n'.join(failed_downloads))


def scan_existing_files(directory, extension):
    """Scan existing files in a directory and check their integrity.

    Args:
        directory: Directory path.
        extension: File extension.

    Returns:
        Set of existing and complete files.
    """
    existing_files = set()
    for filename in os.listdir(directory):
        if filename.endswith(extension):
            filepath = os.path.join(directory, filename)
            if check_file_integrity(filepath):
                accession = filename.split('_')[1]
                existing_files.add(accession)
    return existing_files


def download_and_check_file(assembly_accession, assembly_name, file_type, output_dir, log_file, download_list, failed_downloads):
    """Download and check the integrity of a file.

    Args:
        assembly_accession: Assembly Accession of the genome.
        assembly_name: Assembly Name of the genome.
        file_type: File type (genomic or protein).
        output_dir: Directory to save the file.
        log_file: Path to the wget log file.
        download_list: Set of downloaded files.
        failed_downloads: Set of failed downloads.

    Returns:
        Boolean indicating whether the download and integrity check were successful.
    """
    ftp_path = f'/genomes/all/{assembly_accession[:3]}/{assembly_accession[4:7]}/{assembly_accession[7:10]}/{assembly_accession[10:13]}/{assembly_accession}_{assembly_name}'
    filename = f'{assembly_accession}_{assembly_name}_{file_type}.fna.gz' if file_type == 'genomic' else f'{assembly_accession}_{assembly_name}_{file_type}.faa.gz'
    ftp_url = f'https://ftp.ncbi.nlm.nih.gov{ftp_path}/{filename}'

    file_path = os.path.join(output_dir, filename)

    # Check if the file is already downloaded and complete
    if os.path.exists(file_path):
        if not check_file_integrity(file_path):
            os.remove(file_path)
    if not os.path.exists(file_path):
        for attempt in range(3):
            try:
                print(f"Downloading {filename} from {ftp_url} (Attempt {attempt + 1})")
                download_genome_wget(ftp_url, output_dir, filename, log_file)

                if check_file_integrity(file_path):
                    print(f"Download and integrity check passed for {filename}")
                    download_list.add(assembly_accession)
                    return True
                else:
                    print(f"Download or integrity check failed for {filename}")
            except subprocess.CalledProcessError:
                print(f"Failed to download {filename} (Attempt {attempt + 1})")
        else:
            failed_downloads.add(assembly_accession)
            return False
    return True


def decompress_gz_file(gz_file, output_dir):
    """Decompress a .gz file.

    Args:
        gz_file: Path to the .gz file.
        output_dir: Directory to save the decompressed file.

    Returns:
        Path to the decompressed file.
    """
    decompressed_file = os.path.join(output_dir, os.path.basename(gz_file).replace('.gz', ''))
    with gzip.open(gz_file, 'rb') as f_in:
        with open(decompressed_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return decompressed_file


def create_blast_db(fasta_file, db_type='nucl'):
    """Create a BLAST database from a FASTA file.

    Args:
        fasta_file: Path to the FASTA file.
        db_type: Type of the database ('nucl' for nucleotide, 'prot' for protein).

    Returns:
        None
    """
    db_name = os.path.splitext(fasta_file)[0]
    if not os.path.exists(f"{db_name}.nhr") and not os.path.exists(f"{db_name}.phr"):
        cmd = ['makeblastdb', '-in', fasta_file, '-dbtype', db_type]
        subprocess.run(cmd, check=True)
    else:
        print(f"BLAST database for {fasta_file} already exists. Skipping creation.")


def download_data(bacteria_info, genome_dir, protein_dir, log_file, genome_complete_downloads, protein_complete_downloads, failed_downloads):
    """Download genome and protein data.

    Args:
        bacteria_info: DataFrame containing bacterial genome information.
        genome_dir: Directory to save genome files.
        protein_dir: Directory to save protein files.
        log_file: Path to the wget log file.
        genome_complete_downloads: Set of completed genome downloads.
        protein_complete_downloads: Set of completed protein downloads.
        failed_downloads: Set of failed downloads.
    """
    downloaded_accessions = genome_complete_downloads.intersection(protein_complete_downloads)

    for index, row in bacteria_info.iterrows():
        assembly_accession = row['Assembly Accession']
        assembly_name = row['Assembly Name']

        # Skip already downloaded GCA and GCF with the same ID
        base_accession = assembly_accession.split('_')[1]
        if base_accession in downloaded_accessions:
            continue

        # Download and check genome file
        if download_and_check_file(assembly_accession, assembly_name, 'genomic', genome_dir, log_file, genome_complete_downloads, failed_downloads):
            downloaded_accessions.add(base_accession)

        # Download and check protein file
        download_and_check_file(assembly_accession, assembly_name, 'protein', protein_dir, log_file, protein_complete_downloads, failed_downloads)


def decompress_and_create_db(directory, extension, db_type):
    """Decompress files and create BLAST databases.

    Args:
        directory: Directory containing the files.
        extension: File extension to look for.
        db_type: Type of the database ('nucl' for nucleotide, 'prot' for protein).
    """
    for filename in os.listdir(directory):
        if filename.endswith(extension):
            gz_file = os.path.join(directory, filename)
            decompressed_file = decompress_gz_file(gz_file, directory)
            create_blast_db(decompressed_file, db_type=db_type)


def main():
    tsv_file = './bacteria20.tsv'
    output_dir = os.getcwd()
    genome_dir = os.path.join(output_dir, 'bacteria_genome')
    protein_dir = os.path.join(output_dir, 'bacteria_protein')
    genome_complete_file = 'genome_complete_downloads.txt'
    protein_complete_file = 'protein_complete_downloads.txt'
    faalist_file = 'faalist.txt'
    failed_file = 'failed_downloads.txt'
    log_file = 'wget.log'

    # Create directories to save genome and protein files
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)
    if not os.path.exists(protein_dir):
        os.makedirs(protein_dir)

    bacteria_info = read_bacteria_info(tsv_file)

    # Read the lists of completed and incomplete downloads
    if os.path.exists(genome_complete_file):
        with open(genome_complete_file, 'r') as f:
            genome_complete_downloads = set(f.read().splitlines())
    else:
        genome_complete_downloads = set()

    if os.path.exists(protein_complete_file):
        with open(protein_complete_file, 'r') as f:
            protein_complete_downloads = set(f.read().splitlines())
    else:
        protein_complete_downloads = set()

    if os.path.exists(faalist_file):
        with open(faalist_file, 'r') as f:
            faa_list = set(f.read().splitlines())
    else:
        faa_list = set()

    if os.path.exists(failed_file):
        with open(failed_file, 'r') as f:
            failed_downloads = set(f.read().splitlines())
    else:
        failed_downloads = set()

    # Scan existing genome and protein files
    genome_complete = scan_existing_files(genome_dir, '_genomic.fna.gz')
    protein_complete = scan_existing_files(protein_dir, '_protein.faa.gz')

    genome_complete_downloads.update(genome_complete)
    protein_complete_downloads.update(protein_complete)

    # Download data
    #download_data(bacteria_info, genome_dir, protein_dir, log_file, genome_complete_downloads, protein_complete_downloads, failed_downloads)

    # Save the lists of completed and incomplete downloads
    #update_download_lists(genome_complete_file, protein_complete_file, genome_complete_downloads, protein_complete_downloads)

    # Save the list of Assembly Accessions with faa.gz files
    #update_faa_list(faalist_file, protein_complete_downloads)

    # Save the list of failed downloads
    #update_failed_list(failed_file, failed_downloads)

    # Decompress and create BLAST databases for genome files
    decompress_and_create_db(genome_dir, '_genomic.fna.gz', db_type='nucl')

    # Decompress and create BLAST databases for protein files
    decompress_and_create_db(protein_dir, '_protein.faa.gz', db_type='prot')


if __name__ == '__main__':
    main()
    print("OK\n")