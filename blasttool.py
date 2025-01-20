import os
import pandas as pd
import subprocess
import gzip


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


def main():
    tsv_file = './bacteria20.tsv'
    output_dir = os.getcwd()
    genome_dir = os.path.join(output_dir, 'bacteria_genome')
    protein_dir = os.path.join(output_dir, 'bacteria_protein')
    complete_file = 'complete_downloads.txt'
    incomplete_file = 'incomplete_downloads.txt'
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
    if os.path.exists(complete_file):
        with open(complete_file, 'r') as f:
            complete_downloads = set(f.read().splitlines())
    else:
        complete_downloads = set()

    if os.path.exists(incomplete_file):
        with open(incomplete_file, 'r') as f:
            incomplete_downloads = set(f.read().splitlines())
    else:
        incomplete_downloads = set()

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
    complete_downloads.update(scan_existing_files(genome_dir, '_genomic.fna.gz'))
    faa_list.update(scan_existing_files(protein_dir, '_protein.faa.gz'))

    downloaded_accessions = set()
    downloaded_accessions = complete_downloads.intersection(faa_list)

    for index, row in bacteria_info.iterrows():
        assembly_accession = row['Assembly Accession']
        assembly_name = row['Assembly Name']

        # Skip already downloaded GCA and GCF with the same ID
        base_accession = assembly_accession.split('_')[1]
        if base_accession in downloaded_accessions:
            continue

        # Download and check genome file
        if download_and_check_file(assembly_accession, assembly_name, 'genomic', genome_dir, log_file, complete_downloads, failed_downloads):
            downloaded_accessions.add(base_accession)

        # Download and check protein file
        download_and_check_file(assembly_accession, assembly_name, 'protein', protein_dir, log_file, faa_list, failed_downloads)

    # Save the lists of completed and incomplete downloads
    update_download_lists(complete_file, incomplete_file, complete_downloads, incomplete_downloads)

    # Save the list of Assembly Accessions with faa.gz files
    update_faa_list(faalist_file, faa_list)

    # Save the list of failed downloads
    update_failed_list(failed_file, failed_downloads)


if __name__ == '__main__':
    main()
    print("OK\n")