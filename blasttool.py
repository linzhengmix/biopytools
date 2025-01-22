import os
import pandas as pd
import subprocess
import gzip
import shutil
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
import logging
from logging.handlers import QueueHandler, QueueListener
import queue

# Configure logging
log_queue = queue.Queue()
queue_handler = QueueHandler(log_queue)
formatter = logging.Formatter('%(asctime)s %(message)s')
handler = logging.FileHandler('wget.log')
handler.setFormatter(formatter)
listener = QueueListener(log_queue, handler)
listener.start()

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(queue_handler)


def read_bacteria_info(tsv_file):
    """Read bacterial genome information from a TSV file.

    Args:
        tsv_file: Path to the TSV file.

    Returns:
        DataFrame containing genome information.
    """
    return pd.read_csv(tsv_file, sep='\t')


def download_genome_wget(ftp_url, output_dir, filename):
    """Download bacterial genome file using wget.

    Args:
        ftp_url: FTP URL.
        output_dir: Directory to save the downloaded file.
        filename: Name of the downloaded file.

    Returns:
        Path to the downloaded file.
    """
    local_filepath = os.path.join(output_dir, filename)
    cmd = ['wget', '-O', local_filepath, ftp_url]
    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.info(result.stdout)
    if result.stderr:
        logger.error(result.stderr)
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
        logger.error(f"File integrity check failed for {filepath}: {e}")
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


def download_and_check_file(assembly_accession, assembly_name, file_type, output_dir, download_list, failed_downloads, lock):
    """Download and check the integrity of a file.

    Args:
        assembly_accession: Assembly Accession of the genome.
        assembly_name: Assembly Name of the genome.
        file_type: File type (genomic or protein).
        output_dir: Directory to save the file.
        download_list: Set of downloaded files.
        failed_downloads: Set of failed downloads.
        lock: Threading lock for thread safety.

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
                logger.info(f"Downloading {filename} from {ftp_url} (Attempt {attempt + 1})")
                download_genome_wget(ftp_url, output_dir, filename)

                if check_file_integrity(file_path):
                    logger.info(f"Download and integrity check passed for {filename}")
                    with lock:
                        download_list.add(assembly_accession)
                    return True
                else:
                    logger.error(f"Download or integrity check failed for {filename}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to download {filename} (Attempt {attempt + 1}): {e}")
        else:
            with lock:
                failed_downloads.add(assembly_accession)
            return False
    return True


def decompress_gz_file(gz_file, output_dir):
    """Decompress a .gz file into a separate directory.

    Args:
        gz_file: Path to the .gz file.
        output_dir: Directory to save the decompressed file.

    Returns:
        Path to the decompressed file.
    """
    base_name = os.path.basename(gz_file).replace('.gz', '')
    decompressed_dir = os.path.join(output_dir, base_name)
    os.makedirs(decompressed_dir, exist_ok=True)
    decompressed_file = os.path.join(decompressed_dir, base_name)
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
    logger.info(f"Creating BLAST database for {os.path.basename(fasta_file)}")
    db_files = [f"{db_name}.fna.{ext}" for ext in ('nhr', 'nin', 'nsq', 'ndb', 'not', 'ntf', 'nto')] if db_type == 'nucl' else [f"{db_name}.faa.{ext}" for ext in ('phr', 'pin', 'psq', 'pdb', 'pot', 'ptf', 'pto')]
    if not all(os.path.exists(db_file) for db_file in db_files):
        cmd = ['makeblastdb', '-in', fasta_file, '-dbtype', db_type]
        subprocess.run(cmd, check=True)
    else:
        logger.info(f"BLAST database for {os.path.basename(fasta_file)} already exists. Skipping creation.")


def download_data(bacteria_info, genome_dir, protein_dir, genome_complete_downloads, protein_complete_downloads, failed_downloads, lock):
    """Download genome and protein data.

    Args:
        bacteria_info: DataFrame containing bacterial genome information.
        genome_dir: Directory to save genome files.
        protein_dir: Directory to save protein files.
        genome_complete_downloads: Set of completed genome downloads.
        protein_complete_downloads: Set of completed protein downloads.
        failed_downloads: Set of failed downloads.
        lock: Threading lock for thread safety.
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
        if download_and_check_file(assembly_accession, assembly_name, 'genomic', genome_dir, genome_complete_downloads, failed_downloads, lock):
            downloaded_accessions.add(base_accession)

        # Download and check protein file
        download_and_check_file(assembly_accession, assembly_name, 'protein', protein_dir, protein_complete_downloads, failed_downloads, lock)


def decompress_and_create_db_for_file(gz_file, output_dir, db_type):
    """Decompress a file and create a BLAST database.

    Args:
        gz_file: Path to the .gz file.
        output_dir: Directory to save the decompressed file.
        db_type: Type of the database ('nucl' for nucleotide, 'prot' for protein).

    Returns:
        None
    """
    decompressed_file = decompress_gz_file(gz_file, output_dir)
    create_blast_db(decompressed_file, db_type=db_type)


def decompress_and_create_db(directory, extension, db_type):
    """Decompress files and create BLAST databases.

    Args:
        directory: Directory containing the files.
        extension: File extension to look for.
        db_type: Type of the database ('nucl' for nucleotide, 'prot' for protein).
    """
    max_threads = os.cpu_count() // 3
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = []
        for filename in os.listdir(directory):
            if filename.endswith(extension):
                gz_file = os.path.join(directory, filename)
                futures.append(executor.submit(decompress_and_create_db_for_file, gz_file, directory, db_type))
        for future in futures:
            future.result()


def run_blast_for_genome(query_file, db_subdir, db_dir, output_dir, blast_type):
    """Run BLAST for a specific genome.

    Args:
        query_file: Path to the query file.
        db_subdir: Subdirectory of the BLAST database.
        db_dir: Directory containing the BLAST databases.
        output_dir: Directory to save the output files.
        blast_type: Type of BLAST to run ('tblastn' or 'blastp').

    Returns:
        None
    """
    db_name = os.path.join(db_dir, db_subdir, db_subdir)
    output_file = os.path.join(output_dir, f"{db_subdir}_{blast_type}_results.txt")

    # Check if the output file already has the completion marker
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            if f"#{blast_type} run completed" in f.read():
                logger.info(f"Skipping {db_subdir} as it is already completed.")
                return

    logger.info(f"Running {blast_type} for {db_subdir}...")
    cmd = [
        blast_type,
        '-query', query_file,
        '-db', db_name,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid qstart qend sstart send evalue qcovhsp qcovs sframe sseq qseq pident'
    ]
    subprocess.run(cmd, check=True)
    with open(output_file, 'a') as f:
        f.write(f"\n# {blast_type} run completed\n")


def run_tblastn(query_file, db_dir, output_dir):
    """Run tblastn using the sequences from the query file.

    Args:
        query_file: Path to the query file.
        db_dir: Directory containing the BLAST databases.
        output_dir: Directory to save the output files.

    Returns:
        None
    """
    db_dirs = [d for d in os.listdir(db_dir) if os.path.isdir(os.path.join(db_dir, d))]
    threads = []
    for db_subdir in db_dirs:
        thread = threading.Thread(target=run_blast_for_genome, args=(query_file, db_subdir, db_dir, output_dir, 'tblastn'))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()


def run_blastp(query_file, db_dir, output_dir):
    """Run blastp using the sequences from the query file.

    Args:
        query_file: Path to the query file.
        db_dir: Directory containing the BLAST databases.
        output_dir: Directory to save the output files.

    Returns:
        None
    """
    db_dirs = [d for d in os.listdir(db_dir) if os.path.isdir(os.path.join(db_dir, d))]
    threads = []
    for db_subdir in db_dirs:
        thread = threading.Thread(target=run_blast_for_genome, args=(query_file, db_subdir, db_dir, output_dir, 'blastp'))
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()


def main():
    tsv_file = './bacteria20.tsv'
    output_dir = os.getcwd()
    genome_dir = os.path.join(output_dir, 'bacteria_genome')
    protein_dir = os.path.join(output_dir, 'bacteria_protein')
    genome_complete_file = 'genome_complete_downloads.txt'
    protein_complete_file = 'protein_complete_downloads.txt'
    failed_file = 'failed_downloads.txt'
    query_file = './copper_protein_as_tblastn_query.txt'
    tblastn_output_dir = os.path.join(output_dir, 'tblastn_results')
    blastp_output_dir = os.path.join(output_dir, 'blastp_results')

    # Create directories to save genome and protein files
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)
    if not os.path.exists(protein_dir):
        os.makedirs(protein_dir)
    if not os.path.exists(tblastn_output_dir):
        os.makedirs(tblastn_output_dir)
    if not os.path.exists(blastp_output_dir):
        os.makedirs(blastp_output_dir)

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

    lock = threading.Lock()

    # Download data
    download_data(bacteria_info, genome_dir, protein_dir, genome_complete_downloads, protein_complete_downloads, failed_downloads, lock)

    # Save the lists of completed and incomplete downloads
    update_download_lists(genome_complete_file, protein_complete_file, genome_complete_downloads, protein_complete_downloads)

    # Save the list of failed downloads
    update_failed_list(failed_file, failed_downloads)

    # Decompress and create BLAST databases for genome files
    decompress_and_create_db(genome_dir, '_genomic.fna.gz', db_type='nucl')

    # Decompress and create BLAST databases for protein files
    decompress_and_create_db(protein_dir, '_protein.faa.gz', db_type='prot')

    # Run tblastn and blastp using threading
    tblastn_thread = threading.Thread(target=run_tblastn, args=(query_file, genome_dir, tblastn_output_dir))
    blastp_thread = threading.Thread(target=run_blastp, args=(query_file, protein_dir, blastp_output_dir))

    tblastn_thread.start()
    blastp_thread.start()

    tblastn_thread.join()
    blastp_thread.join()

if __name__ == '__main__':
    main()
    print("\nOK\n")