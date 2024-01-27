from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import os

directory_path = "C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_08_16_Calm_mem_2/"  # Replace with your directory path
calm_db="C:/blast/BLAST_DBs/calmodulin/calm_db"

output_path=f'{directory_path}/Calmodulin' #Where should the results go
fastqpath = f'{directory_path}/pass'

# Convert FastQ to Fasta
def fastq_to_fasta(fastq_files, output_fasta):
    with open(output_fasta, "w") as output_handle:
        for fastq_file in fastq_files:
            sequences = SeqIO.parse(fastq_file, "fastq")
            SeqIO.write(sequences, output_handle, "fasta")

# BLAST the sequences
def run_blast(query_path):
    blastn = NcbiblastnCommandline(cmd="C:/blast/bin/blastn.exe", query=query_path, db=calm_db, outfmt=6)
    stdout, stderr = blastn()
    return stdout

def main():
    
    fastq_files = [os.path.join(fastqpath, file) for file in os.listdir(fastqpath) if file.endswith(".fastq")]
    output_fasta = "output.fasta"
    print("fastqfiles_found")
    fastq_to_fasta(fastq_files, output_fasta)
    print("fasta_files created")
    blast_results = run_blast(output_fasta)
    print("blast_done")
    # Convert BLAST results to CSV format
    csv_header = "read_id,alignment_genome,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"
    csv_content = blast_results.replace('\t', ',')
    
    # Save BLAST results to a file
    with open(f"{output_path}/blast_results.csv", "w") as result_file:
        result_file.write(csv_header + '\n' + csv_content)

    print("BLAST results saved to 'blast_results.csv'")

        
if __name__ == "__main__":
    main()


