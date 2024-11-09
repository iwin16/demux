import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def demultiplex_paired_fastq_by_barcode(fastq_r1, fastq_r2, barcode_file, output_dir, log_file):
    # Load barcodes and compile regex patterns
    barcodes = {}
    with open(barcode_file, 'r') as bf:
        for line in bf:
            barcode, sample_name = line.strip("\n").split("\t")
            barcodes[barcode] = sample_name
    barcode_patterns = {barcode: re.compile(barcode) for barcode in barcodes.keys()}

    # Prepare output files and log
    out_files_r1 = {sample: open(f"{output_dir}/{sample}_R1.fastq", 'w') for sample in barcodes.values()}
    out_files_r2 = {sample: open(f"{output_dir}/{sample}_R2.fastq", 'w') for sample in barcodes.values()}
    unmatched_log = open(log_file, 'w')

    # Track written sequences to avoid duplicates
    written_ids = set()

    # Preload records for R2 and R1 for faster access
    r1_records = {rec.id.split()[0]: rec for rec in SeqIO.parse(fastq_r1, "fastq")}
    r2_records = {rec.id.split()[0]: rec for rec in SeqIO.parse(fastq_r2, "fastq")}

    def process_fastq(fastq_handle, read_type):
        unmatched_ids = set()

        for record in SeqIO.parse(fastq_handle, "fastq"):
            record_id = record.id.split()[0]
            sequence = str(record.seq)
            quality = record.letter_annotations["phred_quality"]

            # Try to match each barcode pattern at any position within the sequence
            matched = False
            for barcode, sample_name in barcodes.items():
                pattern = barcode_patterns[barcode]
                match = pattern.search(sequence)

                if match:
                    print(f"Match found for record {record_id} with barcode {barcode} and sample {sample_name}")

                    # Remove barcode from the sequence and quality scores
                    start, end = match.span()
                    trimmed_sequence = sequence[:start] + sequence[end:]
                    trimmed_quality = quality[:start] + quality[end:]

                    # Retrieve the paired record based on read type
                    paired_record = r2_records[record_id] if read_type == "R1" else r1_records[record_id]

                    # Process the paired read sequence
                    paired_sequence = str(paired_record.seq)
                    paired_quality = paired_record.letter_annotations["phred_quality"]
                    paired_match = pattern.search(paired_sequence)
                    if paired_match:
                        p_start, p_end = paired_match.span()
                        paired_sequence = paired_sequence[:p_start] + paired_sequence[p_end:]
                        paired_quality = paired_quality[:p_start] + paired_quality[p_end:]

                    # Write once to prevent duplicates
                    if record_id not in written_ids:
                        # Clear annotations before assigning new sequence
                        record.annotations.clear()
                        record.letter_annotations.clear()
                        paired_record.annotations.clear()
                        paired_record.letter_annotations.clear()

                        # Assign new trimmed sequences
                        record.seq = Seq(trimmed_sequence)
                        paired_record.seq = Seq(paired_sequence)

                        # Update quality annotations
                        record.letter_annotations["phred_quality"] = trimmed_quality
                        paired_record.letter_annotations["phred_quality"] = paired_quality

                        # Write to output files
                        SeqIO.write(record, out_files_r1[sample_name], "fastq")
                        SeqIO.write(paired_record, out_files_r2[sample_name], "fastq")
                        written_ids.add(record_id)
                    matched = True
                    break

            if not matched:
                print(f"No match for record {record_id}")
                unmatched_ids.add(record_id)

        return unmatched_ids

    # Process R1 and R2 files independently
    with open(fastq_r1, 'r') as r1_handle:
        unmatched_r1 = process_fastq(r1_handle, "R1")
    with open(fastq_r2, 'r') as r2_handle:
        unmatched_r2 = process_fastq(r2_handle, "R2")

    # Log unmatched reads
    unmatched_log.write("Unmatched R1:\n" + "\n".join(unmatched_r1) + "\n")
    unmatched_log.write("Unmatched R2:\n" + "\n".join(unmatched_r2) + "\n")

    # Close all files
    for f in out_files_r1.values():
        f.close()
    for f in out_files_r2.values():
        f.close()
    unmatched_log.close()

# Example usage:
demultiplex_paired_fastq_by_barcode(
    "input/path/KJ3_R1.fastq", 
    "input/path/KJ3_R2.fastq", 
    "input/path/KJ3_barcode.txt", 
    "./kj3", 
    "kj3/log_file.txt"
)
