#Author: Joshua Topper

import os
import gzip
import subprocess
import pysam

class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self, filePath, headerSymbols=['@', '+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath, 'rt')  # 'rt' for text mode with gzip
        else:
            self._file = open(filePath, 'r')  # 'r' for text mode
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fastQ data as tuple ++++
        return tuple(elemList)
    
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._file.close()

####
#Creates a dictionary from the provided clinical data
def read_clinical_data(clinical_data_file):
    """Reads the clinical data from the specified file and returns it as a dictionary."""
    clinical_data = {}
    with open(clinical_data_file, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            fields = line.strip().split('\t')
            sample_name = fields[0]
            barcode = fields[2]  # Assuming barcode is in the third column
            color = fields[1]  # Assuming color is in the second column
            clinical_data[sample_name] = {"Barcode": barcode, "Color": color, "Name": sample_name}
    return clinical_data

#Function to trim the beginning of the reads
def trim_beginning(read_sequence, quality_scores):
    #Assuming the barcode length is known (e.g., 5)
    barcode_length = 5
    return read_sequence[barcode_length:], quality_scores[barcode_length:]

#Function to trim the end of the reads. Determines what to trim based on
#continous poor quality reads. 
def trim_ends(read_sequence, quality_scores):
    last_index = len(quality_scores) - 1
    while last_index >= 0 and quality_scores[last_index] in ['D', 'F']:
        last_index -= 1
    
    return read_sequence[:last_index + 1], quality_scores[:last_index + 1]


#Creates the new sorted fastq files and trims them as specified.
def demultiplex_fastq(input_fastq, clinical_data, output_folder):
    #Demultiplexes the fastq file based on the barcodes in the clinical data.
    os.makedirs(output_folder, exist_ok=True)  # Create the output folder if it doesn't exist
    
    # Iterate over each entry in the clinical data
    for sample in clinical_data:
        barcode = clinical_data[sample]["Barcode"]
        output_filename = os.path.join(output_folder, f"{sample.lower()}_trimmed.fastq")
        output_file = open(output_filename, "w")
        print(f"Processing {sample} ({barcode})...")
        
        with ParseFastQ(input_fastq) as fastq_parser:
            for read in fastq_parser:
                sequence = read[1]  # Extract the sequence from the read tuple
                if sequence.startswith(barcode):
                    #print(f"Found matching barcode for {sample}. Writing read to {output_filename}")
                    # Trim the beginning and end of the sequence and quality scores
                    trimmed_sequence, trimmed_quality = trim_beginning(read[1], read[3])
                    trimmed_sequence, trimmed_quality = trim_ends(trimmed_sequence, trimmed_quality)
                    # Write the trimmed read to the output file
                    output_file.write(f"{read[0]}\n{trimmed_sequence}\n{read[2]}\n{trimmed_quality}\n")
        
        output_file.close()


#Call demultpilex function
if __name__ == "__main__":
    # Input file paths
    input_fastq = "hawkins_pooled_sequences.fastq"
    clinical_data_file = "harrington_clinical_data.txt"
    #Output folder path
    output_folder = "fastqs"

    #Read clinical data
    barcodes = read_clinical_data(clinical_data_file)
    print("Clinical data (first few entries):")
    for i, (barcode, data) in enumerate(barcodes.items()):
        print(f"{i+1}. Barcode: {barcode}, Name: {data['Name']}, Color: {data['Color']}")
        if i == 4:  #Print only the first 5 entries
            break

    #Demultiplex and trim fastq file
    demultiplex_fastq(input_fastq, barcodes, output_folder)

####
#Call the bwa command to transform fastq files into samfiles and align them.
def perform_alignment(fastq_folder, output_folder):
    #Index the reference file if not already indexed
    reference_file = "dgorgon_reference.fa"
    if not os.path.exists(reference_file + ".bwt"):
        subprocess.run(["bwa", "index", reference_file])
    
    #Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    #Iterate over each FASTQ file in the fastq_folder
    for fastq_file in os.listdir(fastq_folder):
        if fastq_file.endswith(".fastq"):
            #Generate output SAM filename
            fastq_basename = os.path.splitext(fastq_file)[0]
            output_sam = os.path.join(output_folder, f"{fastq_basename.replace('_trimmed', '')}.sam")
            
            #Perform alignment using BWA mem
            with open(output_sam, "w") as out_sam_file:
                subprocess.run(["bwa", "mem", reference_file, os.path.join(fastq_folder, fastq_file)], stdout=out_sam_file)

#Output to alignments
fastq_folder = "fastqs"
output_folder = "bams"

perform_alignment(fastq_folder, output_folder)

####
#The following functions are used to convert the sam files to bam files,
#index them, and then remove the intermediate files. This will leave just
#the sorted and indexed bam files.

def convert_sam_to_bam(input_sam, output_bam):
    #Convert SAM to BAM using samtools view.
    subprocess.run(["samtools", "view", "-bS", "-o", output_bam, input_sam])

def sort_bam(input_bam, output_sorted_bam, memory_limit="100M"):
    #Sort BAM files using samtools sort.
    subprocess.run(["samtools", "sort", "-m", memory_limit, "-o", output_sorted_bam, input_bam])

def index_bam(input_sorted_bam):
    #Index sorted BAM files using samtools index.
    subprocess.run(["samtools", "index", input_sorted_bam])

def convert_and_sort_bams(bam_folder):
    #Convert SAM to BAM, sort, index, and delete intermediate files.
    for sam_file in os.listdir(bam_folder):
        if sam_file.endswith(".sam"):
            #Convert SAM to BAM
            input_sam = os.path.join(bam_folder, sam_file)
            output_bam = os.path.join(bam_folder, sam_file.replace(".sam", ".bam"))
            convert_sam_to_bam(input_sam, output_bam)

            #Sort BAM
            output_sorted_bam = os.path.join(bam_folder, sam_file.replace(".sam", ".sorted.bam"))
            sort_bam(output_bam, output_sorted_bam)

            #Index sorted BAM
            index_bam(output_sorted_bam)

            #Delete intermediate files
            os.remove(input_sam)
            os.remove(output_bam)

#Output to bams
bam_folder = "bams"
convert_and_sort_bams(bam_folder)

####
#Performs Variant discovery for each individual sample.
def pileup(clinical_data, bam_folder, reference_file):
    mutation_report = {}
    reference_sequence = pysam.FastaFile(reference_file).fetch("Dgorgon")

    for sample, data in clinical_data.items():
        bam_file = os.path.join(bam_folder, f"{sample.lower()}.sorted.bam")
        samfile = pysam.AlignmentFile(bam_file, "rb")
        mutation_counts = {}

        #Count total reads from the FASTQ file
        total_reads = sum(1 for _ in samfile.fetch("Dgorgon"))

        #Build ntdict for bases.
        for pileupcolumn in samfile.pileup(contig="Dgorgon"):
            position = pileupcolumn.pos
            ref_base = reference_sequence[position].upper()
            ntdict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    continue

                base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                ntdict[base] += 1

            #Check for mutations
            for base, count in ntdict.items():
                if base != ref_base and count > 0:
                    mutation_counts.setdefault(position, {}).setdefault(base, 0)
                    mutation_counts[position][base] += count
        #Count total reads.
        if total_reads > 0:
            mutation_report[sample] = {
                "TotalReads": total_reads,
                "Mutations": [{"Position": pos, "Mutation": mutation, "Count": count}
                              for pos, mutations in mutation_counts.items() for mutation, count in mutations.items()]
            }

        samfile.close()

    return mutation_report

#Compile necessary information and write to the final report.txt file.
def write_report(mutation_report, clinical_data, output_file):
    with open(output_file, "w") as report_file:
        for sample, data in mutation_report.items():
            sample_data = clinical_data.get(sample)
            if sample_data is None:
                continue
            
            color = sample_data.get("Color", "Unknown").lower()
            total_reads = data.get("TotalReads", 0)
            
            report_file.write(f"Sample {sample} had a(n) {color} mold, {total_reads} reads, and\n")
            
            if "Mutations" not in data or len(data["Mutations"]) == 0:
                report_file.write("no mutations detected.\n")
            else:
                report_file.write("had") # <-- I thought this was a bit silly on my part
                for mutation_info in data["Mutations"]:
                    position = mutation_info['Position']
                    mutation = mutation_info['Mutation']
                    mutation_count = mutation_info.get('Count', 0)
                    percentage = (mutation_count / total_reads) * 100 if total_reads > 0 else 0
                    report_file.write(f" {percentage:.2f}% of the reads at position {position} had the mutation {mutation}.\n")
            report_file.write("\n")



clinical_data_file = "harrington_clinical_data.txt"  #Adjust file path as needed
bam_folder = "bams"  #Adjust folder name as needed
reference_file = "dgorgon_reference.fa"  #Adjust file name as needed

clinical_data = read_clinical_data(clinical_data_file)
mutation_report = pileup(clinical_data, bam_folder, reference_file)  #Assign the result of pileup function
output_file = "report.txt"  #Adjust file name as needed

write_report(mutation_report,clinical_data, output_file)  #Call write_report with mutation_report and output_file

