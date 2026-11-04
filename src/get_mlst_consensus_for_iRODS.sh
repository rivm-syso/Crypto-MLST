input_dir=$1
output_dir=$2

reference="src/fasta_seq_mlst.fasta"

for barcode in $input_dir/*; do
    if [ -d "$barcode" ]; then
        echo "Processing directory: $barcode"

        # Combine all FASTQ files into a single FASTQ file
        cat $barcode/*.fastq.gz > $output_dir/$(basename $barcode)_raw.fastq.gz

        # Map reads to the reference genes and sort the BAM file
        minimap2 -ax map-ont $reference $output_dir/$(basename $barcode)_raw.fastq.gz | \
        samtools view -u | \
        samtools sort -o $output_dir/$(basename $barcode)_sorted.bam && \

        # Check for mapped reads in the BAM file
        if [ "$(samtools view -c -F 4 $output_dir/$(basename $barcode)_sorted.bam)" -gt 0 ]; then

            # Create first consensus 
            bcftools mpileup -Ou --fasta-ref $reference $output_dir/$(basename $barcode)_sorted.bam | \
            bcftools call -c -Ob -o $output_dir/$(basename $barcode)_call.bcf.gz && \
            bcftools index $output_dir/$(basename $barcode)_call.bcf.gz && \
            bcftools consensus -f $reference -o $output_dir/$(basename $barcode)_consensus.fasta $output_dir/$(basename $barcode)_call.bcf.gz

            # Remap reads to the consensus sequence
            minimap2 -ax map-ont $output_dir/$(basename $barcode)_consensus.fasta $output_dir/$(basename $barcode)_raw.fastq.gz \
                | samtools view -u - \
                | samtools sort -o $output_dir/$(basename $barcode)_consensus_mapped.bam
            samtools index $output_dir/$(basename $barcode)_consensus_mapped.bam

            # Generate pileup on consensus
            samtools mpileup -aa --reference $output_dir/$(basename $barcode)_consensus.fasta \
                $output_dir/$(basename $barcode)_consensus_mapped.bam \
                > $output_dir/$(basename $barcode)_consensus_mapped.pileup

            # Parse mpileup and create base support file
            python3 src/parse_mpileup.py \
            "$output_dir/$(basename $barcode)_consensus_mapped.pileup" \
            -o "$output_dir/$(basename $barcode)_consensus_base_support.txt"

# Use base support file to create a new consensus, this time omitting regions with 0 fraction
python3 <<EOF

import pandas as pd

base_support = pd.read_csv("$output_dir/$(basename $barcode)_consensus_base_support.txt", sep="\t")

file_out = open("$output_dir/$(basename $barcode)_consensus_final.fasta", "w")
for gene in base_support["chrom"].unique():
    print(f">{gene}", file=file_out)
    print("".join(list(base_support.loc[(base_support["chrom"] == gene) & 
                       # Remove all bases that are not supported by reads
                       ~((base_support["fraction"] < 0.01) & (base_support["coverage"] > 0)), 
                       # Remove the first and last 10 bases
                       "cons_base"]))[10:-10], file=file_out)

file_out.close()
                                    # Show all bases with < 80% fraction (but not deletions, faction == 0) 
dubious_variants = base_support.loc[(((base_support["fraction"] < 0.8) & (base_support["fraction"] >= 0.01)) | 
                                    # Show all bases with coverage < 20 (but not first and last 10 bp of each gene)
                                        ((base_support["coverage"] < 20) & (base_support["coverage"] != 0))) &
                                    ~((base_support["pos"] < 11) | 
                                        (((base_support["chrom"] == "CP003930") & (base_support["pos"] > 441)) |
                                        ((base_support["chrom"] == "CP0028230") & (base_support["pos"] > 466)) |
                                        ((base_support["chrom"] == "CP0031960") & (base_support["pos"] > 577)) |
                                        ((base_support["chrom"] == "CP0021750") & (base_support["pos"] > 497)) |
                                        ((base_support["chrom"] == "CP0024650") & (base_support["pos"] > 516)) |
                                        ((base_support["chrom"] == "CP0012400") & (base_support["pos"] > 493)) | 
                                        ((base_support["chrom"] == "CP007000") & (base_support["pos"] > 548)) |
                                        ((base_support["chrom"] == "CP001010") & (base_support["pos"] > 543))))]
dubious_variants.to_csv("$output_dir/$(basename $barcode)_consensus_final_dubious_bases.txt", sep="\t", index=None)

EOF

            # Create a depth plot
            samtools depth -a $output_dir/$(basename $barcode)_sorted.bam > $output_dir/$(basename $barcode)_depth.txt
            python src/create_depth_plot.py \
                $output_dir/$(basename $barcode)_depth.txt \
                $output_dir/$(basename $barcode)_depth_plot.pdf
            
            rm $output_dir/$(basename $barcode)_depth.txt
            rm $output_dir/$(basename $barcode)_call.bcf.gz
            rm $output_dir/$(basename $barcode)_call.bcf.gz.csi
            rm $output_dir/$(basename $barcode)*.bam
            rm $output_dir/$(basename $barcode)*.bam.bai
            rm $output_dir/$(basename $barcode)*.pileup
            rm $output_dir/$(basename $barcode)_raw.fastq.gz
            rm $output_dir/$(basename $barcode)_consensus.fasta
            rm $output_dir/$(basename $barcode)_consensus.fasta.fai
        
        else
            rm -f "$output_dir/$(basename $barcode)_raw.fastq.gz"
            rm -f "$output_dir/$(basename $barcode)_sorted.bam"
            rm -f "$output_dir/$(basename $barcode)_sorted.bam.bai"
            touch "$output_dir/no_reads_mapped_for_barcode_$(basename $barcode).txt"

        fi
    fi
done
