# Go to the folder containing your FASTQ files
cd ~/Desktop/trim

# Create an output folder so original files are not overwritten
mkdir -p trimmed_output

# Loop over all R1 files and automatically pair R2
for R1 in *R1_001.fastq*  # matches both .fastq and .fastq.gz
do
    # Automatically find the corresponding R2
    R2="${R1/R1_001/R2_001}"

    # Define clean output filenames in the output folder
    out1="trimmed_output/${R1%.fastq*}_trimmed.fastq.gz"
    out2="trimmed_output/${R2%.fastq*}_trimmed.fastq.gz"

    # Run cutadapt to trim poly-C and poly-G at ends
    python3 -m cutadapt -g "C{20}" -g "G{20}" -a "C{20}" -a "G{20}" \
        -o "$out1" -p "$out2" \
        "$R1" "$R2"
done

Put into base terminal this is to trim the samples do this first