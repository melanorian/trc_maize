import re
import csv

# Input and output file names
log_file = "/home/melanie/net/virus_melanie/maize/results/t4_all_v1/salmon_quant.log"
output_csv = "/home/melanie/net/virus_melanie/maize/results/t4_all_v1/salmon_summary.csv"

# Define regex patterns for extraction
patterns = {
    "sample": re.compile(r"Processing\s+/home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R1/([\w\d]+)_\d\.fq\.gz"),
    "equivalence_classes": re.compile(r"Computed\s+([\d,]+)\s+rich equivalence classes"),
    "Reads": re.compile(r"Counted\s+([\d,]+)\s+total reads"),
    "nr_discarded_alignment": re.compile(r"Number of mappings discarded because of alignment score\s*:\s*([\d,]+)"),
    "nr_discarded_entirely_alignment": re.compile(r"Number of fragments entirely discarded because of alignment score\s*:\s*([\d,]+)"),
    "nr_discarded_decoy": re.compile(r"Number of fragments discarded because they are best-mapped to decoys\s*:\s*([\d,]+)"),
    "nr_discarded_dovetail": re.compile(r"Number of fragments discarded because they have only dovetail \(discordant\) mappings to valid targets\s*:\s*([\d,]+)"),
    "Mapping_rate": re.compile(r"Mapping rate\s*=\s*([\d.]+)%")
}

# Initialize data storage
data = []
current_sample = {}

# Read and process log file
with open(log_file, "r") as file:
    for line in file:
        # Check for new sample
        match_sample = patterns["sample"].search(line)
        if match_sample:
            # If a previous sample exists, save it before starting a new one
            if current_sample:
                data.append(current_sample)
            # Start a new sample entry
            sample_name = match_sample.group(1)
            current_sample = {"sample": sample_name}

        # Extract relevant metrics
        for key, pattern in patterns.items():
            if key != "sample":
                match = pattern.search(line)
                if match:
                    current_sample[key] = match.group(1).replace(
                        ",", "")  # Remove commas for numeric values

# Append the last sample found
if current_sample:
    data.append(current_sample)

# Write data to CSV file
with open(output_csv, "w", newline="") as csvfile:
    fieldnames = ["sample", "equivalence_classes", "Reads", "nr_discarded_alignment",
                  "nr_discarded_entirely_alignment", "nr_discarded_decoy",
                  "nr_discarded_dovetail", "Mapping_rate"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(data)

print(f"CSV summary saved to {output_csv}")
