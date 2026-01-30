# Nanopore-workflow-user-guide
This is page is meant to serve as a detailed user guide for beginner-coders working with nanopore.

## Access to cluster
The workflow can be best run on the SILS computer cluster Crunchomics HPC to get results fast while not occupying your computer's own computational resources. You may request an account by emailing W.C.deLeeuw@uva.nl mentioning your UvAnetID and purpose of using the cluster. Now follow the steps provided by the Crunchomics documentation to get started https://crunchomics-documentation.readthedocs.io/en/latest. To help yourself understand basic navigation commands in bash please refer to the following cheat sheet:
<img width="1536" height="1086" alt="image" src="https://github.com/user-attachments/assets/800de4d7-ab1c-43b7-9eaa-d3974b44bfac" />


## About NaMeco
NaMeco is the pipeline we will use to analyze our data and has to be installed using conda. If you have followed the Crunchomics documentation and suggested setups you should have conda installed. For more information about NaMeco always refer to https://github.com/timyerg/NaMeco for the latest updates. If you have any issues try to ask an experienced colleague for help first, but otherwise the developer is very quick to respond on GitHub to any issues you may encounter.

## Setting up your data
First I would recommend to make a folder for your nanopore data. In MobaXTerm, this can be done interactively in the navigation panel on the left. Please note that this is for quick navigating and will not set your working directory. You can use 'cd nanopore' to move to your nanopore directory (refer to cheat sheet) and you can go 1 level back by using 'cd ..' . In your nanopore folder, you can now upload all your barcode folders (see icon in navigation pannel); the contents of which should contain multiple .fastq.gz files.

First we will make a new folder titled mergedreads and merge the fastq files of the same barcode and transfer them to the mergedreads folder using the following script:
```
#!/bin/bash

# 1. Create the destination folder
mkdir -p mergedreads

echo "Starting merge process for compressed files (barcode1 - barcode20)..."

# 2. Loop through numbers 1 to 20
for i in {1..20}; do
    SOURCE_DIR="barcode${i}"
    # Note the .gz extension in the output filename
    OUTPUT_FILE="mergedreads/barcode${i}merged.fastq.gz"

    # 3. Check if directory exists
    if [ -d "$SOURCE_DIR" ]; then
        
        # Check if .fastq.gz files exist inside
        if ls "$SOURCE_DIR"/*.fastq.gz >/dev/null 2>&1; then
            echo "Processing $SOURCE_DIR..."
            
            # Merge all .fastq.gz files into one .fastq.gz file
            cat "$SOURCE_DIR"/*.fastq.gz > "$OUTPUT_FILE"
            
            echo "  -> Created: $OUTPUT_FILE"
        else
            echo "  -> Warning: No .fastq.gz files found in $SOURCE_DIR"
        fi
        
    else
        echo "Directory '$SOURCE_DIR' not found. Skipping."
    fi
done

echo "-----------------------------------"
echo "Process complete."
```
You should now find 1 .fastq.gz file for every barcode in the merged reads folder, meaning we can now install and run NaMeco. Note that you tell the loop to check for barcode 1-20. You may adjust this number to however many barcodes you have for your run.

# Installing NaMeco
You can install NaMeco from the console using:
```
wget https://raw.githubusercontent.com/timyerg/NaMeco/main/NaMeco.yaml
conda env create --file NaMeco.yaml
conda activate NaMeco
```
This has also created and activated an environment that contains all the packages you need. You should see (NaMeco) now instead of (base) to the left of your username in the terminal.

# Running NaMeco
NaMeco is very  user friendly in the sense that it can do everything in 1 command. There are a lot of adjustable options that are run at default settings when not explicitely mentioned in your command line. You can find all settings, what they do and how to adjust them on the NaMeco page https://github.com/timyerg/NaMeco. There are a couple settings I recommend to adjust:
- Your input directory, this will be the path to the mergedreads folder (e.g. --inp_dir /zfs/omics/personal/15827127/nanopore/test2/mergedreads/)
- Your output directory, this can be a folder of  your choosing like 'results' (e.g. --out_dir /zfs/omics/personal/15827127/nanopore/test2/results2/)
- The amount of threads (computing power) the cluster will allocate to the analysis. More threads = faster but also a higher burden on the cluster. In my personal experience setting this to 30 threads falls easily within the socially accepted cluster burden while maintaining high speed analysis. You can play around with this setting depending on the size of your dataset and experience in the runtime.
- Forward primer setting. Default primers from NaMeco may be different from the primers we use in the lab so its helpful to specify the primers that were used
- Reverse primer setting idem

Below is the command I have used for my runs:
```
nameco --inp_dir /zfs/omics/personal/15827127/nanopore/test2/mergedreads/ --out_dir /zfs/omics/personal/15827127/nanopore/test2/results2/ --threads 30 --primer_F AGAGTTTGATCCTGGCTCAG --primer_R GGTTACCTTGTTACGACTT 
```
NaMeco will now run for a while. Close to the start it will give you some warnings about pkg_resources being depricated but this can be ignored.

# Data analysis
Congrats! You've done it!

Your results should be in your results folder, then final_output. I recommend to download the Speciescount.tsv for further analysis and figure making in R or python, or results can even be exported to QIIME2 for which I will recommend the NaMeco developer page (https://github.com/timyerg/NaMeco). Below is a quick example of a bar plot so you can easily check what your data looks like in R:

## Bar plot
First we load in the species count table and turn it into a table (dataframe) for R. Note that all of the following codes in this chapter are for R and not to be run in MobaXTerm!! Don't forget to change the file paths to where you have your Species_counts.tsv downloaded.
```
# Load necessary libraries
library(tidyverse)

# 1. LOAD AND PREPROCESS DATA ---------------------------------------------
# Read the data (assuming the file is in your working directory)
df <- read_tsv('C:/Users/rickh/OneDrive - UvA/202511_Rick Hoogendijk/NanoporeWorkflow/barcode19/Species_counts.tsv')

# Rename columns for easier handling if needed, though 'Species' and 'merged' work fine.
# We will calculate Relative Abundance (%) and Log Counts for better visualization.
df_processed <- df %>%
  mutate(
    Relative_Abundance = (merged / sum(merged)) * 100,
    Log_Counts = log10(merged + 1), # +1 to avoid log(0)
    # Reorder Species factor by count for sorted plots
    Species = fct_reorder(Species, merged) 
  )

# Print summary to console
print(paste("Total Reads:", sum(df$merged)))
print(paste("Number of Species Detected:", nrow(df)))

# 2. PROCESS DATA
# Calculate Relative Abundance (%) and order species by abundance
df_abundance <- df %>%
  mutate(
    Relative_Abundance = (merged / sum(merged)) * 100,
    # Reorder the 'Species' factor so the plot is sorted from high to low
    Species = fct_reorder(Species, Relative_Abundance)
  )

# 3. GENERATE PLOT
p <- ggplot(df_abundance, aes(x = Relative_Abundance, y = Species)) +
  geom_col(fill = "steelblue", width = 0.7) + # Create bars
  
  # Add percentage text labels at the end of each bar
  geom_text(aes(label = sprintf("%.1f%%", Relative_Abundance)), 
            hjust = -0.1, size = 3.5) +
  
  # Adjust x-axis to make room for the labels
  scale_x_continuous(limits = c(0, max(df_abundance$Relative_Abundance) * 1.15)) +
  
  # Labels and Theme
  labs(
    title = "Microbial Relative Abundance",
    subtitle = "Proportion of total reads per species",
    x = "Relative Abundance (%)",
    y = NULL # Hides the y-axis label "Species" since it's obvious
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines for cleaner look
    axis.text.y = element_text(size = 10, color = "black")
  )

# 4. VIEW AND SAVE
print(p)
```
