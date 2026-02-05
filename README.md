# Nanopore-workflow-user-guide
This is page is meant to serve as a detailed user guide for beginner-coders working with nanopore.

## Access to cluster
The workflow can be best run on the SILS computer cluster Crunchomics HPC to get results fast while not occupying your computer's own computational resources. You may request an account by emailing W.C.deLeeuw@uva.nl mentioning your UvAnetID and purpose of using the cluster. Now follow the steps provided by the Crunchomics documentation to get started https://crunchomics-documentation.readthedocs.io/en/latest. To help yourself understand basic navigation commands in bash please refer to the following cheat sheet:
<img width="1536" height="1086" alt="image" src="https://github.com/user-attachments/assets/800de4d7-ab1c-43b7-9eaa-d3974b44bfac" />


## About NaMeco
NaMeco is the pipeline we will use to analyze our data and has to be installed using conda. If you have followed the Crunchomics documentation and suggested setups you should have conda installed. For more information about NaMeco always refer to https://github.com/timyerg/NaMeco for the latest updates. If you have any issues try to ask an experienced colleague for help first, but otherwise the developer is very quick to respond on GitHub to any issues you may encounter.

## Setting up your data
First I would recommend to make a folder for your nanopore data. In MobaXTerm, this can be done interactively in the navigation panel on the left. Please note that this is for quick navigating and will not set your working directory. You can use 'cd nanopore' to move to your nanopore directory (refer to cheat sheet) and you can go 1 level back by using 'cd ..' . In your nanopore folder, you can now make a new folder for each barcode as you can't upload folders, only files. Also remove the 0 from barcode names so barcode1 not barcode01. Now upload the contents of all your barcode folders (see icon in navigation pannel) to their respective folder.

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
- The amount of threads (computing power) the cluster will allocate to the analysis. More threads = faster but also a higher burden on the cluster. In my personal experience setting this to 100 threads falls easily within the socially accepted cluster burden while maintaining high speed analysis. You can play around with this setting depending on the size of your dataset and experience in the runtime.
- Forward primer setting. Default primers from NaMeco may be different from the primers we use in the lab so its helpful to specify the primers that were used
- Reverse primer setting idem

During our test runs we found that NaMeco does not function properly when ran from the console, so we have to submit a script to slurm. In a text file copy and paste the following information and adjust file paths to your own:
```
#!/bin/bash
#SBATCH --job-name=NaMeco
#SBATCH --output=logs/%x_%j.out
#SBATCH --cpus-per-task 64
#SBATCH --mem=120G
#SBATCH --time=03:30:00

RunName=Test1Rick
in_dir=/zfs/omics/personal/15827127/nanopore/test1/mergedreads/
out_dir=/zfs/omics/personal/15827127/nanopore/test1/results/${RunName}/

mkdir -p $out_dir

echo "Starting NaMeco:$(date)"

source /home/15827127/miniconda3/etc/profile.d/conda.sh
conda activate /home/15827127/miniconda3/envs/NaMeco/

srun nameco \
  --inp_dir ${in_dir} \
  --out_dir ${out_dir} \
  --threads ${SLURM_CPUS_PER_TASK} \
  --primer_F AGAGTTTGATCCTGGCTCAG \
  --primer_R GGTTACCTTGTTACGACTT

echo "Finished NaMeco:$(date)"   
```
Now save the file as NaMeco.slurm on your computer and upload it to your nanopore folder. To start the analysis, you can submit the script to the cluster by typing.

```
sbatch NaMeco.slurm
```
You can check whether your script is running by typing

```
squeue
```
You should check again after about a minute to make sure your script wasn't terminated. For navigating slurm commands please refer to https://crunchomics-documentation.readthedocs.io/en/latest/slurm_overview.html

# Data analysis
Congrats! You've done it!

Your results should be in your results folder, then final_output. Further analysis is described down below and can be most conveniently run from the Rmarkdown file on this GitHub page.

## Figure making
Lets start with setting up our data for figures. Select the family, genus and species count .tsv files and download them to your computer. Then, run the following script in R. Make sure to adjust the file paths.
```
library(readr)
library(tibble)

Family_counts <- read_tsv(
  'C:/Users/rickh/OneDrive - UvA/202511_RickHoogendijk/NanoporeWorkflow/barcode1/Family_counts.tsv',
  show_col_types = FALSE
) %>%
  column_to_rownames(var = colnames(.)[1])


cat(bold(blue("\n Dimensions Family counts: \n")),dim(Family_counts))

Genus_counts <- read_tsv(
  'C:/Users/rickh/OneDrive - UvA/202511_RickHoogendijk/NanoporeWorkflow/barcode1/Genus_counts.tsv',
  show_col_types = FALSE
) %>%
  column_to_rownames(var = colnames(.)[1])


cat(bold(blue("\n Dimensions Genus counts: \n")),dim(Genus_counts))

Species_counts <- read_tsv(
  'C:/Users/rickh/OneDrive - UvA/202511_RickHoogendijk/NanoporeWorkflow/barcode1/Species_counts.tsv',
  show_col_types = FALSE
) %>%
  column_to_rownames(var = colnames(.)[1])


cat(bold(blue("\n Dimensions species counts: \n")),dim(Species_counts))

mat <- as.matrix(Family_counts)

RA_Family <- sweep(mat, 2, colSums(mat, na.rm = TRUE), "/")
```
Now we will make a basic bar plot. You may need you to install some packages. You can do this interactively in the Rmarkdown file or by using 'install.packages'.
```
library("tidyverse")
library("vegan")
library("kableExtra")

# Convert to long format
ra_long <- RA_Family %>%
  as.data.frame() %>%
  rownames_to_column("Family") %>%
  pivot_longer(
    cols = -Family,
    names_to = "Run",
    values_to = "Abundance"
  )

# Convert to percent
ra_long <- ra_long %>%
  mutate(Percent = Abundance * 100)

# Group low-abundance species as "Other"
ra_long_grouped <- ra_long %>%
  mutate(
    Family = ifelse(Percent < 1, "Other", Family) #change if needed (is in 100% = 100 not = 1 range)
  ) %>%
  group_by(Run, Family) %>%
  summarize(Percent = sum(Percent), .groups = "drop")

# Reorder factor so Unclassified and Other are last
Family_levels <- unique(ra_long_grouped$Family)
Family_levels <- setdiff(Family_levels, c("Unclassified", "Other"))
Family_levels <- c(Family_levels, "Unclassified", "Other")
ra_long_grouped$Family <- factor(ra_long_grouped$Family, levels = Family_levels)

ggplot(ra_long_grouped, aes(Run, Percent, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Spectral") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

```

## BLASting unclassified/unexpected clusters
This is an optional check if you encountered any unclassified species or unexpected species, for example in a SynCom. First, we  will check in the clusters how many/which species could not be classified using the following command in R. Make sure to set the file paths to where you have yours saved.
```
library(crayon)
Clusters <- read_tsv('C:/Users/rickh/OneDrive - UvA/202511_RickHoogendijk/NanoporeWorkflow/barcode1/Taxonomy.tsv')
Unclassified_species <- Clusters[grep("unclassified", Clusters$Species),]

cat(bold(yellow("\n Number of unclassified species: \n")),length(Unclassified_species$Cluster),bold(yellow("\n Head unclassified species: \n")))
print(head(Unclassified_species))

cat("\n\n",
    green("Adjust cluster number in chunk below to extract fasta of cluster/unclassified bacteria of intrest \n"))

fasta_lines <- readLines('C:/Users/rickh/OneDrive - UvA/202511_RickHoogendijk/NanoporeWorkflow/barcode1/rep_seqs.fasta')

get_cluster_fasta <- function(fasta_lines, cluster_id) {
  
  # Find the header line of the cluster
  start_idx <- grep(paste0("^>", cluster_id, "$"), fasta_lines)
  if (length(start_idx) == 0) return(NULL)
  
  # Find all FASTA headers
  header_idx <- grep("^>", fasta_lines)
  
  # Find the next header after the cluster
  next_header <- header_idx[header_idx > start_idx][1]
  
  # Define end index
  end_idx <- ifelse(
    is.na(next_header),
    length(fasta_lines),
    next_header - 1
  )
  
  fasta_lines[start_idx:end_idx]
}
```
Next, we will retrieve the sequences from these clusters using the following commands. Make sure to adjust the cluster number to those in your unclassified table.
```
cluster_fasta <- get_cluster_fasta(fasta_lines, "Cluster_3")

cat(cluster_fasta, sep = "\n")
```
You may blast the resulting sequence to double check the species. An example is that in a trial run we found an unclassified Wolbachia spp. which turned out to be E. coli when BLASTed.
