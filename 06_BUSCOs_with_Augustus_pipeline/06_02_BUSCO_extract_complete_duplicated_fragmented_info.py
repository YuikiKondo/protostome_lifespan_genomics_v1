## Purpouse: gff files are not separately outputted for complete/duplicated/fragmented genes in running BUSCO with Augustus. To use the same scripts as processing BUSCO with Miniprot results in following analyses, create a table about which category each gene is in.
## e.g. output file:/home/username/protostome_lifespan_model_augustus/BUSCO_directory/BUSCO_all_species/BUSCO_GCA_029721355.1_KUTeg_genomic.fna/run_metazoa_odb10/busco_sequences/busco_gene_status_summary.csv)

## example output:
#(base) [username@usernumber busco_sequences]$ head busco_gene_status_summary.csv
#GeneID,Fragmented,Multi_copy,Single_copy
#102804at33208,0,0,1
#103854at33208,0,0,1
#107151at33208,0,0,1
#107574at33208,0,0,1
#108764at33208,0,0,1
#109109at33208,0,0,1
#111730at33208,0,0,1
#111746at33208,0,0,1
#114954at33208,0,0,1

import os
import csv

# Path to the parent directory containing all BUSCO_ directories
parent_dir = "/home/username/protostome_lifespan_model_augustus/BUSCO_directory/BUSCO_all_species"

# Categories to look into
categories = ["fragmented_busco_sequences", "multi_copy_busco_sequences", "single_copy_busco_sequences"]

# Loop through all BUSCO_ directories
for species_dir in os.listdir(parent_dir):
    species_path = os.path.join(parent_dir, species_dir)
    
    if os.path.isdir(species_path) and species_dir.startswith("BUSCO_"):

        busco_seq_path = os.path.join(species_path, "run_metazoa_odb10", "busco_sequences")
        if not os.path.exists(busco_seq_path):
            print(f"Skipping {species_dir} – no busco_sequences directory found.")
            continue

        gene_dict = {}

        for category in categories:
            category_path = os.path.join(busco_seq_path, category)
            if not os.path.exists(category_path):
                print(f"Warning: {category_path} not found.")
                continue
            for file in os.listdir(category_path):
                if file.endswith(".faa"):
                    gene_id = file.replace(".faa", "")
                    if gene_id not in gene_dict:
                        gene_dict[gene_id] = {"Fragmented": 0, "Multi_copy": 0, "Single_copy": 0}
                    if category.startswith("fragmented"):
                        gene_dict[gene_id]["Fragmented"] = 1
                    elif category.startswith("multi_copy"):
                        gene_dict[gene_id]["Multi_copy"] = 1
                    elif category.startswith("single_copy"):
                        gene_dict[gene_id]["Single_copy"] = 1

        # Output CSV to each busco_sequences directory
        output_csv_path = os.path.join(busco_seq_path, "busco_gene_status_summary.csv")
        with open(output_csv_path, mode="w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["GeneID", "Fragmented", "Multi_copy", "Single_copy"])
            for gene_id, status in sorted(gene_dict.items()):
                writer.writerow([gene_id, status["Fragmented"], status["Multi_copy"], status["Single_copy"]])
        
        print(f"✓ Written summary for {species_dir} → {output_csv_path}")
