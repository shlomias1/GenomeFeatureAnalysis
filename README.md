# Genome Feature Analysis

## Overview
This project analyzes genome data, specifically from Bacillus subtilis, using GenBank and UniProt datasets. The script processes genome annotations, calculates statistical metrics, and visualizes results to gain insights into genome composition, gene distributions, and sequence characteristics.

## Features
- **Element Counting**: Counts various genome elements such as genes, CDS, tRNA, and rRNA.
- **Pseudo-Gene Identification**: Detects and counts pseudo-genes.
- **Gene Length Analysis**: Extracts gene lengths and visualizes distributions.
- **CA% Analysis**: Calculates the percentage of C and A nucleotides in the genome and in protein-coding genes.
- **Hydrophobic Amino Acid Analysis**: Examines transmembrane regions and their hydrophobic content.
- **Cross-referencing GenBank and UniProt**: Compares gene names between datasets and visualizes overlaps using a Venn diagram.
- **Data Export**: Saves processed genome data to a CSV file.

## Dependencies
- `biopython`
- `pandas`
- `numpy`
- `matplotlib`
- `matplotlib_venn`
- `requests`
- `gzip`
- `io`

## Installation
Ensure you have the required dependencies installed:
```bash
pip install biopython pandas numpy matplotlib matplotlib-venn requests
```

## Usage
1. Upload a GenBank file (`.gb` format) and a UniProt `.tsv.gz` file when prompted.
2. The script processes the data and performs various analyses.
3. Outputs are displayed, including tables and plots.
4. The final dataset is saved as `gene_info.csv` for further investigation.

## Functions
### `count_genome_elements(file_path)`
Counts occurrences of different genome elements and returns a dictionary.

### `count_pseudo_genes(file_path)`
Identifies and counts pseudo-genes.

### `extract_gene_lengths(file_path, gene_type="CDS")`
Extracts gene lengths for protein-coding genes, pseudo-genes, or other specified types.

### `calculate_ca_percentage(file_path)`
Computes the C and A nucleotide percentage across the entire genome.

### `calculate_ca_percentage_per_gene(file_path)`
Computes CA% per protein-coding gene.

### `compare_ca_percentage(genome_ca, coding_genes_ca)`
Compares CA% of the entire genome vs. protein-coding genes and provides insights.

### `top_bottom_ca_genes(file_path, top_n=5)`
Finds the top 5 and bottom 5 genes based on CA%.

### `plot_histogram(data, title, color="blue")`
Plots histograms for gene length distributions.

### `fetch_all_uniprot_data(query, fields)`
Retrieves data from UniProt API for Bacillus subtilis with relevant fields.

### `fetch_genbank_genes(file_path)`
Extracts gene names from a GenBank file.

### `cross_reference_genes(genbank_genes, df_uniprot)`
Compares gene names between GenBank and UniProt and generates a Venn diagram.

### `analyze_hydrophobic_content(df_uniprot, hydrophobic_amino_acids)`
Analyzes the percentage of hydrophobic amino acids in transmembrane sequences.

### `save_gene_info_to_csv(file_path, output_file="gene_info.csv")`
Saves analyzed gene data into a CSV file.

## Expected Outputs
- Tables displaying genome element counts, gene statistics, and CA% metrics.
- Histograms for gene length distributions.
- Venn diagrams comparing GenBank and UniProt genes.
- CSV file containing gene details and calculated metrics.

## Notes
- Some genes in UniProt may have multiple names, requiring careful cross-referencing.
- Not all proteins have transmembrane regions; analysis accounts for missing data.
- Data retrieval from UniProt API is limited to batches and may require pagination.

## License
This project is open-source under the MIT License.
