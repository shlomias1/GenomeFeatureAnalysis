# Genome Analysis - README

## 📌 Introduction
This project provides an **object-oriented** implementation for **analyzing genomic data** from GenBank files. The script allows users to **extract genomic features**, analyze **pseudo-genes**, compute **nucleotide composition**, and cross-reference **UniProt protein data**.

The project uses **Biopython, Pandas, Matplotlib, and Requests** for processing genomic sequences and generating visualizations.

## 📦 Features
- **Parse GenBank files** to extract genomic elements
- **Analyze pseudo-genes** and coding sequences (CDS)
- **Calculate CA% composition** across the genome and protein-coding genes
- **Cross-reference genes with UniProt database** for further protein analysis
- **Visualize data** with histograms, bar charts, and Venn diagrams
- **Save extracted gene information to CSV**

## 🚀 Installation
Ensure you have the following Python libraries installed:
```bash
pip install biopython pandas matplotlib requests
```
If you're running this in Google Colab, `google.colab` is pre-installed.

## 📂 Project Structure
```
📦 Genome Analysis
├── genome_parser.py        # Handles parsing of GenBank files and genome feature extraction
├── genome_analysis.py      # Contains visualization and statistical analysis functions
├── uniprot_fetcher.py      # Handles fetching and processing UniProt data
├── transmembrane_analysis.py  # Analyzes transmembrane and hydrophobicity properties
├── genome_pipeline.py      # Main execution pipeline for genome analysis
├── README.md               # Documentation (this file)
```

## 🔧 Usage
### **1️⃣ Upload a GenBank File**
When running the script, the user is prompted to upload a GenBank file:
```python
uploaded = files.upload()
file_path = next(iter(uploaded))
```

### **2️⃣ Run the Analysis**
To execute the pipeline:
```python
query = "organism_name:Bacillus AND reviewed:true"
fields = "accession,id,protein_name,gene_names,organism_name,length,ft_transmem,gene_primary,organism_id,sequence"
pipeline = GenomePipeline(query, fields)
pipeline.run_analysis()
```

### **3️⃣ Output Results**
- **Genomic Features Count**: Displays the distribution of element types in the genome
- **Pseudo-Gene Analysis**: Computes the count and length distribution of pseudo-genes
- **CA% Computation**: Calculates C and A nucleotide frequency in coding and non-coding sequences
- **Gene Overlap Analysis**: Uses Venn diagram to compare genes in **GenBank** and **UniProt**
- **Protein Analysis**: Fetches UniProt data, finds longest/shortest proteins, and analyzes transmembrane regions

### **4️⃣ Save & Download Data**
The extracted gene information is saved to `gene_info.csv`:
```python
self.genome_parser.save_gene_info_to_csv()
files.download("gene_info.csv")
```

## 📊 Visualizations
The analysis generates various plots, including:
- **Genome Element Distribution (Bar Chart)**
- **Pseudo-Gene and CDS Length Distributions (Histograms)**
- **CA% Distribution in Protein-Coding Genes**
- **Transmembrane Region Lengths**
- **Hydrophobicity Distribution in Transmembrane Proteins**

## 🔍 Example Output
```
===== CA% Comparison =====
CA% in the entire genome: 54.32%
CA% in protein-coding genes: 50.12%
Conclusion: The CA% in protein-coding genes is lower than in the genome, suggesting that non-coding regions contain more C and A bases.
```

## 🛠 Future Improvements
- Support for multiple sequence file formats (FASTA, EMBL)
- Additional protein structural analysis (e.g., AlphaFold integration)
- More efficient parallelized data fetching from UniProt

## 📝 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
