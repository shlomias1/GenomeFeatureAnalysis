from Bio import SeqIO
from google.colab import files
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import IPython.display as display
from csv import Error
import requests
import io
from matplotlib_venn import venn2
import gzip
import re

def count_genome_elements(file_path):
    element_counts = {}  
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            for feature in record.features:
                element_type = feature.type 
                element_counts[element_type] = element_counts.get(element_type, 0) + 1
    return element_counts

def plot_bar_chart(df, x, y, title):
    plt.figure(figsize=(10, 5))
    plt.bar(df[x], df[y], color='skyblue', edgecolor='black', alpha=0.7)
    plt.title(title, fontsize=14)
    plt.xlabel(x, fontsize=12)
    plt.ylabel(y, fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.gca().spines[['top', 'right']].set_visible(False)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    for index, value in enumerate(df[y]):
        plt.text(index, value + max(df[y]) * 0.02, str(value), ha='center', fontsize=10)
    plt.show()

def count_pseudo_genes(file_path):
    pseudo_gene_count = 0
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            for feature in record.features:
                if feature.type == "gene" and "pseudo" in feature.qualifiers:
                    pseudo_gene_count += 1
    return pseudo_gene_count

def extract_gene_lengths(file_path, gene_type="CDS"):
    gene_lengths = []
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            for feature in record.features:
                if feature.type == gene_type:
                    gene_length = len(feature.location)
                    gene_lengths.append(gene_length)
    return gene_lengths

def plot_histogram(data, title, color="blue"):
    if not data:
        print(f"No data available for {title}")
        return
    plt.figure(figsize=(8, 5))
    plt.hist(data, bins=30, color=color, edgecolor="black", alpha=0.7)
    plt.title(title, fontsize=14)
    plt.xlabel("Length (bp)", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.show()

def calculate_ca_percentage(file_path):
    total_length = 0
    ca_count = 0
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            genome_sequence = str(record.seq).upper() 
            total_length += len(genome_sequence)
            ca_count += genome_sequence.count("C") + genome_sequence.count("A")
    if total_length == 0:
        return 0  
    ca_percentage = (ca_count / total_length) * 100
    return ca_percentage, total_length, ca_count

def calculate_ca_percentage_per_gene(file_path):
    ca_percentages = []
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            for feature in record.features:
                if feature.type == "CDS":  
                    gene_seq = feature.location.extract(record).seq.upper()  
                    gene_length = len(gene_seq)

                    if gene_length > 0:
                        ca_count = gene_seq.count("C") + gene_seq.count("A")
                        ca_percentage = (ca_count / gene_length) * 100
                        ca_percentages.append(ca_percentage)
    if not ca_percentages:
        return 0, []
    avg_ca_percentage = np.mean(ca_percentages)
    return avg_ca_percentage, ca_percentages

def compare_ca_percentage(genome_ca, coding_genes_ca):
    print("\n===== CA% Comparison =====")
    print(f"CA% in the entire genome: {genome_ca:.2f}%")
    print(f"CA% in protein-coding genes: {coding_genes_ca:.2f}%")
    difference = abs(genome_ca - coding_genes_ca)
    if difference < 1:
        conclusion = "The CA% in protein-coding genes is nearly identical to the overall genome, suggesting a uniform nucleotide distribution."
    elif coding_genes_ca > genome_ca:
        conclusion = "The CA% in protein-coding genes is higher than in the genome, indicating a possible preference for C and A bases in coding regions."
    else:
        conclusion = "The CA% in protein-coding genes is lower than in the genome, suggesting that non-coding regions contain more C and A bases."
    print("\nConclusion:")
    print(conclusion)

def top_bottom_ca_genes(file_path, top_n=5):
    gene_info = []
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            for feature in record.features:
                if feature.type == "CDS":  
                    gene_seq = feature.location.extract(record).seq.upper()  
                    gene_length = len(gene_seq)
                    if gene_length > 0:
                        ca_count = gene_seq.count("C") + gene_seq.count("A")
                        ca_percentage = (ca_count / gene_length) * 100
                        gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                        start = feature.location.start
                        end = feature.location.end
                        strand = feature.location.strand
                        gene_info.append({
                            "Gene": gene_name,
                            "Start": start,
                            "End": end,
                            "Strand": strand,
                            "CA%": ca_percentage
                        })
    sorted_genes = sorted(gene_info, key=lambda x: x["CA%"], reverse=True)
    print("\n===== Top 5 Genes with Highest CA% =====")
    df_top = pd.DataFrame(sorted_genes[:top_n])
    display.display(df_top)
    print("\n===== Top 5 Genes with Lowest CA% =====")
    df_bottom = pd.DataFrame(sorted_genes[-top_n:])
    display.display(df_bottom)

def save_gene_info_to_csv(file_path, output_file="gene_info.csv"):
    gene_info = []
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            for feature in record.features:
                if "gene" in feature.qualifiers or feature.type in {"CDS", "rRNA", "tRNA"}:
                    gene_seq = feature.location.extract(record).seq.upper()  
                    gene_length = len(gene_seq)
                    if gene_length > 0:
                        ca_count = gene_seq.count("C") + gene_seq.count("A")
                        ca_percentage = (ca_count / gene_length) * 100
                    else:
                        ca_percentage = 0
                    gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                    gene_type = feature.type 
                    start = feature.location.start
                    end = feature.location.end
                    strand = feature.location.strand
                    gene_info.append({
                        "Gene Name": gene_name,
                        "Type": gene_type,
                        "Start": start,
                        "End": end,
                        "Strand": strand,
                        "Length (bp)": gene_length,
                        "CA%": ca_percentage
                    })
    df = pd.DataFrame(gene_info)
    df = df.sort_values(by="Start") 
    df.to_csv(output_file, index=False)
    print(f"\n✅ Gene information saved to {output_file}")

def fetch_all_uniprot_data(query, fields):
    uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": fields,
        "format": "tsv"    
    }
    response = requests.get(uniprot_url, params=params)
    if response.status_code == 200:
        df_uniprot = pd.read_csv(io.StringIO(response.text), sep="\t")
        total_results = response.headers.get("X-Total-Results")
        batch_size = response.headers.get("Link").split("size=")[-1].split("&")[0].split(">")[0]
        if int(batch_size) < int(total_results):
            return None
        return df_uniprot
    else:
        return None

def fetch_genbank_genes(file_path):
    gene_names = set()
    with open(file_path, "r") as file:
        records = SeqIO.parse(file, "genbank")
        for record in records:
            for feature in record.features:
                if feature.type == "CDS" and "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                    gene_names.add(gene_name)
    return gene_names

def cross_reference_genes(genbank_genes, uniprot_genes):
    plt.figure(figsize=(6,6))
    venn = venn2([genbank_genes, uniprot_genes], ('GenBank', 'UniProt'))
    only_in_genbank = int(venn.get_label_by_id('10').get_text().strip())
    only_in_uniprot = int(venn.get_label_by_id('01').get_text().strip())
    common_genes = int(venn.get_label_by_id('11').get_text().strip())
    plt.title("Comparison of Gene Overlap Between GenBank and UniProt")
    plt.show()
    return only_in_genbank, only_in_uniprot, common_genes

def find_top_bottom_proteins(df_uniprot, top_n=5):
    df_uniprot_sorted = df_uniprot.sort_values(by="Length", ascending=False)  
    top_proteins = df_uniprot_sorted.head(top_n) 
    bottom_proteins = df_uniprot_sorted.tail(top_n)  
    print("\n===== 🔝 Top 5 Longest Proteins =====")
    display(top_proteins[["Protein names", "Length"]])
    print("\n===== 🔻 Top 5 Shortest Proteins =====")
    display(bottom_proteins[["Protein names", "Length"]])
    return top_proteins, bottom_proteins

def analyze_transmembrane_regions(df_uniprot):
    transmembrane_lengths = [] 
    transmembrane_proteins = 0 
    for index, row in df_uniprot.iterrows():
        if pd.notna(row["Transmembrane"]):  
            transmembrane_proteins += 1  
            matches = re.findall(r"(\d+)\.\.(\d+)", row["Transmembrane"]) 
            for match in matches:
                start, end = map(int, match)  
                length = end - start + 1 
                transmembrane_lengths.append(length)
    num_regions = len(transmembrane_lengths)
    avg_length = np.mean(transmembrane_lengths) if transmembrane_lengths else 0
    min_length = np.min(transmembrane_lengths) if transmembrane_lengths else 0
    max_length = np.max(transmembrane_lengths) if transmembrane_lengths else 0
    print("\n===== 📊 Transmembrane Analysis =====")
    print(f"🟢 Proteins with transmembrane region: {transmembrane_proteins}")
    print(f"🟢 Total transmembrane regions: {num_regions}")
    print(f"📏 Minimum length: {min_length}, Maximum: {max_length}, Average: {avg_length:.2f}")
    if transmembrane_lengths:
        plot_histogram(transmembrane_lengths, "Distribution of Transmembrane Region Lengths", "purple")
    return transmembrane_proteins, num_regions, transmembrane_lengths

def analyze_hydrophobic_content(df_uniprot,hydrophobic_amino_acids):
    hydrophobic_percentages = []
    for index, row in df_uniprot.iterrows():
        if pd.notna(row["Transmembrane"]) and pd.notna(row["Sequence"]):
            matches = re.findall(r"(\d+)\.\.(\d+)", row["Transmembrane"])  
            for match in matches:
                start, end = map(int, match)
                transmembrane_seq = row["Sequence"][start-1:end]      
                if transmembrane_seq:
                    hydrophobic_count = sum(1 for aa in transmembrane_seq if aa in hydrophobic_amino_acids)
                    hydrophobic_percentage = (hydrophobic_count / len(transmembrane_seq)) * 100
                    hydrophobic_percentages.append(hydrophobic_percentage)
    avg_hydrophobic = np.mean(hydrophobic_percentages) if hydrophobic_percentages else 0
    min_hydrophobic = np.min(hydrophobic_percentages) if hydrophobic_percentages else 0
    max_hydrophobic = np.max(hydrophobic_percentages) if hydrophobic_percentages else 0
    print("\n===== 📊 Hydrophobic Amino Acid Analysis in Transmembrane Sequences =====")
    print(f"🟢 Average hydrophobic amino acid percentage: {avg_hydrophobic:.2f}%")
    print(f"📏 Minimum: {min_hydrophobic:.2f}%, Maximum: {max_hydrophobic:.2f}%")
    if hydrophobic_percentages:
        plot_histogram(hydrophobic_percentages, "Distribution of Hydrophobic Amino Acid Content in Transmembrane Sequences", "blue")
    return avg_hydrophobic, hydrophobic_percentages

def main():
    uploaded = files.upload()
    file_path = next(iter(uploaded)) 

    element_counts = count_genome_elements(file_path)
    df = pd.DataFrame(list(element_counts.items()), columns=["Element Type", "Count"])

    print("\n===== Counting elements in the Bacillus subtilis genome =====")
    display.display(df)
    plot_bar_chart(df, x="Element Type", y="Count", title="Distribution of Genome Element Counts")

    pseudo_gene_count = count_pseudo_genes(file_path)

    protein_coding_gene_lengths = extract_gene_lengths(file_path, "CDS")
    pseudo_gene_lengths = count_pseudo_genes(file_path)
    total_gene_lengths = extract_gene_lengths(file_path, "gene")
    print(f"\nTotal pseudo-genes: {pseudo_gene_count}")
    print(f"\nTotal genes: {len(total_gene_lengths)}")
    print(f"\nAverage Protein-Coding Gene Length: {np.mean(protein_coding_gene_lengths):.2f} bp")
    print(f"Average Pseudo Gene Length: {np.mean(pseudo_gene_lengths):.2f} bp")
    print(f"Average All Genes Length: {np.mean(total_gene_lengths):.2f} bp")

    plot_histogram(protein_coding_gene_lengths, "Protein-Coding Genes Length Distribution", "blue")
    plot_histogram(pseudo_gene_lengths, "Pseudo Genes Length Distribution", "red")
    plot_histogram(total_gene_lengths, "All Genes Length Distribution", "green")

    ca_percentage, total_length, ca_count = calculate_ca_percentage(file_path)
    print("\n===== CA Percentage in Bacillus subtilis Genome =====")
    print(f"Genome Length: {total_length} bp")
    print(f"Total C + A Count: {ca_count}")
    print(f"CA Percentage: {ca_percentage:.2f}%")

    avg_ca_percentage, ca_percentages = calculate_ca_percentage_per_gene(file_path)
    print("\n===== CA Percentage for Protein-Coding Genes =====")
    print(f"Total Protein-Coding Genes: {len(ca_percentages)}")
    print(f"Average CA Percentage: {avg_ca_percentage:.2f}%")
    compare_ca_percentage(ca_percentage, avg_ca_percentage)
    plot_histogram(ca_percentages, "CA% Distribution in Protein-Coding Genes", "purple")
    top_bottom_ca_genes(file_path)

    save_gene_info_to_csv(file_path, output_file="gene_info.csv")
    files.download("gene_info.csv")

    query = "organism_name:Bacillus AND reviewed:true"
    fields = "accession,id,protein_name,gene_names,organism_name,length,ft_transmem,gene_primary,organism_id,sequence"
    df_uniprot = fetch_all_uniprot_data(query,fields)
    if df_uniprot is None:
        uploaded = files.upload()
        uniprot_file_path = next(iter(uploaded))
        with gzip.open(uniprot_file_path, "rt") as f:
          df_uniprot = pd.read_csv(f, sep="\t")
    genbank_genes = fetch_genbank_genes(file_path)
    only_in_genbank, only_in_uniprot, common_genes = cross_reference_genes(genbank_genes, set(df_uniprot["Gene Names"].dropna()))
    print(f"Total genes in GenBank: {len(genbank_genes)}")
    print(f"Total genes in UniProt: {len(set(df_uniprot['Gene Names'].dropna()))}")
    print(f"Genes only in GenBank: {only_in_genbank}")
    print(f"Genes only in UniProt: {only_in_uniprot}")
    print(f"Common genes in both: {common_genes}")
    find_top_bottom_proteins(df_uniprot)
    analyze_transmembrane_regions(df_uniprot)
    hydrophobic_amino_acids = {"A", "P", "L", "I", "V", "M", "P", "T"}
    analyze_hydrophobic_content(df_uniprot, hydrophobic_amino_acids)

if __name__ == "__main__":
    main()
