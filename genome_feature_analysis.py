import io
import gzip
import re
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from google.colab import files
from matplotlib_venn import venn2
from IPython.display import display

class GenomeParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.records = self._load_records()

    def _load_records(self):
        with open(self.file_path, "r") as file:
            return list(SeqIO.parse(file, "genbank"))

    def count_genome_elements(self):
        element_counts = {}
        for record in self.records:
            for feature in record.features:
                element_type = feature.type
                element_counts[element_type] = element_counts.get(element_type, 0) + 1
        return element_counts

    def count_pseudo_genes(self):
        return sum(1 for record in self.records for feature in record.features 
                   if feature.type == "gene" and "pseudo" in feature.qualifiers)

    def extract_pseudo_genes(self, gene_type="gene"):
        return [len(feature.location) for record in self.records for feature in record.features if feature.type == gene_type and "pseudo" in feature.qualifiers]

    def extract_gene_lengths(self, gene_type="CDS"):
        return [len(feature.location) for record in self.records for feature in record.features if feature.type == gene_type]

    def fetch_gene_names(self):
        return {feature.qualifiers.get("gene", ["Unknown"])[0] for record in self.records for feature in record.features if feature.type == "CDS"}

    def calculate_ca_percentage(self):
        total_length = 0
        ca_count = 0
        for record in self.records:
            genome_sequence = str(record.seq).upper()
            total_length += len(genome_sequence)
            ca_count += genome_sequence.count("C") + genome_sequence.count("A")
        ca_percentage = (ca_count / total_length) * 100 if total_length > 0 else 0
        return ca_percentage, total_length, ca_count
    
    def calculate_ca_percentage_per_gene(self):
        ca_percentages = []
        for record in self.records:
            for feature in record.features:
                if feature.type == "CDS":
                    gene_seq = feature.location.extract(record).seq.upper()
                    gene_length = len(gene_seq)
                    if gene_length > 0:
                        ca_count = gene_seq.count("C") + gene_seq.count("A")
                        ca_percentages.append((ca_count / gene_length) * 100)
        avg_ca_percentage = np.mean(ca_percentages) if ca_percentages else 0
        return avg_ca_percentage, ca_percentages
    
    def compare_ca_percentage(self, genome_ca, coding_genes_ca):
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
    
    def top_bottom_ca_genes(self, top_n=5):
        gene_info = []
        for record in self.records:
            for feature in record.features:
                if feature.type == "CDS" and feature.qualifiers.get("gene", [None])[0] is not None:
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
        display(df_top)
        print("\n===== Top 5 Genes with Lowest CA% =====")
        df_bottom = pd.DataFrame(sorted_genes[-top_n:])
        display(df_bottom)

    def save_gene_info_to_csv(self, output_file="gene_info.csv"):
        gene_info = []
        for record in self.records:
            for feature in record.features:
                if "gene" in feature.qualifiers or feature.type in {"CDS", "rRNA", "tRNA"}:
                    gene_seq = feature.location.extract(record).seq.upper()
                    gene_length = len(gene_seq)
                    ca_count = gene_seq.count("C") + gene_seq.count("A") if gene_length > 0 else 0
                    ca_percentage = (ca_count / gene_length) * 100 if gene_length > 0 else 0
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

class GenomeAnalysis:
    @staticmethod
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

    @staticmethod
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

class UniProtFetcher:
    def __init__(self, query, fields):
        self.query = query
        self.fields = fields
        self.base_url = "https://rest.uniprot.org/uniprotkb/search"
        self.df_uniprot = None

    def fetch_data(self):
        params = {"query": self.query, "fields": self.fields, "format": "tsv", "size": 500}
        response = requests.get(self.base_url, params=params)
        if response.status_code == 200:
            total_results = int(response.headers.get("X-Total-Results", 0))
            if total_results > 500:
                print(f"Total records in UniProt: {total_results}. Please upload the dataset manually.")
                uploaded = files.upload()
                file_path = next(iter(uploaded))
                with gzip.open(file_path, "rt") as f:
                    return pd.read_csv(f, sep="\t")
            return pd.read_csv(io.StringIO(response.text), sep="\t")
        return self.df_uniprot

    def find_top_bottom_proteins(self, top_n=5):
        if self.df_uniprot is None or "Length" not in self.df_uniprot.columns:
            print("No valid UniProt data available.")
            return None, None
        df_uniprot_sorted = self.df_uniprot.sort_values(by="Length", ascending=False)
        top_proteins = df_uniprot_sorted.head(top_n)
        bottom_proteins = df_uniprot_sorted.tail(top_n)
        print("\n===== 🔝 Top 5 Longest Proteins =====")
        display(top_proteins[["Protein names", "Length"]])
        print("\n===== 🔻 Top 5 Shortest Proteins =====")
        display(bottom_proteins[["Protein names", "Length"]])
        return top_proteins, bottom_proteins
    
class GeneComparison:
    @staticmethod
    def cross_reference_genes(genbank_genes, uniprot_genes):
        plt.figure(figsize=(6,6))
        venn = venn2([genbank_genes, uniprot_genes], ('GenBank', 'UniProt'))
        only_in_genbank = int(venn.get_label_by_id('10').get_text().strip())
        only_in_uniprot = int(venn.get_label_by_id('01').get_text().strip())
        common_genes = int(venn.get_label_by_id('11').get_text().strip())
        plt.title("Comparison of Gene Overlap Between GenBank and UniProt")
        plt.show()
        return only_in_genbank, only_in_uniprot, common_genes

class TransmembraneAnalysis:
    @staticmethod
    def analyze_transmembrane_regions(df_uniprot):
        transmembrane_lengths = []
        transmembrane_proteins = 0
        for _, row in df_uniprot.iterrows():
            if pd.notna(row["Transmembrane"]):
                transmembrane_proteins += 1
                matches = re.findall(r"(\d+)\.\.(\d+)", row["Transmembrane"])
                for match in matches:
                    start, end = map(int, match)
                    transmembrane_lengths.append(end - start + 1)
        num_regions = len(transmembrane_lengths)
        avg_length = np.mean(transmembrane_lengths) if transmembrane_lengths else 0
        min_length = np.min(transmembrane_lengths) if transmembrane_lengths else 0
        max_length = np.max(transmembrane_lengths) if transmembrane_lengths else 0
        print("\n===== 📊 Transmembrane Analysis =====")
        print(f"🟢 Proteins with transmembrane region: {transmembrane_proteins}")
        print(f"🟢 Total transmembrane regions: {num_regions}")
        print(f"📏 Minimum length: {min_length}, Maximum: {max_length}, Average: {avg_length:.2f}")
        return transmembrane_proteins, num_regions, transmembrane_lengths

class TransmembraneAnalysis:
    @staticmethod
    def analyze_transmembrane_regions(df_uniprot):
        transmembrane_lengths = []
        transmembrane_proteins = 0
        for _, row in df_uniprot.iterrows():
            if pd.notna(row["Transmembrane"]):
                transmembrane_proteins += 1
                matches = re.findall(r"(\d+)\.\.(\d+)", row["Transmembrane"])
                for match in matches:
                    start, end = map(int, match)
                    transmembrane_lengths.append(end - start + 1)
        num_regions = len(transmembrane_lengths)
        avg_length = np.mean(transmembrane_lengths) if transmembrane_lengths else 0
        min_length = np.min(transmembrane_lengths) if transmembrane_lengths else 0
        max_length = np.max(transmembrane_lengths) if transmembrane_lengths else 0
        print("\n===== 📊 Transmembrane Analysis =====")
        print(f"🟢 Proteins with transmembrane region: {transmembrane_proteins}")
        print(f"🟢 Total transmembrane regions: {num_regions}")
        print(f"📏 Minimum length: {min_length}, Maximum: {max_length}, Average: {avg_length:.2f}")
        return transmembrane_proteins, num_regions, transmembrane_lengths

    @staticmethod
    def analyze_hydrophobic_content(df_uniprot, hydrophobic_amino_acids):
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
        return avg_hydrophobic, hydrophobic_percentages

class GenomePipeline:
    def __init__(self, query, fields):
        self.file_path = self.upload_file()
        self.query = query
        self.fields = fields
        self.genome_parser = GenomeParser(self.file_path)
        self.uniprot_data = None

    def upload_file(self):
        uploaded = files.upload()
        return next(iter(uploaded))

    def run_analysis(self):
        element_counts = self.genome_parser.count_genome_elements()
        df = pd.DataFrame(list(element_counts.items()), columns=["Element Type", "Count"])
        display(df)
        GenomeAnalysis.plot_bar_chart(df, "Element Type", "Count", "Genome Element Distribution")

        pseudo_gene_count = self.genome_parser.count_pseudo_genes()
        pseudo_gene = self.genome_parser.extract_pseudo_genes()
        print(f"Total pseudo-genes: {pseudo_gene_count}")
        total_gene = self.genome_parser.extract_gene_lengths("gene")
        print(f"Total genes: {len(total_gene)}")

        protein_coding_gene_lengths = self.genome_parser.extract_gene_lengths("CDS")
        print(f"\nAverage Protein-Coding Gene Length: {np.mean(protein_coding_gene_lengths):.2f} bp")
        print(f"Average Pseudo Gene Length: {np.mean(pseudo_gene_count):.2f} bp")
        print(f"Average All Genes Length: {np.mean(total_gene):.2f} bp")
        GenomeAnalysis.plot_histogram(protein_coding_gene_lengths, "Protein-Coding Genes Length Distribution", "blue")
        GenomeAnalysis.plot_histogram(pseudo_gene, "Pseudo Genes Length Distribution", "orange")
        GenomeAnalysis.plot_histogram(total_gene, "All Genes Length Distribution", "purple")

        ca_percentage, total_length, ca_count = self.genome_parser.calculate_ca_percentage()
        print("\n===== CA Percentage in Genome =====")
        print(f"Genome Length: {total_length} bp")
        print(f"Total C + A Count: {ca_count}")
        print(f"CA Percentage: {ca_percentage:.2f}%")   
        avg_ca_percentage, ca_percentages = self.genome_parser.calculate_ca_percentage_per_gene()
        print("\n===== CA Percentage for Protein-Coding Genes =====")
        print(f"Total Protein-Coding Genes: {len(ca_percentages)}")
        print(f"Average CA Percentage: {avg_ca_percentage:.2f}%")
        self.genome_parser.compare_ca_percentage(ca_percentage, avg_ca_percentage)
        GenomeAnalysis.plot_histogram(ca_percentages, "CA% Distribution in Protein-Coding Genes", "skyblue")
        self.genome_parser.top_bottom_ca_genes()   
        self.genome_parser.save_gene_info_to_csv()  
        files.download("gene_info.csv")   

        uniprot_fetcher = UniProtFetcher(self.query, self.fields)
        self.uniprot_data = uniprot_fetcher.fetch_data()

        genbank_genes = self.genome_parser.fetch_gene_names()
        uniprot_genes = set(self.uniprot_data["Gene Names"].dropna())
        only_in_genbank, only_in_uniprot, common_genes = GeneComparison.cross_reference_genes(genbank_genes, uniprot_genes)
        print(f"Total genes in GenBank: {len(genbank_genes)}")
        print(f"Total genes in UniProt: {len(uniprot_genes)}")
        print(f"Genes only in GenBank: {only_in_genbank}")
        print(f"Genes only in UniProt: {only_in_uniprot}")
        print(f"Common genes in both: {common_genes}")
        top_proteins, bottom_proteins = uniprot_fetcher.find_top_bottom_proteins(top_n=5)

        transmembrane_proteins, num_regions, transmembrane_lengths = TransmembraneAnalysis.analyze_transmembrane_regions(self.uniprot_data)
        GenomeAnalysis.plot_histogram(transmembrane_lengths, "Transmembrane Region Lengths", "purple")
        hydrophobic_amino_acids = {"A", "P", "L", "I", "V", "M", "F", "W"}
        avg_hydrophobic, hydrophobic_percentages = TransmembraneAnalysis.analyze_hydrophobic_content(self.uniprot_data, hydrophobic_amino_acids)
        GenomeAnalysis.plot_histogram(hydrophobic_percentages, "Distribution of Hydrophobic Amino Acid Content in Transmembrane Sequences", "yellow")

if __name__ == "__main__":
    query = "organism_name:Bacillus AND reviewed:true"
    fields = "accession,id,protein_name,gene_names,organism_name,length,ft_transmem,gene_primary,organism_id,sequence"
    pipeline = GenomePipeline(query, fields)
    pipeline.run_analysis()
