# -*- coding: utf-8 -*-
from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import sys
import argparse
from dotenv import load_dotenv

# Загружаем переменные из .env файла (не попадает в GitHub)
load_dotenv()
ncbi_email = os.getenv("NCBI_EMAIL", "")
if not ncbi_email:
    print("[!] WARNING: NCBI_EMAIL not set in .env file.")
    print("    Copy .env.example to .env and add your email.")
    print("    NCBI may block requests without a valid email.\n")
Entrez.email = ncbi_email


def download_genome(accession, output_dir="data"):
    """
    Downloads bacterial genome from NCBI by accession number.
    Saves as GenBank file to data/ folder.
    Returns SeqRecord object.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{accession}.gb")

    if os.path.exists(output_path):
        print(f"[OK] Genome file found locally: {output_path}")
        record = SeqIO.read(output_path, "genbank")
        return record

    print(f"[-->] Downloading genome {accession} from NCBI...")
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        SeqIO.write(record, output_path, "genbank")
        print(f"[OK] Saved: {output_path}")
        print(f"[OK] Organism: {record.annotations.get('organism', 'unknown')}")
        print(f"[OK] Genome length: {len(record.seq):,} nucleotides")
        return record
    except Exception as e:
        print(f"[!] Could not download genome: {e}")
        return None


def generate_demo_data(output_dir="results"):
    """
    Generates realistic simulated AMR gene data for demo mode.
    Uses real gene names from CARD database for S. aureus.
    Returns DataFrame with simulated resistance genes.
    """
    os.makedirs(output_dir, exist_ok=True)

    demo_hits = [
        {"subject": "mecA|beta-lactam|methicillin resistance",
         "identity": 99.2, "length": 2007, "evalue": 0.0, "bitscore": 3701},
        {"subject": "blaZ|beta-lactam|penicillinase",
         "identity": 97.8, "length": 846, "evalue": 0.0, "bitscore": 1561},
        {"subject": "tetM|tetracycline|ribosomal protection",
         "identity": 94.1, "length": 1920, "evalue": 0.0, "bitscore": 3421},
        {"subject": "aacA-aphD|aminoglycoside|bifunctional enzyme",
         "identity": 91.5, "length": 1449, "evalue": 0.0, "bitscore": 2601},
        {"subject": "ermC|macrolide|rRNA methyltransferase",
         "identity": 89.3, "length": 735, "evalue": 1e-180, "bitscore": 1301},
        {"subject": "dfrA|trimethoprim|dihydrofolate reductase",
         "identity": 85.7, "length": 474, "evalue": 1e-120, "bitscore": 861},
        {"subject": "fusB|fusidic acid|elongation factor protection",
         "identity": 82.4, "length": 579, "evalue": 1e-100, "bitscore": 1041},
        {"subject": "norA|fluoroquinolone|efflux pump",
         "identity": 79.6, "length": 1203, "evalue": 1e-90, "bitscore": 1561},
        {"subject": "vanA|glycopeptide|D-Ala-D-Lac ligase",
         "identity": 76.3, "length": 1032, "evalue": 1e-80, "bitscore": 1301},
        {"subject": "cfr|phenicol|rRNA methyltransferase",
         "identity": 72.1, "length": 996, "evalue": 1e-60, "bitscore": 901},
    ]

    df = pd.DataFrame(demo_hits)

    blast_path = os.path.join(output_dir, "blast_hits.tsv")
    df_blast = pd.DataFrame({
        "query": ["NC_007795.1"] * len(df),
        "subject": df["subject"],
        "identity": df["identity"],
        "length": df["length"],
        "mismatch": [0] * len(df),
        "gapopen": [0] * len(df),
        "qstart": [1] * len(df),
        "qend": df["length"],
        "sstart": [1] * len(df),
        "send": df["length"],
        "evalue": df["evalue"],
        "bitscore": df["bitscore"],
    })
    df_blast.to_csv(blast_path, sep="\t", header=False, index=False)
    print(f"[OK] Demo BLAST data generated: {blast_path}")
    return df


def find_resistance_genes(blast_results_file, min_identity=70.0):
    """
    Reads BLAST output file and filters resistance genes by identity.
    Accepts TSV file in BLAST -outfmt 6 format.
    Returns sorted DataFrame with genes above identity threshold.
    Raises FileNotFoundError if file does not exist.
    Raises ValueError if file is empty or has wrong format.
    """
    if not os.path.exists(blast_results_file):
        raise FileNotFoundError(
            f"BLAST results file not found: {blast_results_file}\n"
            "  --> Run in demo mode: python project.py --demo"
        )

    columns = [
        "query", "subject", "identity", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"
    ]

    try:
        df = pd.read_csv(blast_results_file, sep="\t", names=columns)
    except Exception as e:
        raise ValueError(f"Could not read BLAST file: {e}") from e

    if df.empty:
        raise ValueError("BLAST results file is empty.")

    hits = df[df["identity"] >= min_identity].copy()
    hits = hits.sort_values("identity", ascending=False)

    print(f"[OK] Total BLAST hits: {len(df)}")
    print(f"[OK] After filtering (>= {min_identity}%): {len(hits)}")
    return hits


def parse_subject(subject_str):
    """
    Parses BLAST subject field in format 'gene|class|mechanism'.
    Returns (gene_name, antibiotic_class, mechanism) tuple.
    Gracefully handles subjects without pipe separators.
    """
    parts = subject_str.split("|")
    gene = parts[0].strip() if len(parts) > 0 else subject_str
    ab_class = parts[1].strip() if len(parts) > 1 else "unknown"
    mechanism = parts[2].strip() if len(parts) > 2 else ""
    return gene, ab_class, mechanism


def visualize_results(hits, output_dir="results", top_n=10):
    """
    Builds horizontal bar chart of top-N resistance genes.
    Colors bars by ANTIBIOTIC CLASS (not just identity level).
    Labels show 'gene_name (class)' format for readability.
    Saves publication-ready PNG to results/ folder.
    Returns path to saved file or None if no data.
    """
    os.makedirs(output_dir, exist_ok=True)

    if hits.empty:
        print("[!] No data to visualize.")
        return None

    # Цветовая палитра по классам антибиотиков
    CLASS_COLORS = {
        "beta-lactam":      "#c0392b",   # красный
        "tetracycline":     "#e67e22",   # оранжевый
        "aminoglycoside":   "#8e44ad",   # фиолетовый
        "macrolide":        "#2980b9",   # синий
        "fluoroquinolone":  "#16a085",   # бирюзовый
        "glycopeptide":     "#d35400",   # тёмно-оранжевый
        "trimethoprim":     "#27ae60",   # зелёный
        "phenicol":         "#7f8c8d",   # серый
        "fusidic acid":     "#f39c12",   # жёлтый
    }
    DEFAULT_COLOR = "#34495e"  # тёмно-серый для неизвестных классов

    top_hits = hits.head(top_n).copy()

    # Парсим subject и строим красивые метки
    parsed = top_hits["subject"].apply(parse_subject)
    top_hits["gene_name"] = parsed.apply(lambda x: x[0])
    top_hits["ab_class"] = parsed.apply(lambda x: x[1])
    top_hits["mechanism"] = parsed.apply(lambda x: x[2])
    top_hits["label"] = top_hits.apply(
        lambda r: f"{r['gene_name']}  ({r['ab_class']})", axis=1
    )

    colors = [
        CLASS_COLORS.get(row["ab_class"].lower(), DEFAULT_COLOR)
        for _, row in top_hits.iterrows()
    ]

    fig, ax = plt.subplots(figsize=(14, 8))
    fig.patch.set_facecolor("#f0f4f8")
    ax.set_facecolor("#f0f4f8")

    bars = ax.barh(
        top_hits["label"],
        top_hits["identity"],
        color=colors,
        edgecolor="white",
        linewidth=1.0,
        height=0.6
    )

    for bar, val, mech in zip(bars, top_hits["identity"], top_hits["mechanism"]):
        # Процент внутри бара
        ax.text(
            bar.get_width() - 1.2,
            bar.get_y() + bar.get_height() / 2,
            f"{val:.1f}%",
            va="center", ha="right",
            color="white", fontweight="bold", fontsize=10
        )
        # Механизм справа от бара (маленький текст)
        if mech:
            short_mech = mech[:30] + "..." if len(mech) > 30 else mech
            ax.text(
                bar.get_width() + 0.4,
                bar.get_y() + bar.get_height() / 2,
                short_mech,
                va="center", ha="left",
                color="#555555", fontsize=8, style="italic"
            )

    ax.set_xlabel("Identity to CARD reference (%)", fontsize=12, labelpad=10)
    ax.set_title(
        f"Top-{top_n} Antibiotic Resistance Genes Detected\n"
        f"Staphylococcus aureus NCTC 8325 vs. CARD Database",
        fontsize=13, fontweight="bold", pad=15
    )
    ax.set_xlim(60, 130)
    ax.invert_yaxis()
    ax.axvline(x=90, color="#cc3333", linestyle="--",
               linewidth=1, alpha=0.6, label="90% threshold")
    ax.axvline(x=70, color="#cc8800", linestyle="--",
               linewidth=1, alpha=0.6, label="70% threshold")
    ax.legend(loc="lower right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=10)
    plt.tight_layout()

    output_path = os.path.join(output_dir, "top_resistance_genes.png")
    plt.savefig(output_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()

    print(f"[OK] Chart saved: {output_path}")
    return output_path


def parse_args():
    """
    Parses command-line arguments using argparse.
    Supports --demo mode, custom accession number and identity threshold.
    """
    parser = argparse.ArgumentParser(
        prog="project.py",
        description="AMR Genome Analyzer — Veterinary Antimicrobial Resistance Profiler",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python project.py --demo\n"
            "  python project.py --accession NC_007795.1\n"
            "  python project.py --accession NC_007795.1 --identity 80"
        )
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run with simulated data (no BLAST or NCBI needed)"
    )
    parser.add_argument(
        "--accession",
        type=str,
        default="NC_007795.1",
        metavar="ACCESSION",
        help="NCBI accession number (default: NC_007795.1)"
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=70.0,
        metavar="THRESHOLD",
        help="Minimum identity %% for filtering (default: 70.0)"
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        metavar="N",
        help="Number of top genes to display on chart (default: 10)"
    )
    return parser.parse_args()


def main():
    """
    Main entry point. Supports two modes:
    --demo : runs with simulated realistic data (no BLAST needed)
    normal : downloads real genome and reads BLAST results file
    Use --help to see all available options.
    """
    args = parse_args()

    print("=" * 55)
    print("   AMR Genome Analyzer v1.0")
    print("   Veterinary Antimicrobial Resistance Profiler")
    print("=" * 55)

    if args.demo:
        print("\n[--> DEMO MODE]")
        print("[-->] Simulating AMR profile for")
        print("    Staphylococcus aureus NCTC 8325\n")
        generate_demo_data()
        blast_file = os.path.join("results", "blast_hits.tsv")

    else:
        print(f"\n[-->] Using genome accession: {args.accession}")
        record = download_genome(args.accession)
        if record is None:
            print("[!] Genome download failed. Try: python project.py --demo")
            sys.exit(1)

        blast_file = os.path.join("results", "blast_hits.tsv")
        if not os.path.exists(blast_file):
            print(f"\n[!] BLAST results not found: {blast_file}")
            print("[-->] Use demo mode: python project.py --demo")
            sys.exit(0)

    try:
        hits = find_resistance_genes(blast_file, min_identity=args.identity)
    except (FileNotFoundError, ValueError) as e:
        print(f"[!] Error: {e}")
        sys.exit(1)

    if hits.empty:
        print("[!] No resistance genes found above threshold.")
        sys.exit(0)

    output_path = visualize_results(hits, top_n=args.top)

    print("\n" + "=" * 55)
    print("   Analysis complete!")
    print(f"   Resistance genes found: {len(hits)}")
    print(f"   Chart saved to: {output_path}")
    print("=" * 55)


if __name__ == "__main__":
    main()
