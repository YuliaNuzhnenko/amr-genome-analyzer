# -*- coding: utf-8 -*-
"""
VetAMR Scan — Streamlit Web Interface
Wraps the core AMR analysis functions from project.py
into a professional, interactive web application.
"""

import streamlit as st
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io
import os
from dotenv import load_dotenv

# ── Load environment ─────────────────────────────────────────────────────────
load_dotenv()

# ── Page config (MUST be first Streamlit call) ───────────────────────────────
st.set_page_config(
    page_title="VetAMR Scan",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Import core functions from project.py ────────────────────────────────────
from project import generate_demo_data, find_resistance_genes, parse_subject

# ── Custom CSS ───────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', sans-serif;
}

/* Dark gradient header */
.main-header {
    background: linear-gradient(135deg, #0f2027, #203a43, #2c5364);
    padding: 2rem 2.5rem;
    border-radius: 16px;
    margin-bottom: 2rem;
    color: white;
}
.main-header h1 {
    font-size: 2.4rem;
    font-weight: 700;
    margin: 0;
    letter-spacing: -0.5px;
}
.main-header p {
    font-size: 1rem;
    opacity: 0.8;
    margin: 0.4rem 0 0 0;
}
.badge {
    display: inline-block;
    background: rgba(255,255,255,0.15);
    border: 1px solid rgba(255,255,255,0.25);
    border-radius: 20px;
    padding: 3px 12px;
    font-size: 0.75rem;
    margin-top: 0.6rem;
    margin-right: 0.4rem;
}

/* Metric cards */
.metric-card {
    background: #f8fafc;
    border: 1px solid #e2e8f0;
    border-radius: 12px;
    padding: 1.2rem 1.5rem;
    text-align: center;
}
.metric-card .value {
    font-size: 2rem;
    font-weight: 700;
    color: #1e293b;
}
.metric-card .label {
    font-size: 0.8rem;
    color: #64748b;
    text-transform: uppercase;
    letter-spacing: 0.05em;
}

/* Gene table styling */
.gene-table { border-radius: 8px; overflow: hidden; }

/* Info box */
.info-box {
    background: #eff6ff;
    border-left: 4px solid #3b82f6;
    padding: 0.8rem 1rem;
    border-radius: 0 8px 8px 0;
    font-size: 0.9rem;
    color: #1e40af;
    margin-bottom: 1rem;
}
</style>
""", unsafe_allow_html=True)

# ── Header ───────────────────────────────────────────────────────────────────
st.markdown("""
<div class="main-header">
    <h1>🦠 VetAMR Scan</h1>
    <p>Automated Antimicrobial Resistance Gene Profiler for Bacterial Genomes</p>
    <span class="badge">NCBI Entrez</span>
    <span class="badge">BLAST + CARD</span>
    <span class="badge">CS50P Project</span>
    <span class="badge">v1.0</span>
</div>
""", unsafe_allow_html=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/0/07/US-NLM-NCBI-Logo.svg/200px-US-NLM-NCBI-Logo.svg.png", width=120)
    st.markdown("---")
    st.markdown("### ⚙️ Analysis Settings")

    st.markdown("""
    <div style="background-color: #f1f5f9; padding: 15px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid #3b82f6;">
        <h4 style="margin-top: 0; color: #0f172a; font-size: 1.1rem;">🤝 Contact & Collaboration</h4>
        <p style="font-size: 0.85rem; color: #334155; line-height: 1.4;">
            This is a demo version. For detailed real-world analysis, custom pipelines, or collaboration inquiries, please reach out to me!
        </p>
        <a href="mailto:yuliannanuzhnenko@gmail.com" target="_blank" style="text-decoration: none;">
            <div style="background-color: #2563eb; color: white; padding: 8px 10px; border-radius: 6px; text-align: center; margin-bottom: 8px; font-weight: 600; font-size: 0.9rem; transition: 0.2s;">
                📧 Email Me
            </div>
        </a>
        <a href="https://github.com/YuliaNuzhnenko" target="_blank" style="text-decoration: none;">
            <div style="background-color: #1e293b; color: white; padding: 8px 10px; border-radius: 6px; text-align: center; margin-bottom: 8px; font-weight: 600; font-size: 0.9rem; transition: 0.2s;">
                🐙 GitHub Portfolio
            </div>
        </a>
    </div>
    """, unsafe_allow_html=True)
    
    mode = "Demo"


    st.markdown("---")
    min_identity = st.slider(
        "Min Identity Threshold (%)",
        min_value=50.0,
        max_value=100.0,
        value=70.0,
        step=1.0,
        help="Filter genes with identity below this threshold"
    )
    top_n = st.slider(
        "Top N Genes to Display",
        min_value=3,
        max_value=20,
        value=10,
        step=1
    )

    st.markdown("---")
    run_button = st.button("🚀 Run Analysis", use_container_width=True, type="primary")

    st.markdown("---")
    st.markdown("""
    <small style='color:#94a3b8'>
    Built with Python + Streamlit<br>
    Prototype of <b>VetAMR Scan</b><br>
    CS50P Final Project
    </small>
    """, unsafe_allow_html=True)


# ── Color palette (same as project.py) ───────────────────────────────────────
CLASS_COLORS = {
    "beta-lactam":    "#c0392b",
    "tetracycline":   "#e67e22",
    "aminoglycoside": "#8e44ad",
    "macrolide":      "#2980b9",
    "fluoroquinolone":"#16a085",
    "glycopeptide":   "#d35400",
    "trimethoprim":   "#27ae60",
    "phenicol":       "#7f8c8d",
    "fusidic acid":   "#f39c12",
}
DEFAULT_COLOR = "#34495e"


def make_chart(hits, top_n):
    """Builds the matplotlib figure and returns it (for Streamlit display)."""
    top_hits = hits.head(top_n).copy()
    parsed = top_hits["subject"].apply(parse_subject)
    top_hits["gene_name"] = parsed.apply(lambda x: x[0])
    top_hits["ab_class"]  = parsed.apply(lambda x: x[1])
    top_hits["mechanism"] = parsed.apply(lambda x: x[2])
    top_hits["label"] = top_hits.apply(
        lambda r: f"{r['gene_name']}  ({r['ab_class']})", axis=1
    )
    colors = [
        CLASS_COLORS.get(row["ab_class"].lower(), DEFAULT_COLOR)
        for _, row in top_hits.iterrows()
    ]

    fig, ax = plt.subplots(figsize=(13, 7))
    fig.patch.set_facecolor("#f0f4f8")
    ax.set_facecolor("#f0f4f8")

    bars = ax.barh(top_hits["label"], top_hits["identity"],
                   color=colors, edgecolor="white", linewidth=1.0, height=0.6)

    for bar, val, mech in zip(bars, top_hits["identity"], top_hits["mechanism"]):
        ax.text(bar.get_width() - 1.2, bar.get_y() + bar.get_height() / 2,
                f"{val:.1f}%", va="center", ha="right",
                color="white", fontweight="bold", fontsize=10)
        if mech:
            short = mech[:28] + "..." if len(mech) > 28 else mech
            ax.text(bar.get_width() + 0.4, bar.get_y() + bar.get_height() / 2,
                    short, va="center", ha="left",
                    color="#555555", fontsize=8, style="italic")

    ax.set_xlabel("Identity to CARD reference (%)", fontsize=12, labelpad=10)
    ax.set_title(
        f"Top-{top_n} Antibiotic Resistance Genes Detected\n"
        f"Staphylococcus aureus NCTC 8325 vs. CARD Database",
        fontsize=13, fontweight="bold", pad=15
    )
    ax.set_xlim(60, 130)
    ax.invert_yaxis()
    ax.axvline(x=90, color="#cc3333", linestyle="--", linewidth=1, alpha=0.6, label="90% threshold")
    ax.axvline(x=70, color="#cc8800", linestyle="--", linewidth=1, alpha=0.6, label="70% threshold")
    ax.legend(loc="lower right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=10)
    plt.tight_layout()
    return fig, top_hits


# ── Welcome state (before analysis) ─────────────────────────────────────────
if not run_button:
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("""
        <div class="metric-card">
            <div class="value">🔬</div>
            <div class="label">NCBI Entrez Integration</div>
        </div>""", unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="metric-card">
            <div class="value">🧬</div>
            <div class="label">BLAST + CARD Database</div>
        </div>""", unsafe_allow_html=True)
    with col3:
        st.markdown("""
        <div class="metric-card">
            <div class="value">📊</div>
            <div class="label">Publication-Ready Charts</div>
        </div>""", unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("""
    <div class="info-box">
        👈 Configure settings in the sidebar and click <b>Run Analysis</b> to start.
        Use <b>Demo mode</b> for instant results without internet access.
    </div>
    """, unsafe_allow_html=True)

    with st.expander("📖 How to use this tool"):
        st.markdown("""
        1. **Demo version** — click 'Run Analysis' to see instant results with realistic simulated *S. aureus* data.
        2. **Real genome analysis** is available as a professional service. Contact me to analyze your own `.gb`, `.fasta`, or raw sequence files.
        3. Adjust **identity threshold** and **top N genes** to see how the chart dynamically updates.

        **Antibiotic class color legend:**
        """)
        cols = st.columns(3)
        items = list(CLASS_COLORS.items())
        for i, (cls, color) in enumerate(items):
            with cols[i % 3]:
                st.markdown(f"<span style='background:{color};color:white;padding:2px 10px;border-radius:4px;font-size:0.85rem'>{cls}</span>", unsafe_allow_html=True)

# ── Analysis run ─────────────────────────────────────────────────────────────
else:
    blast_file = os.path.join("results", "blast_hits.tsv")

    with st.spinner("Running AMR analysis..."):
        try:
            if "Demo" in mode:
                generate_demo_data(output_dir="results")
            else:
                # In real mode, require BLAST file to already exist
                from Bio import Entrez
                if ncbi_email:
                    Entrez.email = ncbi_email
                if not os.path.exists(blast_file):
                    st.error(
                        "BLAST results file not found at `results/blast_hits.tsv`.\n\n"
                        "Please run BLAST locally and place the `-outfmt 6` output file there, "
                        "or switch to **Demo mode** for instant results."
                    )
                    st.stop()

            hits = find_resistance_genes(blast_file, min_identity=min_identity)

        except FileNotFoundError as e:
            st.error(f"File error: {e}")
            st.stop()
        except ValueError as e:
            st.error(f"Data error: {e}")
            st.stop()
        except Exception as e:
            st.error(f"Unexpected error: {e}")
            st.stop()

    if hits.empty:
        st.warning(f"No resistance genes found above {min_identity}% identity threshold. Try lowering the threshold.")
        st.stop()

    # ── Metrics row ──────────────────────────────────────────────────────────
    st.markdown("### Results")
    m1, m2, m3, m4 = st.columns(4)
    parsed_all = hits["subject"].apply(parse_subject)
    unique_classes = parsed_all.apply(lambda x: x[1]).nunique()
    top_gene_name  = parsed_all.iloc[0][0] if len(hits) > 0 else "—"
    top_identity   = hits.iloc[0]["identity"] if len(hits) > 0 else 0

    with m1:
        st.metric("Genes Found", len(hits))
    with m2:
        st.metric("Antibiotic Classes", unique_classes)
    with m3:
        st.metric("Top Gene", top_gene_name)
    with m4:
        st.metric("Top Identity", f"{top_identity:.1f}%")

    st.markdown("---")

    # ── Chart ────────────────────────────────────────────────────────────────
    fig, top_hits = make_chart(hits, top_n)
    st.pyplot(fig, use_container_width=True)

    # ── Download chart ───────────────────────────────────────────────────────
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    buf.seek(0)
    st.download_button(
        label="⬇️ Download Chart (PNG)",
        data=buf,
        file_name="amr_resistance_genes.png",
        mime="image/png",
    )

    st.markdown("---")

    # ── Gene table ───────────────────────────────────────────────────────────
    st.markdown("### Top Resistance Genes — Detail Table")

    display_df = top_hits[["gene_name", "ab_class", "mechanism", "identity", "length", "evalue", "bitscore"]].copy()
    display_df.columns = ["Gene", "Antibiotic Class", "Resistance Mechanism", "Identity (%)", "Alignment Length", "E-value", "Bitscore"]
    display_df = display_df.reset_index(drop=True)
    display_df.index += 1

    st.dataframe(
        display_df.style.background_gradient(subset=["Identity (%)"], cmap="RdYlGn"),
        use_container_width=True,
        height=400
    )

    # ── Download CSV ─────────────────────────────────────────────────────────
    csv_data = display_df.to_csv(index=True).encode("utf-8")
    st.download_button(
        label="⬇️ Download Table (CSV)",
        data=csv_data,
        file_name="amr_resistance_genes.csv",
        mime="text/csv",
    )

    st.markdown("---")

    # ── Scientific summary ───────────────────────────────────────────────────
    st.markdown("### Clinical Interpretation Summary")
    high   = hits[parsed_all.apply(lambda x: x[1]).isin(["beta-lactam", "glycopeptide"])]
    medium = hits[parsed_all.apply(lambda x: x[1]).isin(["aminoglycoside", "fluoroquinolone"])]

    if not high.empty:
        st.error(
            f"🔴 **High-priority resistance detected:** "
            f"beta-lactam and/or glycopeptide genes found ({len(high)} hits). "
            f"Potential MRSA or VRE profile. Consult a specialist."
        )
    if not medium.empty:
        st.warning(
            f"🟡 **Moderate resistance:** aminoglycoside or fluoroquinolone genes detected ({len(medium)} hits). "
            f"Review treatment options carefully."
        )
    if high.empty and medium.empty:
        st.success("🟢 No high-priority resistance genes detected above threshold.")

    plt.close(fig)
