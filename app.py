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
@import url('https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Outfit', sans-serif !important;
}

/* Glassmorphism Header */
.main-header {
    background: rgba(15, 23, 42, 0.85);
    backdrop-filter: blur(12px);
    -webkit-backdrop-filter: blur(12px);
    border: 1px solid rgba(255, 255, 255, 0.1);
    padding: 2.5rem 3rem;
    border-radius: 24px;
    margin-bottom: 2.5rem;
    color: white;
    box-shadow: 0 20px 40px rgba(0,0,0,0.15);
    position: relative;
    overflow: hidden;
}
.main-header::before {
    content: '';
    position: absolute;
    top: -50%; left: -50%;
    width: 200%; height: 200%;
    background: radial-gradient(circle, rgba(56,189,248,0.1) 0%, transparent 60%);
    z-index: 0;
}
.main-header > * { position: relative; z-index: 1; }

.main-header h1 {
    font-size: 3.2rem;
    font-weight: 700;
    margin: 0;
    background: linear-gradient(90deg, #38bdf8, #818cf8);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    letter-spacing: -1px;
}
.main-header p {
    font-size: 1.15rem;
    color: #94a3b8;
    font-weight: 400;
    margin: 0.5rem 0 1rem 0;
}
.badge {
    display: inline-block;
    background: rgba(56, 189, 248, 0.1);
    border: 1px solid rgba(56, 189, 248, 0.2);
    border-radius: 20px;
    padding: 4px 14px;
    font-size: 0.75rem;
    font-weight: 500;
    margin-top: 0.6rem;
    margin-right: 0.5rem;
    color: #38bdf8;
    backdrop-filter: blur(4px);
}

/* Premium Metric Cards */
.metric-card {
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 16px;
    padding: 1.5rem;
    text-align: center;
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05), 0 2px 4px -1px rgba(0, 0, 0, 0.03);
    transition: all 0.3s ease;
}
.metric-card:hover {
    transform: translateY(-5px);
    box-shadow: 0 12px 20px -5px rgba(0, 0, 0, 0.1);
    border-color: #cbd5e1;
}
.metric-card .value {
    font-size: 2.5rem;
    line-height: 1.2;
}
.metric-card .label {
    font-size: 0.85rem;
    font-weight: 600;
    color: #64748b;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    margin-top: 0.5rem;
}

/* Info box upgrade */
.info-box {
    background: linear-gradient(to right, #eff6ff, #ffffff);
    border-left: 4px solid #3b82f6;
    padding: 1rem 1.2rem;
    border-radius: 8px;
    font-size: 0.95rem;
    font-weight: 500;
    color: #1e3a8a;
    margin-bottom: 1.5rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
}

/* Sidebar button hover */
.sidebar-btn {
    transition: transform 0.2s, box-shadow 0.2s !important;
}
.sidebar-btn:hover {
    transform: translateY(-2px);
    box-shadow: 0 6px 15px rgba(59, 130, 246, 0.25);
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

    mode = st.radio(
        "Mode",
        ["🎭 Demo (simulated data)", "🔬 Real genome (NCBI)"],
        help="Demo mode uses simulated data. Real mode requires professional setup."
    )

    st.markdown("---")

    if "Real genome" in mode:
        accession = st.text_input(
            "NCBI Accession Number",
            value="NC_007795.1",
            help="Example: NC_007795.1 (S. aureus NCTC 8325)"
        )
        ncbi_email = st.text_input(
            "Your Email (for NCBI)",
            value="",
            help="NCBI requires a valid email for API access",
            type="default"
        )
    else:
        accession = "NC_007795.1"
        ncbi_email = ""


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
                hits = find_resistance_genes(blast_file, min_identity=min_identity)
            else:
                st.markdown("""
                <div style="background: linear-gradient(145deg, #ffffff, #f8fafc); padding: 40px; border-radius: 20px; margin-top: 20px; border: 1px solid #e2e8f0; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.1); text-align: center;">
                    <div style="font-size: 3rem; margin-bottom: 10px;">🔬</div>
                    <h2 style="color: #0f172a; font-weight: 700; margin-bottom: 15px; font-family: 'Outfit', sans-serif;">Need Real Genomic Analysis?</h2>
                    <p style="font-size: 1.1rem; color: #475569; line-height: 1.6; max-width: 650px; margin: 0 auto 25px auto;">
                        This web app is a <b>limited demonstration</b> to showcase my bioinformatics capabilities.<br>
                        If you need professional AMR profiling for real genomes, raw sequence processing, or custom bioinformatics pipelines, I provide these services on a freelance basis.
                    </p>
                    <p style="font-size: 1.2rem; color: #0f172a; margin-bottom: 30px; font-weight: 700;">
                        Hire me to analyze your data!
                    </p>
                    <div style="display: flex; justify-content: center; gap: 15px; flex-wrap: wrap;">
                        <a href="mailto:yuliannanuzhnenko@gmail.com" target="_blank" style="text-decoration: none;">
                            <div class="sidebar-btn" style="background: linear-gradient(to right, #2563eb, #3b82f6); color: white; padding: 12px 24px; border-radius: 8px; font-weight: 600; font-size: 1rem;">
                                📧 Contact via Email
                            </div>
                        </a>
                        <a href="https://github.com/YuliaNuzhnenko" target="_blank" style="text-decoration: none;">
                            <div class="sidebar-btn" style="background: #0f172a; color: white; padding: 12px 24px; border-radius: 8px; font-weight: 600; font-size: 1rem;">
                                🐙 View GitHub Profile
                            </div>
                        </a>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                st.stop()

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
