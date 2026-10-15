import streamlit as st
import pandas as pd
from engine import GenomeFirewall
from predictor import generate_susceptibility_report

# Page Configuration
st.set_page_config(
    page_title="Genome Firewall — Clinical Engine",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Minimalist Black Terminal CSS + Interactive Hover DNA Helix
st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;700&family=Inter:wght@300;400;600;700&display=swap');

    /* Pure Black Theme Overrides */
    .stApp {
        background-color: #000000;
        color: #E5E7EB;
        font-family: 'Inter', -apple-system, sans-serif;
    }

    /* Monospace for technical data */
    code, pre, .stDataFrame {
        font-family: 'JetBrains Mono', monospace !important;
    }

    /* Table text wrapping */
    .stDataFrame table { white-space: normal !important; }
    .stDataFrame td { word-wrap: break-word !important; }

    /* Minimalist Dark Border Cards */
    .card-metric {
        background: #0A0A0A;
        border: 1px solid #1F2937;
        border-radius: 4px;
        padding: 14px 18px;
    }
    
    .metric-hdr {
        color: #6B7280;
        font-size: 11px;
        text-transform: uppercase;
        letter-spacing: 1px;
        font-weight: 600;
        font-family: 'JetBrains Mono', monospace;
    }
    
    .metric-num {
        font-size: 24px;
        font-weight: 700;
        font-family: 'JetBrains Mono', monospace;
        margin-top: 4px;
    }

    /* INTERACTIVE 3D DNA STRAND WITH HOVER LABELS */
    .dna-stage {
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        height: 280px;
        background: #000000;
        margin: 20px 0;
    }

    .dna-legend {
        display: flex;
        gap: 24px;
        margin-bottom: 20px;
        font-family: 'JetBrains Mono', monospace;
        font-size: 11px;
        text-transform: uppercase;
        letter-spacing: 1px;
    }

    .legend-a { color: #06B6D4; }
    .legend-b { color: #A855F7; }

    .dna-viewport {
        perspective: 1000px;
        width: 200px;
        height: 180px;
    }

    .dna-helix {
        position: relative;
        width: 100%;
        height: 100%;
        transform-style: preserve-3d;
        animation: spinHelix 10s linear infinite;
    }

    @keyframes spinHelix {
        0% { transform: rotateY(0deg) rotateX(15deg); }
        100% { transform: rotateY(360deg) rotateX(15deg); }
    }

    .dna-rung {
        position: absolute;
        width: 160px;
        height: 1px;
        left: 20px;
        background: rgba(255, 255, 255, 0.05);
        transform-style: preserve-3d;
    }

    /* Node Strand A - Cyan */
    .node-a {
        position: absolute;
        left: -6px;
        top: -6px;
        width: 13px;
        height: 13px;
        border-radius: 50%;
        background: #06B6D4;
        box-shadow: 0 0 8px rgba(6, 182, 212, 0.6);
        cursor: pointer;
    }

    /* Node Strand B - Violet */
    .node-b {
        position: absolute;
        right: -6px;
        top: -6px;
        width: 13px;
        height: 13px;
        border-radius: 50%;
        background: #A855F7;
        box-shadow: 0 0 8px rgba(168, 85, 247, 0.6);
        cursor: pointer;
    }

    /* Tooltip styling on hover */
    .node-a:hover::after {
        content: "Strand 1: Sense (5' -> 3')";
        position: absolute;
        left: 20px;
        top: -20px;
        background: #000000;
        border: 1px solid #06B6D4;
        color: #06B6D4;
        padding: 3px 8px;
        font-size: 10px;
        font-family: 'JetBrains Mono', monospace;
        white-space: nowrap;
        border-radius: 2px;
        z-index: 100;
    }

    .node-b:hover::after {
        content: "Strand 2: Antisense (3' -> 5')";
        position: absolute;
        right: 20px;
        top: -20px;
        background: #000000;
        border: 1px solid #A855F7;
        color: #A855F7;
        padding: 3px 8px;
        font-size: 10px;
        font-family: 'JetBrains Mono', monospace;
        white-space: nowrap;
        border-radius: 2px;
        z-index: 100;
    }

    /* Layer positioning for helical twist */
    .r1  { top: 0px;   transform: rotateY(0deg); }
    .r2  { top: 18px;  transform: rotateY(36deg); }
    .r3  { top: 36px;  transform: rotateY(72deg); }
    .r4  { top: 54px;  transform: rotateY(108deg); }
    .r5  { top: 72px;  transform: rotateY(144deg); }
    .r6  { top: 90px;  transform: rotateY(180deg); }
    .r7  { top: 108px; transform: rotateY(216deg); }
    .r8  { top: 126px; transform: rotateY(252deg); }
    .r9  { top: 144px; transform: rotateY(288deg); }
    .r10 { top: 162px; transform: rotateY(324deg); }
    </style>
""", unsafe_allow_html=True)

# Minimal Header Bar
hdr_left, hdr_right = st.columns([3, 1])
with hdr_left:
    st.markdown("<h3 style='margin:0; padding:0; font-weight:700; letter-spacing:1px; font-family:monospace;'>GENOME FIREWALL</h3>", unsafe_allow_html=True)
    st.markdown("<p style='color:#4B5563; font-size:12px; font-family:monospace; margin:0;'>ANTIMICROBIAL RESISTANCE INFERENCE SYSTEM // VERSION 3.3</p>", unsafe_allow_html=True)

with hdr_right:
    st.markdown("<div style='text-align:right; font-family:monospace; color:#10B981; font-size:11px; padding-top:6px;'>[STATUS: ONLINE]</div>", unsafe_allow_html=True)

st.divider()

# Upload Control
uploaded_file = st.file_uploader("Upload FASTA Assembly (.fasta, .fa, .fna)", type=["fasta", "fa", "fna"])

if not uploaded_file:
    # ---------------------------------------------------------
    # INITIAL STATE: LABELED 3D ROTATING DNA HELIX
    # ---------------------------------------------------------
    st.markdown("""
        <div class="dna-stage">
            <div class="dna-legend">
                <span class="legend-a">• Strand 1: Sense (5' → 3')</span>
                <span class="legend-b">• Strand 2: Antisense (3' → 5')</span>
            </div>
            <div class="dna-viewport">
                <div class="dna-helix">
                    <div class="dna-rung r1"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r2"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r3"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r4"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r5"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r6"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r7"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r8"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r9"><div class="node-a"></div><div class="node-b"></div></div>
                    <div class="dna-rung r10"><div class="node-a"></div><div class="node-b"></div></div>
                </div>
            </div>
            <div style="color: #4B5563; font-size: 11px; font-family: monospace; margin-top: 20px;">
                HOVER OVER NODES TO INSPECT STRANDS // AWAITING FASTA INPUT...
            </div>
        </div>
    """, unsafe_allow_html=True)

else:
    # ---------------------------------------------------------
    # ACTIVE DASHBOARD
    # ---------------------------------------------------------
    firewall = GenomeFirewall()
    analysis = firewall.analyze_genome(uploaded_file)
    report = generate_susceptibility_report(analysis)
    df = pd.DataFrame(report)
    
    # Quantitative Status Indicators
    total_screened = len(df)
    resistant_count = len(df[df["Prediction"].str.contains("Resistant")])
    susceptible_count = len(df[df["Prediction"].str.contains("Susceptible")])
    intrinsic_count = len(df[df["Prediction"].str.contains("Intrinsic")])

    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.markdown(f'<div class="card-metric"><div class="metric-hdr">SCREENED PANEL</div><div class="metric-num" style="color:#06B6D4;">{total_screened}</div></div>', unsafe_allow_html=True)
    with m2:
        st.markdown(f'<div class="card-metric"><div class="metric-hdr">PREDICTED RESISTANT</div><div class="metric-num" style="color:#EF4444;">{resistant_count}</div></div>', unsafe_allow_html=True)
    with m3:
        st.markdown(f'<div class="card-metric"><div class="metric-hdr">PREDICTED SUSCEPTIBLE</div><div class="metric-num" style="color:#10B981;">{susceptible_count}</div></div>', unsafe_allow_html=True)
    with m4:
        st.markdown(f'<div class="card-metric"><div class="metric-hdr">NON-SPECTRUM BARRIER</div><div class="metric-num" style="color:#F59E0B;">{intrinsic_count}</div></div>', unsafe_allow_html=True)

    st.write("")
    
    # Sequence Metadata Summary
    meta_1, meta_2, meta_3, meta_4 = st.columns(4)
    meta_1.caption(f"ACCESSION: `{analysis['accession_id']}`")
    meta_2.caption(f"SPECIES: `{analysis['species']}`")
    meta_3.caption(f"LENGTH: `{analysis['total_bp']:,} BP`")
    meta_4.caption(f"TYPE: `{analysis['assembly_type']}`")

    st.divider()

    # Section 1: Detected Resistance Genes
    st.markdown("##### DETECTED RESISTANCE GENE LOCI")
    if analysis["detected_genes"]:
        gene_df = pd.DataFrame(analysis["detected_genes"])[["name", "drug_class", "start_pos", "end_pos", "identity_pct", "coverage_pct", "ref_db"]]
        gene_df.columns = ["Gene Marker", "Target Class", "Start Locus (bp)", "End Locus (bp)", "Identity %", "Coverage %", "Ref DB"]
        st.dataframe(gene_df, width="stretch")
    else:
        st.info("No primary resistance determinants detected in reference catalog.")

    st.write("")

    # Section 2: Clinical Phenotype Matrix
    st.markdown("##### ANTIMICROBIAL SUSCEPTIBILITY MATRIX")
    matrix_df = df[["Antibiotic", "Prediction", "Confidence", "Primary Determinant", "Mechanism Explanation"]]
    st.dataframe(matrix_df, width="stretch")

    # Synchronized CSV Export
    csv_df = df[["Antibiotic", "Prediction", "Confidence", "Primary Determinant", "Mechanism Explanation", "Alternatives_JSON"]]
    csv_data = csv_df.to_csv(index=False).encode("utf-8")
    
    st.download_button(
        label="EXPORT SYNCHRONIZED REPORT (CSV)",
        data=csv_data,
        file_name=f"{analysis['accession_id']}_report.csv",
        mime="text/csv",
    )

    # Section 3: Verified Rescue Candidates
    resistant_rows = [r for r in report if "Resistant" in r["Prediction"]]
    if resistant_rows:
        st.divider()
        st.markdown("##### EVALUATED RESCUE ALTERNATIVES")
        
        for r in resistant_rows:
            if r["Alternatives"]:
                st.markdown(f"**Target Alternatives for {r['Antibiotic']}** (*Primary Cause: `{r['Primary Determinant']}`*)")
                a_cols = st.columns(len(r["Alternatives"]))
                for idx, alt in enumerate(r["Alternatives"]):
                    with a_cols[idx]:
                        st.success(f"**Option {idx+1}: {alt['drug']}**\n\n• **Class:** {alt['class']}\n\n• **Rationale:** {alt['mechanism']}")
            else:
                st.error(f"**{r['Antibiotic']} (`{r['Primary Determinant']}`):** No viable rescue alternative identified in standard panel due to co-occurring resistance mechanisms.")