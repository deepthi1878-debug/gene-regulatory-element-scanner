import json
from io import StringIO
import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scanner import analyze_record, clean_sequence

st.set_page_config(page_title="Gene Regulatory Element Scanner", page_icon="🧬", layout="wide")
st.title("🧬 Gene Regulatory Element Scanner")
st.caption("Promoter-associated DNA motif screening and basic protein sequence analysis")
st.info("Promoters are DNA regulatory regions. Choose DNA for motif screening. Protein sequences are accepted for protein-property analysis only.")

sequence_type = st.radio("Sequence type", ["DNA", "Protein"], horizontal=True)
input_mode = st.radio("Input method", ["Paste sequence", "Upload FASTA"], horizontal=True)
raw_text = ""
if input_mode == "Paste sequence":
    raw_text = st.text_area("Paste a raw sequence or FASTA-formatted sequence", height=170, placeholder=">sample_sequence\nATGCGTATAAAGCGTACG...")
else:
    uploaded = st.file_uploader("Choose a FASTA file", type=["fa", "fasta", "fna", "faa", "txt"])
    if uploaded is not None:
        raw_text = uploaded.getvalue().decode("utf-8", errors="replace")

if st.button("Analyze sequence", type="primary"):
    if not raw_text.strip():
        st.error("Please paste a sequence or upload a FASTA file.")
        st.stop()
    try:
        if raw_text.lstrip().startswith(">"):
            records = list(SeqIO.parse(StringIO(raw_text), "fasta"))
            if not records:
                st.error("No valid FASTA records were found.")
                st.stop()
        else:
            records = [SeqRecord(Seq(clean_sequence(raw_text)), id="input_sequence", description="User supplied sequence")]
        results = [analyze_record(record, sequence_type) for record in records]
        st.success(f"Successfully analyzed {len(results)} sequence record(s).")
        for result in results:
            st.header(result["record_id"])
            st.caption(result["description"])
            if sequence_type == "DNA":
                c1, c2, c3 = st.columns(3)
                c1.metric("Sequence length", f'{result["length"]:,} bp')
                c2.metric("GC content", f'{result["gc_percent"]}%')
                c3.metric("Motif hits", len(result["motif_hits"]))
                st.subheader("Promoter-associated motif hits")
                if result["motif_hits"]:
                    st.dataframe(pd.DataFrame(result["motif_hits"]), use_container_width=True, hide_index=True)
                else:
                    st.warning("No selected motifs were found. This does not prove that the sequence has no promoter.")
                st.subheader("Promoter candidate regions")
                if result["candidate_regions"]:
                    st.dataframe(pd.DataFrame(result["candidate_regions"]), use_container_width=True, hide_index=True)
                else:
                    st.write("No candidate regions were generated because no motifs were detected.")
                st.caption("Coordinates are 1-based and inclusive. Both strands are scanned. The score is a motif-density heuristic, not a probability.")
            else:
                c1, c2 = st.columns(2)
                c1.metric("Protein length", f'{result["length"]:,} amino acids')
                mw = result.get("molecular_weight_da")
                c2.metric("Estimated molecular weight", f"{mw:,.2f} Da" if mw is not None else "Not calculated")
                if result.get("isoelectric_point") is not None:
                    st.write(f'**Estimated isoelectric point (pI):** {result["isoelectric_point"]}')
                    st.write(f'**GRAVY score:** {result["gravy"]}')
                st.subheader("Amino-acid composition")
                st.dataframe(pd.DataFrame([{"Amino acid": aa, "Count": count} for aa, count in result["amino_acid_counts"].items()]), use_container_width=True, hide_index=True)
                st.caption(result["note"])
                st.warning("Promoter detection was not performed: promoters are DNA regulatory regions, not amino-acid sequences.")
        exportable = [{
            "record_id": r["record_id"], "sequence_type": r["sequence_type"], "length": r["length"],
            "gc_percent": r.get("gc_percent", ""),
            "motif_hits": json.dumps(r.get("motif_hits", [])),
            "candidate_regions": json.dumps(r.get("candidate_regions", [])),
            "protein_properties": json.dumps({k: r.get(k) for k in ["molecular_weight_da", "isoelectric_point", "gravy", "amino_acid_counts"] if k in r}),
        } for r in results]
        st.download_button("Download results as CSV", pd.DataFrame(exportable).to_csv(index=False), file_name="scanner_results.csv", mime="text/csv")
        st.download_button("Download detailed results as JSON", json.dumps(results, indent=2), file_name="scanner_results.json", mime="application/json")
    except Exception as exc:
        st.error(f"Could not analyze sequence: {exc}")
        st.caption("Check that the selected sequence type matches the input alphabet and that the FASTA file is correctly formatted.")

with st.expander("Which motifs are scanned?"):
    st.markdown("- **TATA box:** simplified pattern TATA[AT]A[AT][AG]\n- **CAAT box:** CCAAT\n- **GC box:** GGGCGG\n\nThese are simplified eukaryotic promoter-associated motifs. Motif presence alone cannot confirm promoter function.")
