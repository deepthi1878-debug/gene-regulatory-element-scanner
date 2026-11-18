# Gene Regulatory Element Scanner – Promoter Detection

A beginner-friendly BioPython course project that accepts **DNA and protein sequences**.

- **DNA:** calculates sequence length and GC content, scans both strands for selected promoter-associated motifs (TATA box, CAAT box and GC box), and groups nearby hits into candidate promoter-like regions.
- **Protein:** calculates sequence length and amino-acid composition, plus estimated molecular weight, isoelectric point (pI) and GRAVY score when the sequence contains only standard amino acids.
- **Important:** promoter screening is performed only on DNA. Protein sequences do not contain the original upstream DNA regulatory region.

## Scientific limitation

This is an educational, rule-based motif scanner, not a validated promoter prediction model. Motifs can occur by chance, real promoters may not contain these exact motifs, and promoter activity depends on biological context. The candidate score is a simple motif-density heuristic, not a probability or confidence score.

## Run locally

Requires Python 3.9 or newer.

```bash
python -m venv .venv
```

Activate the environment:

- Windows: `.venv\\Scripts\\activate`
- macOS/Linux: `source .venv/bin/activate`

Install and launch:

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Use the app

1. Select **DNA** or **Protein**.
2. Paste a raw sequence / FASTA text, or upload a FASTA file.
3. Click **Analyze sequence**.
4. Review the results and export them as CSV or JSON.

Included examples are in `examples/`. The DNA example is synthetic and is not an experimentally validated promoter.

## Deploy online

This project uses Streamlit, so GitHub stores the code and **Streamlit Community Cloud** runs the Python app:

1. Push this repository to GitHub.
2. Open [Streamlit Community Cloud](https://share.streamlit.io/).
3. Sign in with GitHub and choose **Create app**.
4. Select this repository, branch `main`, and main file path `app.py`.
5. Deploy. Streamlit will install packages from `requirements.txt`.

## Project structure

```text
app.py
scanner.py
requirements.txt
README.md
.gitignore
examples/
  sample_promoter_dna.fasta
  sample_protein.fasta
```

## Concepts demonstrated

Python strings, functions, loops, conditionals, dictionaries, file handling, regular expressions, FASTA parsing with BioPython, reverse complements, sequence validation, and basic protein sequence analysis.
