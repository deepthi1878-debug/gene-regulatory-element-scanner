"""Educational promoter-motif scanner and basic protein sequence analyzer."""
from collections import Counter
import re
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

DNA_ALPHABET = set("ACGTN")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
PROMOTER_MOTIFS = {"TATA box": r"TATA[AT]A[AT][AG]", "CAAT box": r"CCAAT", "GC box": r"GGGCGG"}

def clean_sequence(raw):
    lines = []
    for line in raw.splitlines():
        line = line.strip()
        if line and not line.startswith(">"):
            lines.append(re.sub(r"\s+", "", line))
    return "".join(lines).upper()

def validate_sequence(sequence, sequence_type):
    if not sequence:
        return False, "No sequence was provided."
    allowed = DNA_ALPHABET if sequence_type == "DNA" else PROTEIN_ALPHABET
    if sequence_type not in ("DNA", "Protein"):
        return False, "Choose DNA or Protein as the sequence type."
    invalid = sorted(set(sequence) - allowed)
    if invalid:
        return False, "Invalid sequence character(s): " + ", ".join(invalid)
    if sequence_type == "Protein" and "*" in sequence[:-1]:
        return False, "A stop symbol (*) is only allowed at the end of a protein sequence."
    return True, "Valid sequence."

def find_motifs(sequence):
    seq = str(sequence).upper()
    hits = []
    for strand, strand_seq in [("+", seq), ("-", str(Seq(seq).reverse_complement()))]:
        for motif_name, pattern in PROMOTER_MOTIFS.items():
            for match in re.finditer(r"(?=(" + pattern + r"))", strand_seq):
                start0, end0 = match.start(), match.start() + len(match.group(1))
                if strand == "+":
                    start, end = start0 + 1, end0
                else:
                    start, end = len(seq) - end0 + 1, len(seq) - start0
                hits.append({"motif": motif_name, "start": start, "end": end,
                             "strand": strand, "matched_sequence": match.group(1)})
    return sorted(hits, key=lambda h: (h["start"], h["end"], h["strand"], h["motif"]))

def group_candidate_regions(hits, max_gap=100):
    if not hits:
        return []
    groups, current = [], [hits[0]]
    current_end = hits[0]["end"]
    for hit in sorted(hits[1:], key=lambda h: h["start"]):
        if hit["start"] - current_end <= max_gap:
            current.append(hit)
            current_end = max(current_end, hit["end"])
        else:
            groups.append(current)
            current, current_end = [hit], hit["end"]
    groups.append(current)
    regions = []
    for index, group in enumerate(groups, 1):
        names = sorted({h["motif"] for h in group})
        regions.append({
            "region": f"Candidate {index}",
            "start": min(h["start"] for h in group),
            "end": max(h["end"] for h in group),
            "motif_count": len(group),
            "motifs_found": ", ".join(names),
            "score": min(100, len(group) * 25),
            "interpretation": "Higher motif density; investigate as a promoter candidate." if len(group) >= 2 else "Single motif only; weak evidence, needs further analysis.",
        })
    return regions

def analyze_dna(sequence):
    seq = str(sequence).upper()
    length = len(seq)
    hits = find_motifs(seq)
    gc_count = seq.count("G") + seq.count("C")
    return {
        "length": length,
        "gc_percent": round(gc_count / length * 100, 2) if length else 0.0,
        "a_count": seq.count("A"), "t_count": seq.count("T"),
        "g_count": seq.count("G"), "c_count": seq.count("C"),
        "motif_hits": hits, "candidate_regions": group_candidate_regions(hits),
    }

def analyze_protein(sequence):
    seq = str(sequence).upper().rstrip("*")
    if not seq:
        raise ValueError("Protein sequence contains no amino acids.")
    standard = set("ACDEFGHIKLMNPQRSTVWY")
    result = {"length": len(seq), "amino_acid_counts": dict(sorted(Counter(seq).items()))}
    if set(seq) - standard:
        result.update({"molecular_weight_da": None, "note": "Sequence contains ambiguous/non-standard amino-acid symbols; molecular weight was not calculated."})
        return result
    analysis = ProteinAnalysis(seq)
    result.update({
        "molecular_weight_da": round(analysis.molecular_weight(), 2),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "gravy": round(analysis.gravy(), 3),
        "note": "Protein properties are sequence-based estimates.",
    })
    return result

def analyze_record(record, sequence_type):
    sequence = str(record.seq).upper()
    valid, message = validate_sequence(sequence, sequence_type)
    if not valid:
        raise ValueError(f"{record.id}: {message}")
    result = analyze_dna(sequence) if sequence_type == "DNA" else analyze_protein(sequence)
    result.update({"record_id": record.id, "description": record.description, "sequence_type": sequence_type})
    return result
