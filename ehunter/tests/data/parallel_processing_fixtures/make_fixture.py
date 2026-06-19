#!/usr/bin/env python3
"""Generate a deterministic synthetic multi-contig ExpansionHunter test fixture.

Produces, under this script's directory:
  reference.fa (+ .fai)
  variant_catalog.json
  reads.bam (+ .bai)
  reads.cram (+ .crai)

The fixture is designed to validate byte-identical EH output across thread
counts. It deliberately covers a set of edge-case scenarios (see README.md).

Catalog coordinate convention (verified against EH source):
  ReferenceRegion "chr:S-E" is 0-based half-open. getSequence() does
  faidx_fetch_seq(start, end-1), so the repeat tract occupies the reference
  0-based half-open interval [S, E), i.e. 1-based [S+1, E].

Two non-obvious requirements that EH's full graph genotyper imposes (both
learned by iterating against the binary; see README.md "How the reads are
built"):

  1. A reverse-strand read's SEQ in the BAM must be stored in FORWARD-reference
     orientation. pysam/htslib store query_sequence verbatim and use the
     is_reverse flag for strand, so we must NOT reverse-complement the bytes we
     hand to pysam (doing so double-reverses them and the graph aligner rejects
     the mate, leaving the pair unclassifiable).

  2. The genotyper's isLowDepth() filter (VariantAnalyzer.cpp) discards a locus
     whose estimated diploid depth is below --min-locus-coverage (default 10).
     With the default --region-extension-length of 1000, the graph flank nodes
     are 1000 bp each, so
         depth = readLen * readCount / (2*1000 - readLen)
     and we need enough read pairs per locus (~90) to clear depth >= 10. The
     optimized-streaming fast path does not apply this filter, but seeking and
     low-mem-streaming full genotyping do.

Everything is seeded so the output bytes are identical on every run.
"""

import json
import os
import random

import pysam

# ----------------------------------------------------------------------------
# Determinism
# ----------------------------------------------------------------------------
RNG = random.Random(20260619)  # fixed seed; never use time-based randomness

HERE = os.path.dirname(os.path.abspath(__file__))
REF_FA = os.path.join(HERE, "reference.fa")
CATALOG = os.path.join(HERE, "variant_catalog.json")
BAM = os.path.join(HERE, "reads.bam")
CRAM = os.path.join(HERE, "reads.cram")

BASES = "ACGT"

# Number of spanning read pairs per genotyped locus. Chosen so that, with the
# default 1000 bp graph flanks, the estimated diploid depth comfortably exceeds
# --min-locus-coverage (10). depth ~= readLen * pairs / (2000 - readLen), so
# shorter reads need more pairs; pairs_for_read_len() sizes each locus to clear
# depth ~= 13-15 in full genotyping. See the module docstring.
def pairs_for_read_len(read_len):
    # solve readLen * N / (2000 - readLen) >= 14  ->  N >= 14*(2000-readLen)/readLen
    import math
    return int(math.ceil(14.0 * (2000 - read_len) / read_len))


def rand_seq(n):
    return "".join(RNG.choice(BASES) for _ in range(n))


# ----------------------------------------------------------------------------
# Reference design
#
# Five contigs chr1..chr5 (this exact order; EH derives contig index from the
# FASTA/BAM header order). Each is a few kb of random ACGT with embedded perfect
# repeat tracts where the catalog loci sit. Each tract has >= 1000 bp of flank on
# both sides so the full-extension (1000 bp) graph flank nodes are well defined.
# ----------------------------------------------------------------------------

tract_index = {}  # locus_id -> (contig_name, record dict)
contigs = {}  # name -> sequence


def build_contig(plan):
    """Build a contig sequence from a plan of segments.

    plan: list of ('random', n) or ('repeat', motif, count, locus_meta_dict).
    Returns (sequence, list_of_tract_records); each tract record carries the
    0-based half-open [start, end) of the embedded repeat tract plus metadata.
    """
    seq_parts = []
    pos = 0  # 0-based offset of next base
    tracts = []
    for item in plan:
        if item[0] == "random":
            seq_parts.append(rand_seq(item[1]))
            pos += item[1]
        elif item[0] == "repeat":
            _, motif, count, meta = item
            start = pos
            seq_parts.append(motif * count)
            pos += len(motif) * count
            rec = dict(meta)
            rec.update(motif=motif, count=count, start=start, end=pos)
            tracts.append(rec)
        else:
            raise ValueError(item)
    return "".join(seq_parts), tracts


# --- chr1: simple repeat locus (CAG) + far-pair scenario (b) ----------------
# Reads on chr1 are 100 bp. Repeat tract long enough to span and flank.
contigs["chr1"], t = build_contig([
    ("random", 1500),
    ("repeat", "CAG", 15, dict(locus_id="CHR1_CAG")),  # 45 bp tract
    ("random", 4000),
])
for r in t:
    tract_index[r["locus_id"]] = ("chr1", r)

# --- chr2: simple repeat locus (CCG) + cross-contig mate scenario (a) -------
# Reads on chr2 are 150 bp (mixed-length vs chr1, scenario (f)).
contigs["chr2"], t = build_contig([
    ("random", 1200),
    ("repeat", "CCG", 12, dict(locus_id="CHR2_CCG")),  # 36 bp tract
    ("random", 3800),
])
for r in t:
    tract_index[r["locus_id"]] = ("chr2", r)

# --- chr3: locus with ZERO coverage (scenario d) ----------------------------
contigs["chr3"], t = build_contig([
    ("random", 1200),
    ("repeat", "GAA", 10, dict(locus_id="CHR3_GAA")),  # 30 bp tract, no reads
    ("random", 3000),
])
for r in t:
    tract_index[r["locus_id"]] = ("chr3", r)

# --- chr4: multi-variant locus (scenario c) ---------------------------------
# LocusStructure "(A|T)CCCC(ATTCT)*": a SmallVariant followed by a Repeat. The
# variantId string order ("ATTCT" vs the SNV) differs from positional order,
# exercising the optimized fast path's unordered_map iteration. This locus goes
# through FULL genotyping (the fast path only handles single-repeat loci).
chr4_pre = rand_seq(1400)
snv_base = "A"          # deterministic SNV reference base (A allele)
linker = "CCCC"
atct_motif = "ATTCT"
atct_count = 14         # 70 bp tract
chr4_repeat = atct_motif * atct_count
chr4_post = rand_seq(3500)
contigs["chr4"] = chr4_pre + snv_base + linker + chr4_repeat + chr4_post
snv_start = len(chr4_pre)
repeat_start = snv_start + 1 + len(linker)
repeat_end = repeat_start + len(chr4_repeat)
tract_index["CHR4_MULTI"] = (
    "chr4",
    dict(
        locus_id="CHR4_MULTI", motif=atct_motif, count=atct_count,
        start=repeat_start, end=repeat_end,
        snv_start=snv_start, snv_base=snv_base, linker=linker,
    ),
)

# --- chr5: NO catalog loci, but physical home of cross-contig mates (e) ------
contigs["chr5"], _ = build_contig([("random", 5000)])

ORDER = ["chr1", "chr2", "chr3", "chr4", "chr5"]


# ----------------------------------------------------------------------------
# Write reference.fa and index it
# ----------------------------------------------------------------------------
def write_reference():
    with open(REF_FA, "w") as f:
        for name in ORDER:
            f.write(">%s\n" % name)
            s = contigs[name]
            for i in range(0, len(s), 60):
                f.write(s[i : i + 60] + "\n")
    pysam.faidx(REF_FA)


# ----------------------------------------------------------------------------
# Catalog
# ----------------------------------------------------------------------------
def write_catalog():
    catalog = []

    r = tract_index["CHR1_CAG"][1]
    catalog.append({
        "LocusId": "CHR1_CAG", "LocusStructure": "(CAG)*",
        "ReferenceRegion": "chr1:%d-%d" % (r["start"], r["end"]),
        "VariantType": "Repeat",
    })

    r = tract_index["CHR2_CCG"][1]
    catalog.append({
        "LocusId": "CHR2_CCG", "LocusStructure": "(CCG)*",
        "ReferenceRegion": "chr2:%d-%d" % (r["start"], r["end"]),
        "VariantType": "Repeat",
    })

    r = tract_index["CHR3_GAA"][1]
    catalog.append({
        "LocusId": "CHR3_GAA", "LocusStructure": "(GAA)*",
        "ReferenceRegion": "chr3:%d-%d" % (r["start"], r["end"]),
        "VariantType": "Repeat",
    })

    r = tract_index["CHR4_MULTI"][1]
    s = r["snv_start"]
    catalog.append({
        "LocusId": "CHR4_MULTI", "LocusStructure": "(A|T)CCCC(ATTCT)*",
        "ReferenceRegion": [
            "chr4:%d-%d" % (s, s + 1),
            "chr4:%d-%d" % (r["start"], r["end"]),
        ],
        "VariantType": ["SmallVariant", "Repeat"],
    })

    with open(CATALOG, "w") as f:
        json.dump(catalog, f, indent=4)
        f.write("\n")


# ----------------------------------------------------------------------------
# Reads
#
# Read records are accumulated as dicts and emitted coordinate-sorted. Both
# mates' SEQ are stored in FORWARD-reference orientation (see module docstring);
# strand is conveyed only via the is_reverse / mate_is_reverse flags.
# ----------------------------------------------------------------------------
reads = []  # list of pair dicts


def contig_seq_with_repeat(contig, tract_start, tract_end, motif, alt_count):
    """Contig sequence with the repeat tract resized to alt_count copies."""
    s = contigs[contig]
    return s[:tract_start] + (motif * alt_count) + s[tract_end:]


def add_pair(name, tid_name, r1_start, r1_seq, r2_start, r2_seq,
             r1_reversed, r2_reversed, r2_tid=None):
    reads.append(dict(
        name=name, tid=tid_name,
        r1_start=r1_start, r1_seq=r1_seq, r1_rev=r1_reversed,
        r2_tid=r2_tid if r2_tid is not None else tid_name,
        r2_start=r2_start, r2_seq=r2_seq, r2_rev=r2_reversed,
    ))


def tile_spanning_pairs(locus_id, read_len, n_pairs, prefix):
    """Generate FR spanning read pairs over the (reference-length) repeat tract.

    read1 (forward) starts upstream and spans the tract; read2 (reverse-strand
    flag, forward SEQ bytes) sits downstream so the fragment is a proper FR pair.
    Reads match the reference exactly (reference-length allele) so EH genotypes
    the locus at the reference repeat size.
    """
    contig_name, rec = tract_index[locus_id]
    tract_start, tract_end = rec["start"], rec["end"]
    full = contig_seq_with_repeat(
        contig_name, tract_start, tract_end, rec["motif"], rec["count"])

    span_window_lo = max(0, tract_end - read_len)   # latest start still spanning
    span_window_hi = max(0, tract_start - 5)         # keep >= 5 bp left flank
    if span_window_lo > span_window_hi:
        raise ValueError("read_len %d too short to span tract of %d bp at %s"
                         % (read_len, tract_end - tract_start, locus_id))
    window = span_window_hi - span_window_lo + 1

    for i in range(n_pairs):
        r1_start = span_window_lo + (i * 7) % window if window > 1 else span_window_lo
        r1_seq = full[r1_start : r1_start + read_len]
        if len(r1_seq) < read_len:
            continue
        # Mate downstream, clear of the tract; proper FR pair.
        r2_start = r1_start + 400 - read_len
        if r2_start + read_len > len(full):
            r2_start = len(full) - read_len
        r2_seq = full[r2_start : r2_start + read_len]   # forward-ref SEQ bytes
        add_pair("%s_%03d" % (prefix, i), contig_name,
                 r1_start, r1_seq, r2_start, r2_seq, False, True)


# chr1: 100 bp reads, spanning the CAG tract -- fast path in optimized mode
tile_spanning_pairs("CHR1_CAG", read_len=100, n_pairs=pairs_for_read_len(100), prefix="chr1cag")

# chr2: 150 bp reads, spanning the CCG tract -- fast path, mixed read length
tile_spanning_pairs("CHR2_CCG", read_len=150, n_pairs=pairs_for_read_len(150), prefix="chr2ccg")

# chr4: 150 bp reads spanning the ATTCT tract -- full genotyping (multi-variant)
tile_spanning_pairs("CHR4_MULTI", read_len=150, n_pairs=pairs_for_read_len(150), prefix="chr4multi")


# ----------------------------------------------------------------------------
# Scenario (b): same-contig far pair on chr1.
# A read pair whose two ends are far apart on chr1 (mate distance >= 1000 bp):
# read1 sits ~1500 bp upstream of the CHR1_CAG locus, its mate spans the repeat.
# EH's prepass caches such distant mates within the region-extension window.
# ----------------------------------------------------------------------------
def add_far_pair_chr1():
    contig_name, rec = tract_index["CHR1_CAG"]
    tract_start, tract_end = rec["start"], rec["end"]
    full = contig_seq_with_repeat(
        contig_name, tract_start, tract_end, rec["motif"], rec["count"])
    read_len = 100
    span_start = max(0, tract_start - 10)              # mate that spans the tract
    up_start = max(0, span_start - 1500)               # read ~1500 bp upstream
    # read1 = upstream (forward), read2 = downstream spanning mate (reverse flag).
    add_pair("chr1_farpair_000", contig_name,
             up_start, full[up_start : up_start + read_len],
             span_start, full[span_start : span_start + read_len],
             False, True)


add_far_pair_chr1()


# ----------------------------------------------------------------------------
# Scenario (a): cross-contig mate. The CHR2_CCG-spanning mate is on chr2; its
# mate maps to chr5 (which has no catalog loci). Exercises mate handling where
# read and mate are on different contigs.
# ----------------------------------------------------------------------------
def add_cross_contig_pairs():
    contig_name, rec = tract_index["CHR2_CCG"]
    tract_start, tract_end = rec["start"], rec["end"]
    full = contig_seq_with_repeat(
        contig_name, tract_start, tract_end, rec["motif"], rec["count"])
    chr5 = contigs["chr5"]
    read_len = 150
    for i in range(6):
        chr2_start = max(0, tract_start - 20 + i * 3)  # chr2 mate spans the repeat
        chr5_start = 1000 + i * 200                    # mate on chr5
        add_pair("chr2_crosscontig_%03d" % i, contig_name,
                 chr2_start, full[chr2_start : chr2_start + read_len],
                 chr5_start, chr5[chr5_start : chr5_start + read_len],
                 False, True, r2_tid="chr5")


add_cross_contig_pairs()


# ----------------------------------------------------------------------------
# Emit BAM (coordinate-sorted) then CRAM
# ----------------------------------------------------------------------------
def write_bam_and_cram():
    tid_of = {name: i for i, name in enumerate(ORDER)}
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": len(contigs[name])} for name in ORDER],
    }

    aln_specs = []
    for p in reads:
        r1_tid, r2_tid = tid_of[p["tid"]], tid_of[p["r2_tid"]]
        r1_len, r2_len = len(p["r1_seq"]), len(p["r2_seq"])
        tlen = (p["r2_start"] + r2_len) - p["r1_start"] if r1_tid == r2_tid else 0
        aln_specs.append(dict(
            name=p["name"], is_read1=True, tid=r1_tid, pos=p["r1_start"],
            seq=p["r1_seq"], reverse=p["r1_rev"], mate_tid=r2_tid,
            mate_pos=p["r2_start"], mate_reverse=p["r2_rev"], tlen=tlen))
        aln_specs.append(dict(
            name=p["name"], is_read1=False, tid=r2_tid, pos=p["r2_start"],
            seq=p["r2_seq"], reverse=p["r2_rev"], mate_tid=r1_tid,
            mate_pos=p["r1_start"], mate_reverse=p["r1_rev"],
            tlen=-tlen if r1_tid == r2_tid else 0))

    aln_specs.sort(key=lambda a: (a["tid"], a["pos"], 0 if a["is_read1"] else 1, a["name"]))

    with pysam.AlignmentFile(BAM, "wb", header=header) as out:
        for a in aln_specs:
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = a["name"]
            seg.query_sequence = a["seq"]
            seg.flag = 0
            seg.is_paired = True
            seg.is_proper_pair = a["mate_tid"] == a["tid"]
            seg.is_read1 = a["is_read1"]
            seg.is_read2 = not a["is_read1"]
            seg.is_reverse = a["reverse"]
            seg.mate_is_reverse = a["mate_reverse"]
            seg.reference_id = a["tid"]
            seg.reference_start = a["pos"]
            seg.mapping_quality = 60
            seg.cigartuples = [(0, len(a["seq"]))]  # all-M
            seg.next_reference_id = a["mate_tid"]
            seg.next_reference_start = a["mate_pos"]
            seg.template_length = a["tlen"]
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(a["seq"]))
            out.write(seg)

    pysam.index(BAM)

    with pysam.AlignmentFile(BAM, "rb") as inb:
        with pysam.AlignmentFile(CRAM, "wc", template=inb, reference_filename=REF_FA) as outc:
            for seg in inb:
                outc.write(seg)
    pysam.index(CRAM)


# ----------------------------------------------------------------------------
write_reference()
write_catalog()
write_bam_and_cram()

print("Wrote:")
for p in (REF_FA, REF_FA + ".fai", CATALOG, BAM, BAM + ".bai", CRAM, CRAM + ".crai"):
    print("  %s (%d bytes)" % (p, os.path.getsize(p)))

print("\nTract placements (0-based half-open):")
for lid in ("CHR1_CAG", "CHR2_CCG", "CHR3_GAA", "CHR4_MULTI"):
    c, r = tract_index[lid]
    print("  %-12s %s:%d-%d  motif=%s x%d" % (lid, c, r["start"], r["end"], r["motif"], r["count"]))
print("  total read pairs: %d  (alignment records: %d)" % (len(reads), 2 * len(reads)))
