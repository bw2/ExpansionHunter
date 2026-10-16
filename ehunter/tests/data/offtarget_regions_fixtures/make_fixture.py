#!/usr/bin/env python3
"""Generate a deterministic synthetic ExpansionHunter fixture with off-target read pairs.

Produces, under this script's directory:
  reference.fa (+ .fai)
  variant_catalog.json
  reads.bam (+ .bai)

The fixture checks that optimized-streaming mode genotypes a locus that lists OfftargetRegions
exactly as seeking mode does. See README.md for the layout and the read scenarios.

Read-building conventions are the same as in ../parallel_processing_fixtures/make_fixture.py:
a mapped reverse-strand read's SEQ is stored in forward-reference orientation (strand is conveyed
by the is_reverse flag only), and catalog ReferenceRegions are 0-based half-open.

An unmapped read is stored as sequenced, the way BWA stores it: SEQ is the read's own sequence and
is_reverse is not set. Such a read is placed at its mapped mate's position (BWA's convention for a
pair with one unmapped mate), so region queries and the coordinate-sorted stream both reach it.

Everything is seeded so the output bytes are identical on every run.
"""

import json
import os
import random

import pysam

RNG = random.Random(20261015)  # fixed seed; never use time-based randomness

HERE = os.path.dirname(os.path.abspath(__file__))
REF_FA = os.path.join(HERE, "reference.fa")
CATALOG = os.path.join(HERE, "variant_catalog.json")
BAM = os.path.join(HERE, "reads.bam")

READ_LENGTH = 150
FRAGMENT_LENGTH = 400

# C9LIKE: a GGCCCC repeat on chr1 with 5 reference repeat units and an expanded allele of 200
# units (1200 bp, much longer than a fragment), so read pairs from inside the expansion carry no
# flank and can only be counted through the off-target regions.
C9LIKE_MOTIF = "GGCCCC"
C9LIKE_REFERENCE_UNITS = 5
C9LIKE_EXPANDED_UNITS = 200

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def random_sequence(length):
    return "".join(RNG.choice("ACGT") for _ in range(length))


# ----------------------------------------------------------------------------
# Reference
#
# chr1: C9LIKE repeat at [1500, 1530).
# chr2: an off-target region around a short GGCCCC tract at [1200, 1236), plus an ordinary CAG
#       locus at [3236, 3281) that has no off-target regions.
# chr3: an off-target region around a short GGGGCC tract at [1200, 1236). It is the reverse
#       complement of GGCCCC, so pure-repeat reads placed here are stored on the opposite strand.
# ----------------------------------------------------------------------------
CHR1_LEFT_FLANK = random_sequence(1500)
CHR1_RIGHT_FLANK = random_sequence(3000)
C9LIKE_START = len(CHR1_LEFT_FLANK)
C9LIKE_END = C9LIKE_START + C9LIKE_REFERENCE_UNITS * len(C9LIKE_MOTIF)

CHR2_CAG_UNITS = 15
contigs = {
    "chr1": CHR1_LEFT_FLANK + C9LIKE_MOTIF * C9LIKE_REFERENCE_UNITS + CHR1_RIGHT_FLANK,
    "chr2": random_sequence(1200) + "GGCCCC" * 6 + random_sequence(2000) + "CAG" * CHR2_CAG_UNITS
    + random_sequence(1800),
    "chr3": random_sequence(1200) + "GGGGCC" * 6 + random_sequence(1500),
}
CONTIG_ORDER = ["chr1", "chr2", "chr3"]
CHR2_CAG_START = 1200 + 36 + 2000
CHR2_CAG_END = CHR2_CAG_START + 3 * CHR2_CAG_UNITS

# Off-target regions. A read is assigned to one only when it lies entirely inside it, so the
# pure-repeat reads below are placed at [1100, 1350] on both contigs.
OFFTARGET_REGIONS = {"chr2": (1000, 1500), "chr3": (1000, 1500)}


def write_reference():
    with open(REF_FA, "w") as f:
        for name in CONTIG_ORDER:
            f.write(">%s\n" % name)
            sequence = contigs[name]
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i + 60] + "\n")
    pysam.faidx(REF_FA)


def write_catalog():
    catalog = [
        {
            "LocusId": "C9LIKE",
            "LocusStructure": "(%s)*" % C9LIKE_MOTIF,
            "ReferenceRegion": "chr1:%d-%d" % (C9LIKE_START, C9LIKE_END),
            "VariantType": "RareRepeat",
            "OfftargetRegions": ["%s:%d-%d" % (name, start, end) for name, (start, end) in OFFTARGET_REGIONS.items()],
        },
        {
            "LocusId": "CHR2_CAG",
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr2:%d-%d" % (CHR2_CAG_START, CHR2_CAG_END),
            "VariantType": "Repeat",
        },
    ]
    with open(CATALOG, "w") as f:
        json.dump(catalog, f, indent=4)
        f.write("\n")


# ----------------------------------------------------------------------------
# Reads
#
# Each read is a dict with: name, is_read1, contig, pos, seq, is_reverse, cigar (list of
# (operation, length) tuples, or None for an unmapped read), mapq, is_unmapped. Mate fields are
# filled in when the pair is written.
# ----------------------------------------------------------------------------
pairs = []


def mapped_read(contig, pos, seq, is_reverse, cigar=None, mapq=60):
    return dict(contig=contig, pos=pos, seq=seq, is_reverse=is_reverse,
                cigar=cigar if cigar is not None else [(0, len(seq))], mapq=mapq, is_unmapped=False)


def unmapped_read(seq):
    # Contig and position are filled in from the mate when the pair is written.
    return dict(contig=None, pos=None, seq=seq, is_reverse=False, cigar=None, mapq=0, is_unmapped=True)


def add_pair(name, read1, read2):
    if read1["is_unmapped"] and read2["is_unmapped"]:
        raise ValueError("both mates of %s are unmapped" % name)
    for read, mate in ((read1, read2), (read2, read1)):
        if read["is_unmapped"]:
            read["contig"], read["pos"] = mate["contig"], mate["pos"]
    pairs.append((name, read1, read2))


# --- C9LIKE, normal allele (reference length): spanning read pairs -------------------------------
def add_c9like_spanning_pairs(count):
    chr1 = contigs["chr1"]
    earliest_start = C9LIKE_END - READ_LENGTH + 5   # keep >= 5 bp of right flank
    latest_start = C9LIKE_START - 5                 # keep >= 5 bp of left flank
    window = latest_start - earliest_start + 1
    for i in range(count):
        read1_start = earliest_start + (i * 7) % window
        read2_start = read1_start + FRAGMENT_LENGTH - READ_LENGTH
        add_pair("c9like_normal_%03d" % i,
                 mapped_read("chr1", read1_start, chr1[read1_start:read1_start + READ_LENGTH], False),
                 mapped_read("chr1", read2_start, chr1[read2_start:read2_start + READ_LENGTH], True))


# --- C9LIKE, expanded allele ----------------------------------------------------------------------
EXPANDED_ALLELE = CHR1_LEFT_FLANK + C9LIKE_MOTIF * C9LIKE_EXPANDED_UNITS + CHR1_RIGHT_FLANK
EXPANDED_REPEAT_START = len(CHR1_LEFT_FLANK)
EXPANDED_REPEAT_END = EXPANDED_REPEAT_START + C9LIKE_EXPANDED_UNITS * len(C9LIKE_MOTIF)
EXPANSION_LENGTH = EXPANDED_REPEAT_END - EXPANDED_REPEAT_START - (C9LIKE_END - C9LIKE_START)


def place_expanded_allele_read(start, is_reverse):
    """Returns how BWA would place the expanded-allele read covering [start, start + READ_LENGTH).

    A read inside the repeat is returned as None (the caller decides where it goes). A read in a
    flank is mapped there. A read that crosses a repeat edge is mapped by its flank bases plus the
    reference repeat sequence next to them, and the rest is soft-clipped.
    """
    end = start + READ_LENGTH
    seq = EXPANDED_ALLELE[start:end]
    reference_repeat_length = C9LIKE_END - C9LIKE_START
    if start >= EXPANDED_REPEAT_START and end <= EXPANDED_REPEAT_END:
        return None
    if end <= EXPANDED_REPEAT_START:
        return mapped_read("chr1", start, seq, is_reverse)
    if start >= EXPANDED_REPEAT_END:
        return mapped_read("chr1", start - EXPANSION_LENGTH, seq, is_reverse)
    if start < EXPANDED_REPEAT_START:
        # Left flank, then repeat: aligned through the reference repeat sequence, then clipped.
        aligned_length = min(READ_LENGTH, EXPANDED_REPEAT_START - start + reference_repeat_length)
        cigar = [(0, aligned_length)] + ([(4, READ_LENGTH - aligned_length)] if aligned_length < READ_LENGTH else [])
        return mapped_read("chr1", start, seq, is_reverse, cigar)
    # Repeat, then right flank: clipped, then aligned from inside the reference repeat sequence.
    right_flank_bases = end - EXPANDED_REPEAT_END
    repeat_bases_aligned = min(reference_repeat_length, EXPANDED_REPEAT_END - start)
    aligned_length = right_flank_bases + repeat_bases_aligned
    cigar = [(4, READ_LENGTH - aligned_length), (0, aligned_length)]
    return mapped_read("chr1", C9LIKE_END - repeat_bases_aligned, seq, is_reverse, cigar)


def offtarget_pure_repeat_read(contig, pos, allele_start, is_reverse_on_allele):
    """A pure-repeat read from the expanded allele, mapped (MAPQ 0) inside an off-target region.

    chr2's tract has the allele's orientation; chr3's has the opposite one, so a read placed there
    maps to the other strand and its stored SEQ is the reverse complement of the allele's bases.
    """
    seq = EXPANDED_ALLELE[allele_start:allele_start + READ_LENGTH]
    if contig == "chr3":
        return mapped_read(contig, pos, reverse_complement(seq), not is_reverse_on_allele, mapq=0)
    return mapped_read(contig, pos, seq, is_reverse_on_allele, mapq=0)


def add_c9like_expanded_allele_pairs():
    """Tiles forward-reverse fragments across the expanded allele.

    - Both reads in flanks or crossing a repeat edge: mapped at chr1.
    - One read inside the repeat, its mate in a flank: the inside read is unmapped and placed at
      its mate (an anchored in-repeat read).
    - Both reads inside the repeat: placed in an off-target region. Most pairs land with both mates
      in the same off-target region (nearby); every fifth pair is split across chr2 and chr3 (far
      apart).
    """
    pure_repeat_pair_index = 0
    for fragment_number, fragment_start in enumerate(range(EXPANDED_REPEAT_START - 900, EXPANDED_REPEAT_END + 600, 12)):
        read2_start = fragment_start + FRAGMENT_LENGTH - READ_LENGTH
        read1 = place_expanded_allele_read(fragment_start, False)
        read2 = place_expanded_allele_read(read2_start, True)
        name = "c9like_expanded_%04d" % fragment_number

        if read1 is not None and read2 is not None:
            add_pair(name, read1, read2)
            continue
        if read1 is not None or read2 is not None:
            # Sequenced bases of an unmapped read: read2 comes from the reverse strand.
            if read1 is None:
                read1 = unmapped_read(EXPANDED_ALLELE[fragment_start:fragment_start + READ_LENGTH])
            else:
                read2 = unmapped_read(reverse_complement(EXPANDED_ALLELE[read2_start:read2_start + READ_LENGTH]))
            add_pair(name, read1, read2)
            continue

        offset = 1100 + (pure_repeat_pair_index * 13) % 130
        if pure_repeat_pair_index % 5 == 4:
            read1 = offtarget_pure_repeat_read("chr2", offset, fragment_start, False)
            read2 = offtarget_pure_repeat_read("chr3", offset + 20, read2_start, True)
        else:
            contig = "chr2" if pure_repeat_pair_index % 2 == 0 else "chr3"
            read1 = offtarget_pure_repeat_read(contig, offset, fragment_start, False)
            read2 = offtarget_pure_repeat_read(contig, offset + 120, read2_start, True)
            if contig == "chr3":
                # Opposite strand: the read that maps first on chr3 is read2's placement.
                read1["pos"], read2["pos"] = read2["pos"], read1["pos"]
        add_pair(name, read1, read2)
        pure_repeat_pair_index += 1


# --- CHR2_CAG: an ordinary locus with no off-target regions --------------------------------------
def add_chr2_cag_spanning_pairs(count):
    chr2 = contigs["chr2"]
    earliest_start = CHR2_CAG_END - READ_LENGTH + 5
    latest_start = CHR2_CAG_START - 5
    window = latest_start - earliest_start + 1
    for i in range(count):
        read1_start = earliest_start + (i * 7) % window
        read2_start = read1_start + FRAGMENT_LENGTH - READ_LENGTH
        add_pair("chr2cag_%03d" % i,
                 mapped_read("chr2", read1_start, chr2[read1_start:read1_start + READ_LENGTH], False),
                 mapped_read("chr2", read2_start, chr2[read2_start:read2_start + READ_LENGTH], True))


add_c9like_spanning_pairs(90)
add_c9like_expanded_allele_pairs()
add_chr2_cag_spanning_pairs(180)


# ----------------------------------------------------------------------------
# BAM
# ----------------------------------------------------------------------------
def aligned_length(read):
    return sum(length for operation, length in read["cigar"] if operation in (0, 2, 3, 7, 8))


def write_bam():
    contig_index = {name: i for i, name in enumerate(CONTIG_ORDER)}
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": len(contigs[name])} for name in CONTIG_ORDER],
    }

    records = []
    for name, read1, read2 in pairs:
        both_mapped_on_one_contig = (not read1["is_unmapped"] and not read2["is_unmapped"]
                                     and read1["contig"] == read2["contig"])
        if both_mapped_on_one_contig:
            left, right = (read1, read2) if read1["pos"] <= read2["pos"] else (read2, read1)
            template_length = right["pos"] + aligned_length(right) - left["pos"]
        for is_read1, read, mate in ((True, read1, read2), (False, read2, read1)):
            if both_mapped_on_one_contig:
                tlen = template_length if read is left else -template_length
            else:
                tlen = 0
            records.append(dict(name=name, is_read1=is_read1, read=read, mate=mate, tlen=tlen))

    records.sort(key=lambda r: (contig_index[r["read"]["contig"]], r["read"]["pos"],
                                r["read"]["is_unmapped"], r["name"], not r["is_read1"]))

    with pysam.AlignmentFile(BAM, "wb", header=header) as out:
        for record in records:
            read, mate = record["read"], record["mate"]
            segment = pysam.AlignedSegment(out.header)
            segment.query_name = record["name"]
            segment.query_sequence = read["seq"]
            segment.flag = 0
            segment.is_paired = True
            segment.is_read1 = record["is_read1"]
            segment.is_read2 = not record["is_read1"]
            segment.is_unmapped = read["is_unmapped"]
            segment.mate_is_unmapped = mate["is_unmapped"]
            segment.is_reverse = read["is_reverse"]
            segment.mate_is_reverse = mate["is_reverse"]
            segment.is_proper_pair = (not read["is_unmapped"] and not mate["is_unmapped"]
                                      and read["contig"] == mate["contig"])
            segment.reference_id = contig_index[read["contig"]]
            segment.reference_start = read["pos"]
            segment.mapping_quality = read["mapq"]
            if read["cigar"] is not None:
                segment.cigartuples = read["cigar"]
            segment.next_reference_id = contig_index[mate["contig"]]
            segment.next_reference_start = mate["pos"]
            segment.template_length = record["tlen"]
            segment.query_qualities = pysam.qualitystring_to_array("I" * len(read["seq"]))
            out.write(segment)

    pysam.index(BAM)


write_reference()
write_catalog()
write_bam()

print("Wrote:")
for path in (REF_FA, REF_FA + ".fai", CATALOG, BAM, BAM + ".bai"):
    print("  %s (%d bytes)" % (path, os.path.getsize(path)))
print("C9LIKE chr1:%d-%d, CHR2_CAG chr2:%d-%d, %d read pairs"
      % (C9LIKE_START, C9LIKE_END, CHR2_CAG_START, CHR2_CAG_END, len(pairs)))
