#!/usr/bin/env python3
"""Generate a deterministic synthetic ExpansionHunter fixture for the IsInsidePseudoautosomalRegion catalog field.

Produces, under this script's directory:
  reference.fa (+ .fai)
  variant_catalog.json
  reads.bam (+ .bai)

Fourteen CAG loci, each heterozygous for a 12-unit (reference length) and a 16-unit allele:
  PAR_CAG_01 ... PAR_CAG_12  on chrX, marked "IsInsidePseudoautosomalRegion": true in the catalog
  NONPAR_CAG                 on chrX, not marked
  AUTOSOMAL_CAG              on chr1

See README.md.

Read-building conventions are the same as in ../offtarget_regions_fixtures/make_fixture.py: a reverse-strand
read's SEQ is stored in forward-reference orientation, and catalog ReferenceRegions are 0-based half-open.
Everything is seeded so the output bytes are identical on every run.
"""

import json
import os
import random

import pysam

RNG = random.Random(20261108)  # fixed seed; never use time-based randomness

HERE = os.path.dirname(os.path.abspath(__file__))
REF_FA = os.path.join(HERE, "reference.fa")
CATALOG = os.path.join(HERE, "variant_catalog.json")
BAM = os.path.join(HERE, "reads.bam")

READ_LENGTH = 150
FRAGMENT_LENGTH = 400
MIN_FLANK_IN_SPANNING_READ = 20
READ_PAIRS_PER_ALLELE = 50

MOTIF = "CAG"
REFERENCE_UNITS = 12
LONGER_ALLELE_UNITS = 16

# ----------------------------------------------------------------------------
# Reference and catalog
#
# chrX (58 kb): PAR_CAG_01 ... PAR_CAG_12 at 2000, 6000, ..., 46000, and NONPAR_CAG at 54000. The loci are 4 kb apart,
# so the 1 kb flanks ExpansionHunter extracts reads from around each repeat don't overlap.
# chr1 (6 kb): AUTOSOMAL_CAG at 3000.
# ----------------------------------------------------------------------------
PSEUDOAUTOSOMAL_LOCUS_IDS = ["PAR_CAG_%02d" % (i + 1) for i in range(12)]
LOCI = [
    # (locus id, contig, repeat start, is inside a pseudoautosomal region)
    *[(locus_id, "chrX", 2000 + 4000 * i, True) for i, locus_id in enumerate(PSEUDOAUTOSOMAL_LOCUS_IDS)],
    ("NONPAR_CAG", "chrX", 54000, False),
    ("AUTOSOMAL_CAG", "chr1", 3000, False),
]
CONTIG_LENGTHS = {"chr1": 6000, "chrX": 58000}
CONTIG_ORDER = ["chr1", "chrX"]
REPEAT_LENGTH = REFERENCE_UNITS * len(MOTIF)


def random_sequence(length):
    return "".join(RNG.choice("ACGT") for _ in range(length))


def build_contig(name):
    """Returns the contig's sequence: random bases with each of its loci's 12-unit repeat."""
    pieces = []
    position = 0
    for _, contig, start, _ in sorted((locus for locus in LOCI if locus[1] == name), key=lambda locus: locus[2]):
        pieces.append(random_sequence(start - position))
        pieces.append(MOTIF * REFERENCE_UNITS)
        position = start + REPEAT_LENGTH
    pieces.append(random_sequence(CONTIG_LENGTHS[name] - position))
    return "".join(pieces)


contigs = {name: build_contig(name) for name in CONTIG_ORDER}


def write_reference():
    with open(REF_FA, "w") as f:
        for name in CONTIG_ORDER:
            f.write(">%s\n" % name)
            sequence = contigs[name]
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i + 60] + "\n")
    pysam.faidx(REF_FA)


def write_catalog():
    catalog = []
    for locus_id, contig, start, is_inside_pseudoautosomal_region in LOCI:
        record = {
            "LocusId": locus_id,
            "LocusStructure": "(%s)*" % MOTIF,
            "ReferenceRegion": "%s:%d-%d" % (contig, start, start + REPEAT_LENGTH),
            "VariantType": "Repeat",
        }
        if is_inside_pseudoautosomal_region:
            record["IsInsidePseudoautosomalRegion"] = True
        catalog.append(record)
    with open(CATALOG, "w") as f:
        json.dump(catalog, f, indent=4)
        f.write("\n")


# ----------------------------------------------------------------------------
# Reads: for each locus, READ_PAIRS_PER_ALLELE read pairs from each allele. Read 1 spans the repeat with at
# least MIN_FLANK_IN_SPANNING_READ bases of flank on each side; read 2 lies in the right flank. A read from the
# 16-unit allele is aligned with a 12 bp insertion after the 6th repeat unit.
# ----------------------------------------------------------------------------
pairs = []


def add_pairs_for_allele(locus_id, contig, start, units):
    reference = contigs[contig]
    end = start + REPEAT_LENGTH
    inserted_length = (units - REFERENCE_UNITS) * len(MOTIF)
    haplotype = reference[:start] + MOTIF * units + reference[end:]
    haplotype_repeat_end = end + inserted_length
    earliest_start = haplotype_repeat_end - READ_LENGTH + MIN_FLANK_IN_SPANNING_READ
    latest_start = start - MIN_FLANK_IN_SPANNING_READ
    window = latest_start - earliest_start + 1
    insertion_position = start + 6 * len(MOTIF)  # on the reference, also on the haplotype
    for i in range(READ_PAIRS_PER_ALLELE):
        read1_start = earliest_start + (i * 7) % window
        read1_sequence = haplotype[read1_start:read1_start + READ_LENGTH]
        if inserted_length:
            bases_before_insertion = insertion_position - read1_start
            cigar = [(0, bases_before_insertion), (1, inserted_length),
                     (0, READ_LENGTH - bases_before_insertion - inserted_length)]
        else:
            cigar = [(0, READ_LENGTH)]
        # Read 2 starts in the right flank, so its haplotype and reference positions differ by the insertion.
        read2_haplotype_start = read1_start + FRAGMENT_LENGTH - READ_LENGTH
        read2_sequence = haplotype[read2_haplotype_start:read2_haplotype_start + READ_LENGTH]
        read2_reference_start = read2_haplotype_start - inserted_length
        assert reference[read2_reference_start:read2_reference_start + READ_LENGTH] == read2_sequence
        name = "%s_%dunits_%03d" % (locus_id.lower(), units, i)
        pairs.append((name,
                      dict(contig=contig, pos=read1_start, seq=read1_sequence, is_reverse=False, cigar=cigar),
                      dict(contig=contig, pos=read2_reference_start, seq=read2_sequence, is_reverse=True,
                           cigar=[(0, READ_LENGTH)])))


for locus_id, contig, start, _ in LOCI:
    for units in (REFERENCE_UNITS, LONGER_ALLELE_UNITS):
        add_pairs_for_allele(locus_id, contig, start, units)


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
        template_length = read2["pos"] + aligned_length(read2) - read1["pos"]
        records.append(dict(name=name, is_read1=True, read=read1, mate=read2, tlen=template_length))
        records.append(dict(name=name, is_read1=False, read=read2, mate=read1, tlen=-template_length))
    records.sort(key=lambda r: (contig_index[r["read"]["contig"]], r["read"]["pos"], r["name"], not r["is_read1"]))

    with pysam.AlignmentFile(BAM, "wb", header=header) as out:
        for record in records:
            read, mate = record["read"], record["mate"]
            segment = pysam.AlignedSegment(out.header)
            segment.query_name = record["name"]
            segment.query_sequence = read["seq"]
            segment.flag = 0
            segment.is_paired = True
            segment.is_proper_pair = True
            segment.is_read1 = record["is_read1"]
            segment.is_read2 = not record["is_read1"]
            segment.is_reverse = read["is_reverse"]
            segment.mate_is_reverse = mate["is_reverse"]
            segment.reference_id = contig_index[read["contig"]]
            segment.reference_start = read["pos"]
            segment.mapping_quality = 60
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
print("%d read pairs" % len(pairs))
