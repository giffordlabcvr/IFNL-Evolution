#!/usr/bin/env python3

import argparse
import os
import xml.etree.ElementTree as ET


def find_xml_files(directory):
    return sorted(
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".xml")
    )


def get_qualifier(feature, name):
    for qualifier in feature.findall(".//GBQualifier"):
        qname = qualifier.findtext("GBQualifier_name")
        qvalue = qualifier.findtext("GBQualifier_value")
        if qname == name:
            return qvalue
    return None


def parse_record(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    if root.tag == "GBSeq":
        gbseq = root
    else:
        gbseq = root.find(".//GBSeq")

    if gbseq is None:
        raise ValueError("No GBSeq found")

    accession = (
        gbseq.findtext("GBSeq_primary-accession")
        or gbseq.findtext("GBSeq_accession-version")
    )

    sequence_length = int(gbseq.findtext("GBSeq_length"))

    cds_features = []

    for feature in gbseq.findall(".//GBFeature"):
        if feature.findtext("GBFeature_key") != "CDS":
            continue

        location = feature.findtext("GBFeature_location")

        cds_features.append({
            "location": location,
            "gene": get_qualifier(feature, "gene"),
            "protein_id": get_qualifier(feature, "protein_id"),
            "translation": get_qualifier(feature, "translation"),
            "codon_start": get_qualifier(feature, "codon_start"),
        })

    return {
        "accession": accession,
        "length": sequence_length,
        "cds_features": cds_features,
    }

def parse_simple_location(location):
    """
    Parse simple GenBank locations such as:
        79..1281
        <79..1281
        79..>1281

    Returns (start, end), or None for compound/unsupported locations.
    """

    if not location:
        return None

    cleaned = location.replace("<", "").replace(">", "")

    if ".." not in cleaned:
        return None

    if "join" in cleaned or "complement" in cleaned:
        return None

    start_text, end_text = cleaned.split("..", 1)

    try:
        return int(start_text), int(end_text)
    except ValueError:
        return None


def write_glue(records, source_name, output_file):

    with open(output_file, "w") as out:

        out.write("# Generated from GenBank XML CDS annotations\n\n")

        for record in records:

            accession = record["accession"]

            if len(record["cds_features"]) != 1:
                continue

            cds = record["cds_features"][0]
            coords = parse_simple_location(cds["location"])

            if coords is None:
                continue

            cds_start, cds_end = coords

            reference_name = "GB_REF_" + accession

            out.write(
                f"create reference {reference_name} "
                f"{source_name} {accession}\n\n"
            )

            out.write(f"reference {reference_name}\n\n")

            out.write("  add feature-location mRNA\n")
            out.write(
                f"  feature-location mRNA add segment 1 {record['length']}\n\n"
            )

            out.write("  add feature-location orf\n")
            out.write(
                f"  feature-location orf add segment "
                f"{cds_start} {cds_end}\n\n"
            )

            out.write("  exit\n\n")


def main():

    parser = argparse.ArgumentParser(
        description="Generate GLUE reference feature locations from GenBank XML"
    )

    parser.add_argument("source_directory")
    parser.add_argument("output_glue")
    parser.add_argument(
        "--source-name",
        required=True,
        help="GLUE source name containing the sequences"
    )

    args = parser.parse_args()

    xml_files = find_xml_files(args.source_directory)

    records = []

    for xml_file in xml_files:
        try:
            record = parse_record(xml_file)
            records.append(record)

            print(
                f"{record['accession']}: "
                f"{len(record['cds_features'])} CDS feature(s)"
            )

        except Exception as e:
            print(f"ERROR {xml_file}: {e}")

    write_glue(
        records,
        args.source_name,
        args.output_glue
    )

    print()
    print(f"Processed {len(records)} GenBank records")
    print(f"Wrote {args.output_glue}")


if __name__ == "__main__":
    main()



