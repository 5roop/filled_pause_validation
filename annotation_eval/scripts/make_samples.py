import polars as pl
from pathlib import Path
from shutil import copy
from lxml import etree
import datetime

try:
    inpath = snakemake.input[0]
    output = snakemake.output[0]
except NameError:
    inpath = "samples.jsonl"
    output = "brisi"

df = pl.read_ndjson(inpath)
for row in df.iter_rows(named=True):
    wav = Path("TG", row["lang"], row["file"], f"""{row["file"]}.wav""")
    destination_wav = Path(
        "EAF_samples", row["lang"], row["file"], row["file"]
    ).with_suffix(".wav")
    destination_wav.parent.mkdir(exist_ok=True, parents=True)
    copy(wav, destination_wav)
    eaf_path = destination_wav.with_suffix(".eaf")
    nsmap = {"xsi": "http://www.w3.org/2001/XMLSchema-instance"}
    root = etree.Element(
        "ANNOTATION_DOCUMENT",
        {
            "AUTHOR": "Peter Rupnik, JSI",
            "DATE": datetime.datetime.now().isoformat(),
            "FORMAT": "3.0",
            "VERSION": "3.0",
            "{http://www.w3.org/2001/XMLSchema-instance}noNamespaceSchemaLocation": "http://www.mpi.nl/tools/elan/EAFv3.0.xsd",
        },
        nsmap=nsmap,
    )
    counter = 0
    # Add header elements
    header = etree.SubElement(
        root,
        "HEADER",
        {
            "MEDIA_FILE": eaf_path.with_suffix(".wav").name.encode("utf8"),
            "TIME_UNITS": "milliseconds",
        },
    )

    # Add time order (empty for now)
    time_order = etree.SubElement(root, "TIME_ORDER")

    # Add the tier for "machine"
    machine_tier = etree.SubElement(
        root,
        "TIER",
        {"LINGUISTIC_TYPE_REF": "default-lt", "TIER_ID": "machine"},
    )
    # Add the tier for "human"
    human_tier = etree.SubElement(
        root,
        "TIER",
        {"LINGUISTIC_TYPE_REF": "default-lt", "TIER_ID": "human"},
    )

    for y in row["y_true"]:
        start, stop = str(int(1000 * y[0])), str(int(1000 * y[1]))
        ts1, ts2 = f"ts{counter}", f"ts{counter + 1}"
        counter += 2
        annotation_value = etree.Element("ANNOTATION_VALUE")
        annotation_value.text = "eee"
        alignable_annotation = etree.Element("ALIGNABLE_ANNOTATION")
        alignable_annotation.attrib["ANNOTATION_ID"] = f"a{counter}"
        alignable_annotation.attrib["TIME_SLOT_REF1"] = ts1
        alignable_annotation.attrib["TIME_SLOT_REF2"] = ts2
        alignable_annotation.append(annotation_value)

        annotation = etree.Element("ANNOTATION")
        annotation.append(alignable_annotation)
        human_tier.append(annotation)

        time_slot = etree.Element("TIME_SLOT")
        time_slot.attrib["TIME_SLOT_ID"] = ts1
        time_slot.attrib["TIME_VALUE"] = start
        time_order.append(time_slot)

        time_slot = etree.Element("TIME_SLOT")
        time_slot.attrib["TIME_SLOT_ID"] = ts2
        time_slot.attrib["TIME_VALUE"] = stop
        time_order.append(time_slot)
    for y in row["y_pred"]:
        start, stop = str(int(1000 * y[0])), str(int(1000 * y[1]))
        ts1, ts2 = f"ts{counter}", f"ts{counter + 1}"
        counter += 2
        annotation_value = etree.Element("ANNOTATION_VALUE")
        annotation_value.text = "eee"
        alignable_annotation = etree.Element("ALIGNABLE_ANNOTATION")
        alignable_annotation.attrib["ANNOTATION_ID"] = f"a{counter}"
        alignable_annotation.attrib["TIME_SLOT_REF1"] = ts1
        alignable_annotation.attrib["TIME_SLOT_REF2"] = ts2
        alignable_annotation.append(annotation_value)

        annotation = etree.Element("ANNOTATION")
        annotation.append(alignable_annotation)
        machine_tier.append(annotation)

        time_slot = etree.Element("TIME_SLOT")
        time_slot.attrib["TIME_SLOT_ID"] = ts1
        time_slot.attrib["TIME_VALUE"] = start
        time_order.append(time_slot)

        time_slot = etree.Element("TIME_SLOT")
        time_slot.attrib["TIME_SLOT_ID"] = ts2
        time_slot.attrib["TIME_VALUE"] = stop
        time_order.append(time_slot)

    # Add linguistic type
    linguistic_type = etree.SubElement(
        root,
        "LINGUISTIC_TYPE",
        {
            "GRAPHIC_REFERENCES": "false",
            "LINGUISTIC_TYPE_ID": "default-lt",
            "TIME_ALIGNABLE": "true",
        },
    )

    # Create the XML tree
    tree = etree.ElementTree(root)
    etree.indent(tree)
    2 + 2
    # Write the XML to a file
    with open(str(eaf_path), "wb") as file:
        tree.write(file, encoding="utf-8", xml_declaration=True, pretty_print=True)
Path(output).write_text("\n")
