import polars as pl
from pydub import AudioSegment
from pathlib import Path
import datetime
from tqdm import tqdm

for what in "cz pl".split():
    df = pl.read_ndjson(f"4_validation_splits/{what}.jsonl", infer_schema_length=None)
    for row in tqdm(
        df.iter_rows(named=True), total=df.shape[0], desc=f"Segmenting {what}"
    ):
        filename = Path(row["audio"]).with_suffix("").name
        eaf_path = Path("5_final_empty_files", what, filename, f"{filename}.eaf")
        eaf_path.parent.mkdir(
            exist_ok=True,
            parents=True,
        )
        audio = AudioSegment.from_file(row["audio"])
        audio.export(
            eaf_path.with_suffix(".wav"),
            format="wav",
            parameters=["-ac", "1", "-ar", "16000"],
        )

        nsmap = {"xsi": "http://www.w3.org/2001/XMLSchema-instance"}

        from lxml import etree

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

        # Add the tier for "FilledPause"
        tier = etree.SubElement(
            root,
            "TIER",
            {"LINGUISTIC_TYPE_REF": "default-lt", "TIER_ID": "FilledPause"},
        )

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
