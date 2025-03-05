from pathlib import Path
import polars as pl
df = pl.read_ndjson(snakemake.input[0])


def read_eaf(f: Path | str):
    from lxml import etree

    doc = etree.fromstring(Path(f).read_bytes())
    timestamps = {
        i.attrib.get("TIME_SLOT_ID"): 1e-3 * int(i.attrib.get("TIME_VALUE"))
        for i in doc.findall(".//TIME_SLOT")
    }
    annotations = [
        [
            timestamps.get(i.attrib.get("TIME_SLOT_REF1")),
            timestamps.get(i.attrib.get("TIME_SLOT_REF2")),
        ]
        for i in doc.findall(".//ALIGNABLE_ANNOTATION")
    ]

    return annotations

annotations = pl.Series("y_true", [read_eaf(i) for i in df["ann"]])
df = df.insert_column(1, annotations)

df.write_ndjson(snakemake.output[0])
