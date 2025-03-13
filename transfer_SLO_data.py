import polars as pl
from lxml import etree
from pathlib import Path
import datetime


def frames_to_intervals(
    frames: list[int], drop_short=True, drop_initial=True, short_cutoff_s=0.08
) -> list[tuple[float]]:
    """Transforms a list of ones or zeros, corresponding to annotations on frame
    levels, to a list of intervals ([start second, end second]).

    Allows for additional filtering on duration (false positives are often
    short) and start times (false positives starting at 0.0 are often an
    artifact of poor segmentation).

    :param list[int] frames: Input frame labels
    :param bool drop_short: Drop everything shorter than short_cutoff_s,
        defaults to True
    :param bool drop_initial: Drop predictions starting at 0.0, defaults to True
    :param float short_cutoff_s: Duration in seconds of shortest allowable
        prediction, defaults to 0.08
    :return list[tuple[float]]: List of intervals [start_s, end_s]
    """
    from itertools import pairwise
    import pandas as pd

    results = []
    ndf = pd.DataFrame(
        data={
            "time_s": [0.020 * i for i in range(len(frames))],
            "frames": frames,
        }
    )
    ndf = ndf.dropna()
    indices_of_change = ndf.frames.diff()[ndf.frames.diff() != 0].index.values
    for si, ei in pairwise(indices_of_change):
        if ndf.loc[si : ei - 1, "frames"].mode()[0] == 0:
            pass
        else:
            results.append(
                (
                    round(ndf.loc[si, "time_s"], 3),
                    round(ndf.loc[ei, "time_s"], 3),
                )
            )
    if drop_short and (len(results) > 0):
        results = [i for i in results if (i[1] - i[0] >= short_cutoff_s)]
    if drop_initial and (len(results) > 0):
        results = [i for i in results if i[0] != 0.0]
    return results


def process_labels(s: str):
    l = eval(s)
    l = [{"0": 0, "filledPause": 1}[i] for i in l]
    return frames_to_intervals(l, drop_short=False, drop_initial=False)


df = pl.read_csv(
    "/cache/peterr/mezzanine_resources/ML/data/test_filledPause.csv"
).with_columns(
    pl.col("segment_path")
    .map_elements(lambda s: Path(s).name, return_dtype=pl.String)
    .alias("audio"),
    pl.col("labels").map_elements(
        lambda s: process_labels(s), return_dtype=pl.List(pl.List(pl.Float64))
    ),
)
for row in df.iter_rows(named=True):
    eaf_path = Path(
        "/cache",
        "peterr",
        "filled_pause_validation",
        "annotation_eval",
        "anns",
        "SL",
        "gold",
        row["audio"],
    ).with_suffix(".ann.eaf")
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

    # Add the tier for "human"
    human_tier = etree.SubElement(
        root,
        "TIER",
        {"LINGUISTIC_TYPE_REF": "default-lt", "TIER_ID": "human"},
    )

    for y in row["labels"]:
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
2 + 2
