import polars as pl
from pathlib import Path

try:
    ann = snakemake.input.ann
    inferred = snakemake.input.inferred
    out = snakemake.output[0]
except NameError:
    ann = "annotations.jsonl"
    inferred = "inferred.jsonl"
    out = "brisi.jsonl"


def set_file(row):
    if row["country"] in "RS HR".split():
        return row["id"]
    else:
        return Path(row["audio"]).with_suffix("").name


ann = pl.read_ndjson(ann)
inf = pl.read_ndjson(inferred)
metadata = (
    pl.read_ndjson("../data/*/ParlaSpeech-*.jsonl")
    .with_columns(country=pl.col("id").str.extract(r".*ParlaMint-(.{2}).*"))
    .with_columns(
        file=pl.struct(["country", "id", "audio"]).map_elements(
            set_file, return_dtype=pl.String
        )
    )
    .filter(pl.col("file").is_in(ann["file"]))
    .unnest("speaker_info")
    .select("file Speaker_gender".split())
)


ann.join(inf, on="file", how="left", validate="m:1").sort("file").join(
    metadata, on="file", how="left", validate="m:1"
).write_ndjson(out)
2 + 2
