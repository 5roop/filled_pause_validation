import polars as pl

try:
    ann = snakemake.input.ann
    inferred = snakemake.input.inferred
    out = snakemake.output[0]
except NameError:
    ann = "annotations.jsonl"
    inferred = "inferred.jsonl"
    out = "brisi.jsonl"


ann = pl.read_ndjson(ann)
inf = pl.read_ndjson(inferred)

ann.join(inf, on="file", how="left", validate="m:1").sort("file").write_ndjson(out)
2 + 2
