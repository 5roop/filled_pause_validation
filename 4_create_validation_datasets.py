import polars as pl

what = "cz"
for what in "pl cz".split():
    r = []
    df = pl.read_ndjson(f"3_inference/{what}.jsonl").with_columns(
        pl.col("y_pred").eq([]).alias("is_empty")
    )
    for gender in "MF":
        for is_empty in [True, False]:
            subset = df.filter(
                (pl.col("Speaker_gender") == gender) & (pl.col("is_empty") == is_empty)
            ).sample(100)
            r.append(subset)
    df = pl.concat(r, how="vertical_relaxed")
    print(
        df.group_by(pl.col("is_empty Speaker_gender ".split())).agg(
            pl.col("y_pred").count()
        )
    )
    df.write_ndjson(f"4_validation_splits/{what}.jsonl")

2 + 2
