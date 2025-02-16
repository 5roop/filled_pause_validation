import polars as pl
from pathlib import Path


dfcz = (
    pl.read_ndjson("data/ParlaSpeech-CZ/ParlaSpeech-CZ.v1.0.jsonl")
    .unnest("speaker_info")
    .filter(
        pl.col("audio_length").is_between(6, 20)
        & pl.col("Speaker_gender").is_in(["M", "F"])
    )
    .with_columns(
        pl.col("audio")
        .map_elements(lambda s: f"data/ParlaSpeech-CZ/{s}", return_dtype=pl.String)
        .alias("audio")
    )
    .filter(
        pl.col("audio").map_elements(
            lambda s: Path(s).exists(), return_dtype=pl.Boolean
        )
    )
)

print("CZ shape:", dfcz.shape)
dfpl = (
    pl.read_ndjson("data/ParlaSpeech-PL/ParlaSpeech-PL.v1.0.jsonl")
    .unnest("speaker_info")
    .filter(
        pl.col("audio_length").is_between(6, 20)
        & pl.col("Speaker_gender").is_in(["M", "F"])
    )
    .with_columns(
        pl.col("audio")
        .map_elements(lambda s: f"data/ParlaSpeech-PL/{s}", return_dtype=pl.String)
        .alias("audio")
    )
    .filter(
        pl.col("audio").map_elements(
            lambda s: Path(s).exists(), return_dtype=pl.Boolean
        )
    )
)
print("pl shape:", dfpl.shape)
for df, name in ((dfpl, "pl"), (dfcz, "cz")):
    r = []
    for gender in "MF":
        subset = df.filter(pl.col("Speaker_gender") == gender).sample(2500)
        r.append(subset)
    df = pl.concat(r, how="vertical_relaxed").select(
        ["audio", "audio_length", "Speaker_gender"]
    )
    df.write_ndjson(f"2_supersplits/{name}.jsonl")


2 + 2
