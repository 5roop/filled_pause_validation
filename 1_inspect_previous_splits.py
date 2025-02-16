from pathlib import Path
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

old_splits = Path("/cache", "peterr", "ps_fp_validation", "y_pred_y_true.jsonl")
rs_metadata_path = Path("data", "ParlaSpeech-RS", "ParlaSpeech-RS.v1.0.jsonl")
hr_metadata_path = Path("data", "ParlaSpeech-HR", "ParlaSpeech-HR.v2.0.jsonl")
assert all([old_splits.exists(), rs_metadata_path.exists(), hr_metadata_path.exists()])
df = pl.read_ndjson(old_splits)
ids = df["name"].to_list()

rs = (
    pl.read_ndjson(rs_metadata_path)
    .with_columns(
        pl.col("id").is_in(ids).alias("in_sample"), pl.lit("HR").alias("country")
    )
    .unnest("speaker_info")
).with_row_index()
hr = (
    pl.read_ndjson(hr_metadata_path)
    .with_columns(
        pl.col("id").is_in(ids).alias("in_sample"), pl.lit("RS").alias("country")
    )
    .unnest("speaker_info")
).with_row_index()

df = pl.concat([rs, hr], how="vertical_relaxed").filter(pl.col("audio_length") < 50)
print(
    df.filter(pl.col("in_sample").eq(True))
    .group_by("Speaker_gender country".split())
    .agg(pl.col("index").count())
)
dfp = df.to_pandas()
print("Plotting")
for feature in "audio_length index Speaker_gender Speaker_role Speaker_birth Speaker_MP Speaker_party Party_status".split():
    sns.displot(
        dfp,
        x=feature,
        hue="in_sample",
        col="country",
        stat="density",
        common_norm=False,
        multiple="dodge",
    )
    plt.gcf().savefig(f"1_images/{feature}.png")
2 + 2
