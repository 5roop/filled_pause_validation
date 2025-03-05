from pathlib import Path
import polars as pl
from tqdm import tqdm

try:
    available_anns = snakemake.input
except:
    available_anns = list(Path("anns").glob("**/*.ann.eaf"))
available_anns = list(map(str, available_anns))
df = pl.DataFrame({"ann": map(str, available_anns)}).with_columns(
    pl.col("ann")
    .str.extract(r"anns/(.*)/([^/]+)/(.*)\.ann\.eaf", group_index=1)
    .alias("lang"),
    pl.col("ann")
    .str.extract(r"anns/(.*)/([^/]+)/(.*)\.ann\.eaf", group_index=2)
    .alias("who"),
    pl.col("ann")
    .str.extract(r"anns/(.*)/([^/]+)/(.*)\.ann\.eaf", group_index=3)
    .alias("file"),
)


def find_audio(f: str):
    candidates = list(Path("../data/").glob(f"**/{f}.flac")) + list(
        Path("../data/").glob(f"**/{f}.mp3")
    )
    if len(candidates) == 1:
        return str(candidates[0])
    else:
        print("Could not find file")
        raise FileNotFoundError(f"Expected 1 file, found {len(candidates)}")


print("Finding audios...")
to_find = df["file"].unique()

from concurrent.futures import ProcessPoolExecutor

with ProcessPoolExecutor(max_workers=1000) as executor:
    audios = executor.map(find_audio, to_find)
print("Audios found")
audio_mapper = pl.DataFrame({"file": to_find, "audio_path": audios})

df.sort("file").write_ndjson(snakemake.output.files_df)
audio_mapper.sort("file").write_ndjson(snakemake.output.audios)
