try:
    inpath = snakemake.input.master
    lang = snakemake.wildcards["lang"]
    file = snakemake.wildcards["file"]
    wav = snakemake.output.wav
except:
    inpath = "master.jsonl"
    lang = "PL"
    file = "r-CCvzhAo1M_4298.26-4304.28"
    wav = "TG/PL/r-CCvzhAo1M_4298.26-4304.28/r-CCvzhAo1M_4298.26-4304.28.wav"

from pathlib import Path
import polars as pl
from pydub import AudioSegment
import textgrid
from tqdm import tqdm


def find_audio(f: str):
    candidates = (
        list(Path("../data/").glob(f"**/{f}.flac"))
        + list(Path("../data/").glob(f"**/{f}.mp3"))
        + list(Path("../data/").glob(f"**/{f}.wav"))
    )
    if len(candidates) == 1:
        return str(candidates[0])
    else:
        print("Could not find file, continuing searching by id")
        candidates = (
            pl.read_ndjson("../data/*/ParlaSpeech-*.jsonl")
            .select(["audio", "id"])
            .filter(pl.col("id").str.contains(f))["audio"]
        )
        if len(candidates) == 1:
            hit = Path(candidates[0]).with_suffix("").name
            return find_audio(hit)
        else:
            raise FileNotFoundError(f"Expected 1 file, found {len(candidates)}")


df = pl.read_ndjson(inpath)

subset = df.filter(pl.col("lang").eq(lang) & pl.col("file").eq(file))
audio_path = Path(wav)
audio_path.parent.mkdir(parents=True, exist_ok=True)
# if not audio_path.exists():
audio = AudioSegment.from_file(find_audio(file))
audio.export(
    audio_path,
    format="wav",
    parameters=["-ac", "1", "-ar", "16000"],
)
audio_len = len(audio) / 1e3
# else:
# audio_len = len(AudioSegment.from_file(audio_path))
tg = textgrid.TextGrid()
y_pred_tier = textgrid.IntervalTier(
    name="y pred drop short and initial", minTime=0, maxTime=audio_len
)
y_pred = subset.head(1)["y_pred_drop_short_and_initial"][0]
for i in y_pred:
    y_pred_tier.addInterval(textgrid.Interval(round(i[0], 3), round(i[1], 3), "eee"))
tg.append(y_pred_tier)
y_pred_tier = textgrid.IntervalTier(
    name="y pred drop short", minTime=0, maxTime=audio_len
)
y_pred = subset.head(1)["y_pred_drop_short"][0]
for i in y_pred:
    y_pred_tier.addInterval(textgrid.Interval(round(i[0], 3), round(i[1], 3), "eee"))
tg.append(y_pred_tier)
for row in subset.iter_rows(named=True):
    tier = textgrid.IntervalTier(
        name=row["who"].replace("_", " "), minTime=0, maxTime=audio_len
    )
    for i in row["y_true"]:
        tier.addInterval(textgrid.Interval(round(i[0], 3), round(i[1], 3), "eee"))
    tg.append(tier)
tg.write(snakemake.output.tg)
