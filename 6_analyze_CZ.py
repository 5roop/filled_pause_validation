from pathlib import Path
import polars as pl
from tqdm import tqdm

try:
    df = pl.read_ndjson("6_cz_y_true_y_pred.jsonl")
except FileNotFoundError:
    annotated_files = list(Path("6_annotation_results/cz/").glob("*.ann.eaf"))

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

    def find_audio(f: str):
        candidates = list(
            Path("../parlaspeech_cz_fp_annotation/audio").glob(f"**/{f}.flac")
        )
        if len(candidates) == 1:
            return str(candidates[0])
        else:
            print("Could not find file")
            raise FileNotFoundError(f"Expected 1 file, found {len(candidates)}")

    annotations = [read_eaf(i) for i in annotated_files]

    df = (
        pl.DataFrame({"eaf_path": annotated_files, "y_true": annotations})
        .with_columns(
            pl.col("eaf_path")
            .map_elements(
                lambda p: p.with_suffix("").with_suffix("").name, return_dtype=pl.String
            )
            .alias("file"),
            pl.col("eaf_path").map_elements(str, return_dtype=pl.String),
        )
        .with_columns(
            pl.col("eaf_path")
            .map_elements(read_eaf, return_dtype=pl.List(pl.List(pl.Float64)))
            .alias("y_true")
        )
    )
    audios = pl.Series(
        "audio",
        [
            find_audio(i)
            for i in tqdm(df["file"].to_list(), desc="Finding audios, inshallah")
        ],
    )
    df = df.insert_column(1, audios)

    from transformers import (
        AutoFeatureExtractor,
        Wav2Vec2BertForAudioFrameClassification,
    )
    from datasets import Dataset, Audio
    import datasets
    import torch
    import numpy as np
    from pathlib import Path
    import polars as pl

    device = torch.device("cuda:2")
    model_name = "classla/wav2vecbert2-filledPause"
    feature_extractor = AutoFeatureExtractor.from_pretrained(model_name)
    model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(model_name).to(
        device
    )

    def frames_to_intervals(
        frames: list[int], drop_short=True, drop_initial=True, short_cutoff_s=0.08
    ) -> list[tuple[float]]:
        """Transforms a list of ones or zeros, corresponding to annotations on frame
        levels, to a list of intervals ([start second, end second]).

        Allows for additional filtering on duration (false positives are often short)
        and start times (false positives starting at 0.0 are often an artifact of
        poor segmentation).

        :param list[int] frames: Input frame labels
        :param bool drop_short: Drop everything shorter than short_cutoff_s, defaults to True
        :param bool drop_initial: Drop predictions starting at 0.0, defaults to True
        :param float short_cutoff_s: Duration in seconds of shortest allowable prediction, defaults to 0.08
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
                        round(ndf.loc[ei - 1, "time_s"], 3),
                    )
                )
        if drop_short and (len(results) > 0):
            results = [i for i in results if (i[1] - i[0] >= short_cutoff_s)]
        if drop_initial and (len(results) > 0):
            results = [i for i in results if i[0] != 0.0]
        return results

    def evaluator(chunks):
        sampling_rate = chunks["audio"][0]["sampling_rate"]
        with torch.no_grad():
            inputs = feature_extractor(
                [i["array"] for i in chunks["audio"]],
                return_tensors="pt",
                sampling_rate=sampling_rate,
            ).to(device)
            logits = model(**inputs).logits
        y_pred = np.array(logits.cpu()).argmax(axis=-1)
        intervals = [
            frames_to_intervals(i, drop_short=False, drop_initial=False)
            for i in y_pred.tolist()
        ]
        # print(intervals)
        return {"y_pred": intervals}

    ds = Dataset.from_pandas(df.to_pandas()).cast_column(
        "audio", Audio(16000, mono=True)
    )
    ds = ds.map(
        evaluator,
        batched=True,
        batch_size=20,
    )
    y_preds = ds["y_pred"]
    s = pl.Series("y_pred", y_preds)
    df = df.insert_column(1, s)
    df.write_ndjson("6_cz_y_true_y_pred.jsonl")

    2 + 2


from sklearn.metrics import classification_report, confusion_matrix


def is_overlapping(this, other):
    if (this[0] < other[1]) and (this[1] > other[0]):
        return True
    return False


import pandas as pd
import numpy as np

print("Performance on CZ:")
y_pred_events, y_true_events = [], []
FPs = []
FNs = []
for row in df.iter_rows(named=True):
    y_pred = row["y_pred"]
    y_true = row["y_true"]
    if (y_pred == []) and (y_true == []):
        y_pred_events.append(0)
        y_true_events.append(0)
    for l in y_pred + y_true:
        if any([is_overlapping(l, x) for x in y_pred + y_true if x != l]):
            y_pred_events.append(1)
            y_true_events.append(1)
        else:
            if l in y_pred:
                y_pred_events.append(1)
                y_true_events.append(0)
                FPs.append(l)
            else:
                y_true_events.append(1)
                y_pred_events.append(0)
                FNs.append(l)
pl.Config.set_tbl_cols(-1)
pl.Config.set_tbl_rows(-1)
pl.Config.set_tbl_width_chars(600)
pl.Config.set_fmt_str_lengths(50)
print(df.select(["file", "y_true", "y_pred"]))
print(
    "Classification report for Nela vs y_pred on event level: ",
    classification_report(y_true_events, y_pred_events, digits=3),
    sep="\n",
)
print("Confusion matrix for Nela vs y_pred:")
print(confusion_matrix(y_pred=y_pred_events, y_true=y_true_events))


Path("to_reannotate.txt").write_text(
    "\n".join(
        df.filter(pl.col("y_pred").ne([]) | pl.col("y_true").ne([]))["eaf_path"]
        .str.split("/")
        .list.get(2)
        .to_list()
    )
)
master_df = (
    pl.read_ndjson("data/ParlaSpeech-CZ/ParlaSpeech-CZ.v1.0.jsonl")
    .with_columns(
        file=pl.col("audio").str.split("/").list.last().str.replace(".flac", "")
    )
    .filter(pl.col("file").is_in(df["file"]))
    .unnest("speaker_info")
    .select(["file", "Speaker_gender"])
)

sanity_check = df.with_columns(pl.col("y_pred").eq([]).alias("is_empty")).join(
    master_df, on="file"
)

print("Sanity check: composition of sample:")
print(sanity_check.group_by(["Speaker_gender", "is_empty"]).len())


Path("7_inspection_tgs/cz").mkdir(exist_ok=True, parents=True)
from pydub import AudioSegment
from lxml import etree
import datetime
from textgrid import TextGrid, IntervalTier, Interval
from tqdm import tqdm

for row in tqdm(df.iter_rows(named=True), total=df.shape[0]):
    tg_path = Path("7_inspection_tgs/cz/", row["file"], row["file"] + ".TextGrid")
    tg_path.parent.mkdir(parents=True, exist_ok=True)

    audio = AudioSegment.from_file(row["audio"])
    audio.export(
        tg_path.with_suffix(".wav"),
        format="wav",
        parameters=["-ac", "1", "-ar", "16000"],
    )
    audio_len = 1e-3 * len(audio)

    tg = TextGrid()
    y_true_tier = IntervalTier(name="true", minTime=0, maxTime=audio_len)
    y_pred_tier = IntervalTier(name="pred", minTime=0, maxTime=audio_len)

    for event in row["y_true"]:
        y_true_tier.addInterval(Interval(event[0], event[1], "eee"))
    for event in row["y_pred"]:
        if event[0] == event[1]:
            event[1] = event[1] + 0.0001
        y_pred_tier.addInterval(Interval(event[0], event[1], "eee"))
    tg.append(y_true_tier)
    tg.append(y_pred_tier)
    tg.write(tg_path)
2 + 2
