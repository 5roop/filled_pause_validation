from pathlib import Path

import datasets
import numpy as np
import polars as pl
import torch
from datasets import Audio, Dataset
from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification

file = snakemake.wildcards["file"]


def find_audio(f: str):
    candidates = list(Path("../data/").glob(f"**/{f}.flac")) + list(
        Path("../data/").glob(f"**/{f}.mp3")
    )
    if len(candidates) == 1:
        return str(candidates[0])
    else:
        print("Could not find file")
        raise FileNotFoundError(f"Expected 1 file, found {len(candidates)}")


audio = find_audio(file)

df = pl.DataFrame({"file": [file], "audio_path": [audio]})

device = torch.device("cpu")
model_name = "classla/wav2vecbert2-filledPause"
feature_extractor = AutoFeatureExtractor.from_pretrained(model_name)
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(model_name).to(device)


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
                    round(ndf.loc[ei - 1, "time_s"], 3),
                )
            )
    if drop_short and (len(results) > 0):
        results = [i for i in results if (i[1] - i[0] >= short_cutoff_s)]
    if drop_initial and (len(results) > 0):
        results = [i for i in results if i[0] != 0.0]
    return results


def evaluator(chunks):
    sampling_rate = chunks["audio_path"][0]["sampling_rate"]
    with torch.no_grad():
        inputs = feature_extractor(
            [i["array"] for i in chunks["audio_path"]],
            return_tensors="pt",
            sampling_rate=sampling_rate,
        ).to(device)
        logits = model(**inputs).logits
    y_pred = np.array(logits.cpu()).argmax(axis=-1)
    return {
        "y_pred_raw": [
            frames_to_intervals(i, drop_short=False, drop_initial=False)
            for i in y_pred.tolist()
        ],
        "y_pred_drop_short": [
            frames_to_intervals(i, drop_short=True, drop_initial=False)
            for i in y_pred.tolist()
        ],
        "y_pred_drop_initial": [
            frames_to_intervals(i, drop_short=False, drop_initial=True)
            for i in y_pred.tolist()
        ],
        "y_pred_drop_short_and_initial": [
            frames_to_intervals(i, drop_short=True, drop_initial=True)
            for i in y_pred.tolist()
        ],
    }


ds = Dataset.from_pandas(df.to_pandas()).cast_column(
    "audio_path", Audio(16000, mono=True)
)
ds = ds.map(
    evaluator,
    batched=True,
    batch_size=1,
)
df = ds.to_polars().select(["file", pl.selectors.starts_with("y_pred")])
try:
    df.write_ndjson(snakemake.output[0])
except NameError:
    df.write_ndjson("brisi.jsonl")
