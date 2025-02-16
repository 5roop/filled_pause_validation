from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification
from datasets import Dataset, Audio
import datasets
import torch
import numpy as np
from pathlib import Path
import polars as pl


device = torch.device("cuda:6")
model_name = "classla/wav2vecbert2-filledPause"
feature_extractor = AutoFeatureExtractor.from_pretrained(model_name)
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(model_name).to(device)


def frames_to_intervals(frames: list[int]) -> list[tuple[float]]:
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
                (round(ndf.loc[si, "time_s"], 3), round(ndf.loc[ei - 1, "time_s"], 3))
            )
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
    intervals = [frames_to_intervals(i) for i in y_pred.tolist()]
    # print(intervals)
    return {"y_pred": intervals}


for what in "cz pl".split():
    f = f"2_supersplits/{what}.jsonl"
    df = pl.read_ndjson(f)
    ds = Dataset.from_pandas(df.to_pandas()).cast_column(
        "audio", Audio(16000, mono=True)
    )
    ds = ds.map(
        evaluator,
        batched=True,
        batch_size=30,
    )
    y_preds = ds["y_pred"]
    s = pl.Series("y_pred", y_preds)
    df = df.insert_column(1, s)
    df.write_ndjson(f"3_inference/{what}.jsonl")

    2 + 2
