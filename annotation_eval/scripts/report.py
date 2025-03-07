try:
    infile = snakemake.input[0]
except NameError:
    infile = "master.jsonl"

import polars as pl

df = pl.read_ndjson(infile)

langs = df["lang"].unique()


def is_overlapping(this, other):
    if (this[0] < other[1]) and (this[1] > other[0]):
        return True
    return False


def intervals_to_events(y_trues, y_preds):
    from sklearn.metrics import accuracy_score, precision_score, recall_score

    y_pred_events, y_true_events = [], []
    FPs = []
    FNs = []
    TPs = []
    try:
        y_preds, y_trues = y_preds.to_list(), y_trues.to_list()
    except:
        pass
    # Move events a smidge:
    y_preds_shifted = []
    for instance in y_preds:
        instance_shifted = []
        for event in instance:
            instance_shifted.append([event[0] + 1e-6, event[1] + 1e-6])
        y_preds_shifted.append(instance_shifted)
    y_preds = y_preds_shifted
    for y_pred, y_true in zip(y_preds, y_trues):
        if (y_pred == []) and (y_true == []):
            y_pred_events.append(0)
            y_true_events.append(0)
        for event in y_pred + y_true:
            if any([is_overlapping(event, x) for x in y_pred + y_true if x != event]):
                y_pred_events.append(1)
                y_true_events.append(1)
                TPs.append(event)
            else:
                if event in y_pred:
                    y_pred_events.append(1)
                    y_true_events.append(0)
                    FPs.append(event)
                else:
                    y_true_events.append(1)
                    y_pred_events.append(0)
                    FNs.append(event)
    return {
        "y_pred": y_pred_events,
        "y_true": y_true_events,
        "recall": round(recall_score(y_true_events, y_pred_events), 3),
        "precision": round(precision_score(y_true_events, y_pred_events), 3),
        "accuracy": round(accuracy_score(y_true_events, y_pred_events), 3),
    }


file_count = df.group_by(["who", "lang"]).agg(
    pl.col("file").count().alias("file count")
)
print(
    "File count:",
    file_count,
)


# Krippendorff:
def generate_pairs(who_list):
    from itertools import combinations

    return [list(i) for i in combinations(who_list, 2)]


def calculate_krippendorf_alpha(first, second):
    import krippendorff
    import numpy as np

    return krippendorff.alpha(np.array([first, second]))


def extract_results_from_k(row, k):
    who1 = row["who_1"]  # E.g. test
    who2 = row["who_2"]  # E.g. nela_1
    lang = row["lang"]

    subset = k.filter(
        pl.col(f"lang_{who1}").eq(lang) & pl.col(f"lang_{who2}").eq(lang)
    ).filter(
        pl.col(f"y_true_{who1}").is_not_null() & pl.col(f"y_true_{who2}").is_not_null()
    )
    r = intervals_to_events(
        subset[f"y_true_{who1}"],
        subset[f"y_true_{who2}"],
    )
    r_inv = intervals_to_events(
        subset[f"y_true_{who2}"],
        subset[f"y_true_{who1}"],
    )
    return {
        "observed_agreement": r["accuracy"],
        "krippendorff_alpha": calculate_krippendorf_alpha(r["y_pred"], r["y_true"]),
        "common_files": subset.shape[0],
        "precision": f"""{r["precision"]} <-> {r_inv["precision"]}""",
        "recall": f"""{r["recall"]} <-> {r_inv["recall"]}""",
    }


k = df.select(["file", "who", "lang", "y_true"]).pivot("who", index="file")
iaa = (
    (
        df.group_by("lang")
        .agg(pl.col("who").unique().alias("who_unique"))
        .with_columns(
            pl.col("who_unique")
            .map_elements(generate_pairs, return_dtype=pl.List(pl.List(pl.String)))
            .alias("who_pairs")
        )
    )
    .explode("who_pairs")
    .select(
        pl.col("lang"),
        pl.col("who_pairs").list.get(0).alias("who_1"),
        pl.col("who_pairs").list.get(1).alias("who_2"),
    )
    .drop_nulls()
    .with_columns(
        pl.struct(["who_1", "who_2", "lang"])
        .map_elements(
            lambda row: extract_results_from_k(row, k), return_dtype=pl.Struct
        )
        .alias("stats")
    )
    .unnest("stats")
)
print(iaa)


# Stats
how_cols = df.select(pl.selectors.starts_with("y_pred")).columns
other_cols = df.select(~pl.selectors.starts_with("y_pred")).columns

molten = df.unpivot(
    on=how_cols, index=other_cols, variable_name="how", value_name="y_pred"
).with_columns(pl.col("how").str.replace(r"y_pred_", ""))

metrics = (
    molten.group_by(["lang", "who", "how"])
    .agg(
        pl.map_groups(
            exprs=["y_true", "y_pred"],
            function=lambda serieses: intervals_to_events(serieses[0], serieses[1]),
        ).alias("stats"),
        pl.col("ann").n_unique().alias("num_files"),
    )
    .unnest("stats")
    .select("lang who how recall precision num_files".split())
    .sort("lang who how".split())
    .filter(pl.col("how").is_in(["raw", "drop_short_and_initial"]))
    .with_columns((2 / (1 / pl.col("recall") + 1 / pl.col("precision"))).alias("F1"))
)
pl.Config.set_tbl_rows(100)
print(metrics)

biggest_differences = (
    df.with_columns(
        (
            pl.col("y_true").list.len().cast(pl.Int64)
            - pl.col("y_pred_drop_short_and_initial").list.len().cast(pl.Int64)
        )
        .abs()
        .alias("abs_diff")
    )
    .sort("abs_diff file".split(), descending=True)
    .select("lang file abs_diff".split())
    .unique(subset=["file"], maintain_order=True, keep="first")
    .filter(pl.col("abs_diff").ge(1))
)

print(biggest_differences)

from datetime import datetime

import pytz
from jinja2 import Template
from pathlib import Path

when = datetime.now(tz=pytz.timezone("Europe/Ljubljana")).isoformat()
template = """
# Automated evaluation report

Report compiled: {{ when }}

## Composition of available files:

Simple count of available files per language and per annotator:

{{ file_count }}

## Inter-annotator agreement

Displayed for all setups where it makes sense (at least two annotators.)

In case of three or more annotators, display pair-wise comparisons.

{{ iaa }}

## Validation metrics:

{{ metrics }}

## Inspection of the biggest discrepancies

Methodology: Count _number_ of positive events in predictions and annotations.
Order by biggest difference. The TextGrids and audio are available on our
server at `/cache/peterr/filled_pause_validation/annotation_eval/TG`

{{ biggest_differences }}


## Manual analysis:

FP: false positive

FN: false negative

All comments based on auscultative examination

| lang   | file                           |   comment  |
|:-------|:-------------------------------|:-----------|
| CZ     | 2022101815281542_470.52-480.92 |          One FP by model |
| CZ     | 2022061510181032_206.52-222.18 |          One FP by model, 2 FN by annotator|
| CZ     | 2022012613181332_256.23-271.19 |          3 FN by annotator |
| CZ     | 2020030415581612_281.38-288.52 |          3 FN by annotator|
| PL     | 3sjBHcXY-w8_15353.68-15370.68  |          3:1 phenomenon |
| CZ     | 2022092715181532_217.04-231.94 |          twice 2:1 phenomenon |
| CZ     | 2022061415281542_188.70-195.46 |          2 FN by annotator |
| CZ     | 2022050420082022_565.36-579.73 |          One prominent /aam/ perhaps lexical? One probably FN by annotator|
| CZ     | 2022032918581912_236.43-249.13 |          One probable FP by model, one probable FN by annotator |
| CZ     | 2022031113281342_335.85-348.25 |          Probable FN by annotator, one short FN by annotator|
| CZ     | 2020050518181832_677.79-684.86 |          2 consecutive FN by model at the beginning . Raw has them as one filled pause, but we cut it afterwards in post processing|
| CZ     | 2020013112281242_124.92-134.81 |          2 FN by model|
| CZ     | 2014032610281042_361.00-377.77 |          1 FN by model at 0.0 that we cut in post processing. 1 FN by model, detected by 1 annotator only |
| CZ     | 2014021119081922_520.16-535.18 |          1 debatable FP (/bih ə/) by model, one FN by one annotator|
"""

2 + 2
Path(snakemake.output[0]).write_text(
    Template(template).render(
        dict(
            when=when,
            file_count=file_count.to_pandas().to_markdown(index=False),
            iaa=iaa.to_pandas().to_markdown(index=False),
            metrics=metrics.to_pandas().to_markdown(index=False),
            biggest_differences=biggest_differences.to_pandas().to_markdown(
                index=False
            ),
        )
    )
)
