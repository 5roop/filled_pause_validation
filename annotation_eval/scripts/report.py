try:
    infile = snakemake.input[0]
    metrics_plot = snakemake.output.metricsplot
    gender_plot = snakemake.output.genderplot
except NameError:
    infile = "master.jsonl"
    metrics_plot = "images/metrics.png"
    gender_plot = "images/gender.png"

import polars as pl
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

df = pl.read_ndjson(infile)

langs = df["lang"].unique()


def is_overlapping(this, other):
    if (this[0] < other[1]) and (this[1] > other[0]):
        return True
    return False


def intervals_to_events(y_trues, y_preds):
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

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
    inhibited_events = []
    for y_pred, y_true in zip(y_preds, y_trues):
        if (y_pred == []) and (y_true == []):
            y_pred_events.append(0)
            y_true_events.append(0)
            continue
        for event in y_pred:
            overlapping = [
                is_overlapping(event, x) for x in y_true if x not in inhibited_events
            ]
            if any(overlapping):
                y_pred_events.append(1)
                y_true_events.append(1)
                inhibited_events.extend([x for x in y_true if is_overlapping(event, x)])
            else:
                y_pred_events.append(1)
                y_true_events.append(0)
        for event in y_true:
            if event not in inhibited_events:
                y_pred_events.append(0)
                y_true_events.append(1)
    return {
        "y_true": y_true_events,
        "y_pred": y_pred_events,
        "recall": round(
            recall_score(y_true_events, y_pred_events, zero_division=1.0), 3
        ),
        "precision": round(
            precision_score(y_true_events, y_pred_events, zero_division=1.0), 3
        ),
        "accuracy": round(accuracy_score(y_true_events, y_pred_events), 3),
        "f1": round(f1_score(y_true_events, y_pred_events, zero_division=1.0), 3),
    }


y_true = [[[1, 2], [2.2, 2.6], [3, 4]]]
y_pred = [[[0, 5]]]
intervals_to_events(y_true, y_pred)
file_count = df.group_by(
    [
        "lang",
        "who",
    ]
).agg(pl.col("file").count().alias("file count"))
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
).with_columns(
    pl.col("how").str.replace(r"y_pred_", ""),
    # pl.struct(["y_true", "y_pred"])
    # .map_elements(
    #     lambda row: intervals_to_events([row["y_true"]], [row["y_pred"]]),
    #     return_dtype=pl.Struct,
    # )
    # .alias("rowstats"),
)

metrics = (
    molten.group_by(
        [
            "lang",
            "who",
            "how",
        ]
    )
    .agg(
        pl.map_groups(
            exprs=["y_true", "y_pred"],
            function=lambda serieses: intervals_to_events(serieses[0], serieses[1]),
        ).alias("stats"),
        pl.col("ann").n_unique().alias("num_files"),
    )
    .unnest("stats")
    .with_columns((2 / (1 / pl.col("recall") + 1 / pl.col("precision"))).alias("F1"))
    .select("lang who how recall precision F1 num_files".split())
    .filter(
        pl.col("how").is_in(
            ["raw", "drop_short", "drop_short_and_initial", "drop_initial"]
        )
    )
    .sort("lang who F1 how".split(), descending=[False, False, True, False])
)
pl.Config.set_tbl_rows(100)
print(metrics)

g = sns.relplot(
    metrics.unpivot(
        ["recall", "precision", "F1"],
        index="lang who how".split(),
        variable_name="metric",
        value_name="value",
    ).sort(
        pl.col("how").map_elements(
            lambda h: {
                "raw": 0,
                "drop_initial": 1,
                "drop_short": 2,
                "drop_short_and_initial": 3,
            }[h],
            return_dtype=pl.Int64,
        )
    ),
    x="how",
    row="metric",
    style="who",
    hue="lang",
    kind="line",
    y="value",
    facet_kws={"sharey": False, "sharex": False},
)
sns.move_legend(g, "center left", bbox_to_anchor=(1, 0.5))
plt.tight_layout()
g.savefig(metrics_plot)
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
    .filter(pl.col("abs_diff").ge(2))
)

print(biggest_differences)


gender = (
    molten.group_by(["lang", "who", "how", "Speaker_gender"])
    .agg(
        pl.map_groups(
            exprs=["y_true", "y_pred"],
            function=lambda serieses: intervals_to_events(serieses[0], serieses[1]),
        ).alias("stats"),
        pl.col("ann").n_unique().alias("num_files"),
    )
    .unnest("stats")
    .with_columns((2 / (1 / pl.col("recall") + 1 / pl.col("precision"))).alias("F1"))
    .select("lang who how recall precision F1 Speaker_gender num_files".split())
    .filter(
        pl.col("how").is_in(
            [
                "raw",
                "drop_short",
                "drop_short_and_initial",
                "drop_initial",
            ]
        )
    )
    .sort(
        pl.col("how").map_elements(
            lambda h: {
                "raw": 0,
                "drop_initial": 1,
                "drop_short": 2,
                "drop_short_and_initial": 3,
            }[h],
            return_dtype=pl.Int64,
        )
    )
)

g = sns.relplot(
    gender,
    x="Speaker_gender",
    y="F1",
    hue="lang",
    style="who",
    row="how",
    kind="line",
    facet_kws={"sharey": False, "sharex": False},
)
sns.move_legend(g, "center left", bbox_to_anchor=(1, 0.5))
plt.tight_layout()
g.savefig(gender_plot)

from scipy.stats import ranksums

overall_gender_stats = df.group_by("Speaker_gender").agg(
    pl.map_groups(
        exprs=["y_true", "y_pred_drop_short_and_initial"],
        function=lambda serieses: intervals_to_events(serieses[0], serieses[1])["f1"],
    ).alias("F1_score"),
    pl.map_groups(
        exprs=["y_true", "y_pred_drop_short_and_initial"],
        function=lambda serieses: [
            intervals_to_events(pl.Series("bluu", [t]), pl.Series("bluu", [p]))["f1"]
            for t, p in zip(serieses[0], serieses[1])
        ],
    ).alias("individual_f1s"),
)

per_language_gender_stats = (
    df.group_by("lang Speaker_gender".split())
    .agg(
        pl.map_groups(
            exprs=["y_true", "y_pred_drop_short_and_initial"],
            function=lambda serieses: intervals_to_events(serieses[0], serieses[1])[
                "f1"
            ],
        ).alias("F1_score")
    )
    .sort(["lang", "Speaker_gender"])
)

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

## Gender-based analysis:
The tables show results on `drop_short_and_initial`. For other post-processing methods, check out the plot below.

Overall:

{{ overall_gender_stats }}

Per country:

{{  per_language_gender_stats}}

### A graphical representation of gender disparities for individual languages, annotators, and  postprocessing steps:

![](images/gender.png)
## Validation metrics:

{{ metrics }}

### How do post-processing approaches affect the performance?

![plot of postprocessing steps](images/metrics.png)



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
| CZ     | 2022101815281542_470.52-480.92 |          One FP by model, 2 likely FN by annotator |
| CZ     | 2022061510181032_206.52-222.18 |          One FP by model, 2 FN by annotator|
| CZ     | 2022012613181332_256.23-271.19 |          3 FN by annotator |
| CZ     | 2020030415581612_281.38-288.52 |          3 FN by annotator|
| PL     | 3sjBHcXY-w8_15353.68-15370.68  |          3:1 phenomenon |
| CZ     | 2022092715181532_217.04-231.94 |          twice 2:1 phenomenon |
| CZ     | 2022061415281542_188.70-195.46 |          2 FN by annotator |
| CZ     | 2022050420082022_565.36-579.73 |          One FP: prominent /aam/ perhaps lexical? One probably FN by annotator|
| CZ     | 2022032918581912_236.43-249.13 |          One probable FP by model, one probable FN by annotator |
| CZ     | 2022031113281342_335.85-348.25 |          Probable FN by annotator, one short FN by annotator|
| CZ     | 2020050518181832_677.79-684.86 |          2 consecutive FN by model at the beginning . Raw has them as one filled pause, but we cut it afterwards in post processing|
| CZ     | 2020013112281242_124.92-134.81 |          2 FN by model|
| CZ     | 2014032610281042_361.00-377.77 |          1 FN by model at 0.0 that we cut in post processing. 1 FN by model, detected by 1 annotator only |
| CZ     | 2014021119081922_520.16-535.18 |          1 debatable FP (/bih ə/) by model, one FN by one annotator|
"""

2 + 2
Path(snakemake.output.report).write_text(
    Template(template).render(
        dict(
            when=when,
            file_count=file_count.sort("file count", descending=True)
            .to_pandas()
            .to_markdown(index=False),
            iaa=iaa.sort("observed_agreement", descending=True)
            .to_pandas()
            .to_markdown(index=False),
            metrics=metrics.to_pandas().to_markdown(index=False),
            biggest_differences=biggest_differences.to_pandas().to_markdown(
                index=False
            ),
            per_language_gender_stats=per_language_gender_stats.to_pandas().to_markdown(
                index=False
            ),
            overall_gender_stats=overall_gender_stats.select(
                "Speaker_gender F1_score".split()
            )
            .to_pandas()
            .to_markdown(index=False),
        )
    )
)
