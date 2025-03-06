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
    from sklearn.metrics import recall_score, precision_score, accuracy_score

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


print(
    "File count:",
    df.group_by(["who", "lang"]).agg(pl.col("file").count().alias("file count")),
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

pl.Config.set_tbl_rows(300)
print(
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
)
