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
    from sklearn.metrics import recall_score, precision_score

    y_pred_events, y_true_events = [], []
    FPs = []
    FNs = []
    TPs = []
    for y_pred, y_true in zip(y_preds.to_list(), y_trues.to_list()):
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
        # "y_pred": y_pred_events,
        # "y_true": y_true_events,
        "recall": recall_score(y_true_events, y_pred_events),
        "precision": precision_score(y_true_events, y_pred_events),
    }


how_cols = df.select(pl.selectors.starts_with("y_pred")).columns
other_cols = df.select(~pl.selectors.starts_with("y_pred")).columns

molten = df.unpivot(
    on=how_cols, index=other_cols, variable_name="how", value_name="y_pred"
).with_columns(pl.col("how").str.replace(r"y_pred_", ""))

print(
    molten.group_by(["lang", "who", "how"])
    .agg(
        pl.map_groups(
            exprs=["y_true", "y_pred"],
            function=lambda serieses: intervals_to_events(serieses[0], serieses[1]),
        ).alias("stats")
    )
    .unnest("stats")
)

2 + 2
