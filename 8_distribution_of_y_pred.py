import polars as pl
import seaborn as sns


def is_overlapping(this, other):
    if (this[0] < other[1]) and (this[1] > other[0]):
        return True
    return False


df = (
    pl.read_ndjson("/cache/peterr/artur_fp_annotation/3_inferred_FP.jsonl")
    .filter(pl.col("filled_pauses").ne([]))
    .explode("filled_pauses")
    .with_columns(
        duration=(
            pl.col("filled_pauses").list.last() - pl.col("filled_pauses").list.first()
        )
    )
)


g = sns.displot(
    df.select(["duration"]),
    x="duration",
    multiple="stack",
    kind="hist",
    cumulative=True,
    stat="probability",
    binwidth=0.02,
    binrange=[0, 1],
    facet_kws=dict(sharey=False),
    height=4,
    aspect=1,
)

g.savefig("8_distribution_of_durations.png", dpi=300)
2 + 2
