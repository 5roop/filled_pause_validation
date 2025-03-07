
# Automated evaluation report

Report compiled: 2025-03-07T11:07:12.443243+01:00

## Composition of available files:

Simple count of available files per language and per annotator:

| who    | lang   |   file count |
|:-------|:-------|-------------:|
| urban  | PL     |           11 |
| nela_1 | CZ     |          400 |
| peter  | PL     |            2 |
| drejc  | CZ     |           58 |

## Inter-annotator agreement

Displayed for all setups where it makes sense (at least two annotators.)

In case of three or more annotators, display pair-wise comparisons.

| lang   | who_1   | who_2   |   observed_agreement |   krippendorff_alpha |   common_files | precision   | recall      |
|:-------|:--------|:--------|---------------------:|---------------------:|---------------:|:------------|:------------|
| PL     | peter   | urban   |                1     |             1        |              2 | 1.0 <-> 1.0 | 1.0 <-> 1.0 |
| CZ     | nela_1  | drejc   |                0.866 |             0.696154 |             58 | 0.9 <-> 0.9 | 0.9 <-> 0.9 |

## Validation metrics:

| lang   | who    | how                    |   recall |   precision |   num_files |       F1 |
|:-------|:-------|:-----------------------|---------:|------------:|------------:|---------:|
| CZ     | drejc  | drop_short_and_initial |    0.914 |       0.86  |          58 | 0.886178 |
| CZ     | drejc  | raw                    |    0.927 |       0.835 |          58 | 0.878598 |
| CZ     | nela_1 | drop_short_and_initial |    0.958 |       0.874 |         400 | 0.914074 |
| CZ     | nela_1 | raw                    |    0.97  |       0.848 |         400 | 0.904906 |
| PL     | peter  | drop_short_and_initial |    1     |       1     |           2 | 1        |
| PL     | peter  | raw                    |    1     |       1     |           2 | 1        |
| PL     | urban  | drop_short_and_initial |    1     |       1     |          11 | 1        |
| PL     | urban  | raw                    |    1     |       1     |          11 | 1        |

## Inspection of the biggest discrepancies

Methodology: Count _number_ of positive events in predictions and annotations.
Order by biggest difference. The TextGrids and audio are available on our
server at `/cache/peterr/filled_pause_validation/annotation_eval/TG`

| lang   | file                           | who    |   abs_diff |
|:-------|:-------------------------------|:-------|-----------:|
| CZ     | 2020030415581612_281.38-288.52 | nela_1 |          3 |
| CZ     | 2022012613181332_256.23-271.19 | nela_1 |          3 |
| CZ     | 2022061510181032_206.52-222.18 | nela_1 |          3 |
| CZ     | 2022101815281542_470.52-480.92 | nela_1 |          3 |
| CZ     | 2014021119081922_520.16-535.18 | drejc  |          2 |
| CZ     | 2014032610281042_361.00-377.77 | drejc  |          2 |
| CZ     | 2020013112281242_124.92-134.81 | nela_1 |          2 |
| CZ     | 2020050518181832_677.79-684.86 | nela_1 |          2 |
| CZ     | 2022031113281342_335.85-348.25 | nela_1 |          2 |
| CZ     | 2022032918581912_236.43-249.13 | nela_1 |          2 |
| CZ     | 2022050420082022_565.36-579.73 | nela_1 |          2 |
| CZ     | 2022061415281542_188.70-195.46 | nela_1 |          2 |
| CZ     | 2022092715181532_217.04-231.94 | nela_1 |          2 |
| PL     | 3sjBHcXY-w8_15353.68-15370.68  | urban  |          2 |
| CZ     | 2014012215481602_466.41-483.03 | drejc  |          1 |
| CZ     | 2014020509080922_307.05-313.06 | drejc  |          1 |
| CZ     | 2014021119081922_520.16-535.18 | nela_1 |          1 |
| CZ     | 2014021210481102_250.80-263.65 | drejc  |          1 |
| CZ     | 2014021210481102_250.80-263.65 | nela_1 |          1 |
| CZ     | 2014021310481102_128.35-135.65 | drejc  |          1 |
| CZ     | 2014021310481102_128.35-135.65 | nela_1 |          1 |
| CZ     | 2014031814281442_319.60-338.99 | nela_1 |          1 |
| CZ     | 2014032010281042_437.15-450.97 | drejc  |          1 |
| CZ     | 2014032010281042_437.15-450.97 | nela_1 |          1 |
| CZ     | 2014032017381752_433.49-440.63 | nela_1 |          1 |
| CZ     | 2014042914181432_465.17-473.93 | drejc  |          1 |
| CZ     | 2014042914181432_465.17-473.93 | nela_1 |          1 |
| CZ     | 2014050710281042_140.71-152.91 | drejc  |          1 |
| CZ     | 2014050710281042_140.71-152.91 | nela_1 |          1 |
| CZ     | 2014050714081422_696.82-711.99 | drejc  |          1 |