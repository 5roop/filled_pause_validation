
# Automated evaluation report

Report compiled: 2025-03-07T09:26:49.793590+01:00

## Composition of available files:

|    | who    | lang   |   file count |
|---:|:-------|:-------|-------------:|
|  0 | drejc  | CZ     |           58 |
|  1 | nela_1 | CZ     |          400 |
|  2 | peter  | PL     |            2 |

## Inter-annotator agreement

|    | lang   | who_1   | who_2   |   observed_agreement |   krippendorff_alpha |   common_files | precision   | recall      |
|---:|:-------|:--------|:--------|---------------------:|---------------------:|---------------:|:------------|:------------|
|  0 | CZ     | drejc   | nela_1  |                0.866 |             0.696154 |             58 | 0.9 <-> 0.9 | 0.9 <-> 0.9 |

## Validation metrics:

|    | lang   | who    | how                    |   recall |   precision |   num_files |
|---:|:-------|:-------|:-----------------------|---------:|------------:|------------:|
|  0 | CZ     | drejc  | drop_short_and_initial |    0.914 |       0.86  |          58 |
|  1 | CZ     | drejc  | raw                    |    0.927 |       0.835 |          58 |
|  2 | CZ     | nela_1 | drop_short_and_initial |    0.958 |       0.874 |         400 |
|  3 | CZ     | nela_1 | raw                    |    0.97  |       0.848 |         400 |
|  4 | PL     | peter  | drop_short_and_initial |    1     |       1     |           2 |
|  5 | PL     | peter  | raw                    |    1     |       1     |           2 |
