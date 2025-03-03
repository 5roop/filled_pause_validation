# filled_pause_validation
Scripts for validating the Filled Pause Detector on RS, HR, PL, CZ ParlaSpeech data


* 0_download_data.sh: downloads, concatenates, decompresses RS, HR, PL, CZ data
* 1_inspect_previous_splits.py: compares the overall distributions with the distributions of the previous (HR and RS) splits. Insights:
the validation splits were 50-50 male-female, with durations ranging from 6 to 20 seconds. Other parameters seem to be randomized and not filtered in any way.
* 2_prepare_splits.py: Prepares a supersplit of 5000 sentences (2500 M, 2500 F) to be annotated with our model. Only instances that have audio available will be taken into account.
* 3_run_inference.py: Runs all the data through our model.
* 4_create_validation_datasets.py: From the available annotated sentences, it randomly selects 100 male sentences with automatically predicted Filled Pauses, 100 male sentences without Filled Pauses, and the same for female sentences.
* 5_create_ean_files.py: From the final selection, the audio files are converted to 16kHz mono wav, an empty ELAN file is created, and both are placed in their own directory.

## 2025-03-01T09:06:39

Due to a miscomunication, an older version of the CZ files were annotated. This will not be a problem, as the sample used is still valid, it just means that additional steps will have to be taken to analyze it.

## 2025-03-03T10:09:55

Results of the analysis:

```
Classification report for human vs y_pred on event level:
              precision    recall  f1-score   support

           0      0.948     0.682     0.793       292
           1      0.850     0.980     0.910       540

    accuracy                          0.875       832
   macro avg      0.899     0.831     0.852       832
weighted avg      0.885     0.875     0.869       832

Confusion matrix for human vs y_pred:
[[199  93]
 [ 11 529]]

```