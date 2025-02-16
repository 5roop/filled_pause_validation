# filled_pause_validation
Scripts for validating the Filled Pause Detector on RS, HR, PL, CZ ParlaSpeech data


* 0_download_data.sh: downloads, concatenates, decompresses RS, HR, PL, CZ data
* 1_inspect_previous_splits.py: compares the overall distributions with the distributions of the previous (HR and RS) splits. Insights:
the validation splits were 50-50 male-female, with durations ranging from 6 to 20 seconds. Other parameters seem to be randomized and not filtered in any way.
* 2_prepare_splits.py: Prepares a supersplit of 5000 sentences (2500 M, 2500 F) to be annotated with our model. Only instances that have audio available will be taken into account.
* 3_run_inference.py: Runs all the data through our model.
* 4_create_validation_datasets.py: From the available annotated sentences, it randomly selects 100 male sentences with automatically predicted Filled Pauses, 100 male sentences without Filled Pauses, and the same for female sentences.
* 5_create_ean_files.py: From the final selection, the audio files are converted to 16kHz mono wav, an empty ELAN file is created, and both are placed in their own directory.

