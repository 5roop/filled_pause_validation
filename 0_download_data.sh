mkdir -p data/ParlaSpeech-PL/
cd data/ParlaSpeech-PL
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1686{/ParlaSpeech-PL.v1.0.jsonl.gz,/ParlaSpeech-PL.v1.0.part1.tgz,/ParlaSpeech-PL.v1.0.part2.tgz,/README.txt}
cat ParlaSpeech-PL*part* > ParlaSpeech-PL-v1.0.tgz
tar -xzvf ParlaSpeech-PL-v1.0.tgz
gzip -dk ParlaSpeech-PL.v1.0.jsonl.gz
# rm *.tgz *.gz
cd ../..


mkdir -p data/ParlaSpeech-RS
cd data/ParlaSpeech-RS
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1834{/ParlaSpeech-RS.v1.0.jsonl.gz,/ParlaSpeech-RS.v1.0.part1.tgz,/ParlaSpeech-RS.v1.0.part2.tgz,/README.txt}
cat ParlaSpeech-RS-v1.0.part*.tgz > ParlaSpeech-RS-v1.0.tgz
tar -xzvf ParlaSpeech-RS-v1.0.tgz
gzip -dk ParlaSpeech-RS.v1.0.jsonl.gz
# rm *.tgz
cd ../..


mkdir -p data/ParlaSpeech-CZ/
cd data/ParlaSpeech-CZ
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1785{/ParlaSpeech-CZ.v1.0.jsonl.gz,/ParlaSpeech-CZ.v1.0.part1.tgz,/ParlaSpeech-CZ.v1.0.part2.tgz,/ParlaSpeech-CZ.v1.0.part3.tgz,/ParlaSpeech-CZ.v1.0.part4.tgz}
cat ParlaSpeech-CZ*part* > ParlaSpeech-CZ-v1.0.tgz
tar -xzvf ParlaSpeech-CZ-v1.0.tgz
gzip -dk ParlaSpeech-CZ.v1.0.jsonl.gz
# rm *.tgz
cd ../..



mkdir -p data/ParlaSpeech-HR
cd data/ParlaSpeech-HR
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1914{/ParlaSpeech-HR.v2.0.jsonl.gz,/ParlaSpeech-HR.v2.0.part1.tgz,/ParlaSpeech-HR.v2.0.part2.tgz,/ParlaSpeech-HR.v2.0.part3.tgz,/ParlaSpeech-HR.v2.0.part4.tgz,/ParlaSpeech-HR.v2.0.part5.tgz,/ParlaSpeech-HR.v2.0.part6.tgz,/README.txt}
cat ParlaSpeech-HR*part* > ParlaSpeech-HR-v2.0.tgz
tar -xzvf ParlaSpeech-HR-v2.0.tgz
gzip -dk ParlaSpeech-HR.v2.0.jsonl.gz
# rm *.tgz
cd ../..
