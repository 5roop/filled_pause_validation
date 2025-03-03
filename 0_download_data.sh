mkdir -p data/ParlaSpeech-HR
cd data/ParlaSpeech-HR
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1914/ParlaSpeech-HR.v2.0.jsonl.gz
gzip -dk ParlaSpeech-HR.v2.0.jsonl.gz
rm *.tgz



cd ../..
mkdir -p data/ParlaSpeech-RS
cd data/ParlaSpeech-RS
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1834/ParlaSpeech-RS.v1.0.jsonl.gz
gzip -dk ParlaSpeech-RS.v1.0.jsonl.gz
rm *.tgz
cd ../..


mkdir -p data/ParlaSpeech-PL/
cd data/ParlaSpeech-PL
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1686{/ParlaSpeech-PL.v1.0.jsonl.gz,/ParlaSpeech-PL.v1.0.part1.tgz,/ParlaSpeech-PL.v1.0.part2.tgz,/README.txt}
ls *part* | xargs -I % sh -c 'tar -kxf %'
gzip -dk ParlaSpeech-PL.v1.0.jsonl.gz
rm *.tgz *.gz
cd ../..




mkdir -p data/ParlaSpeech-CZ/
cd data/ParlaSpeech-CZ
pwd
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1785{/ParlaSpeech-CZ.v1.0.jsonl.gz,/ParlaSpeech-CZ.v1.0.part1.tgz,/ParlaSpeech-CZ.v1.0.part2.tgz,/ParlaSpeech-CZ.v1.0.part3.tgz,/ParlaSpeech-CZ.v1.0.part4.tgz}
ls *part* | xargs -I % sh -c 'tar -kxf %'
gzip -dk ParlaSpeech-CZ.v1.0.jsonl.gz
rm *.tgz
cd ../..



