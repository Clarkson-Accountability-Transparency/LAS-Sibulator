# Notes on sibulator

## Requirements
1. pandas
2. numpy
3. beautifulsoup (Only really used with "Strider" formats.) also needs lxml

Here is the list of what I was running with 
`pip freeze`
beautifulsoup4==4.11.1
bs4==0.0.1
lxml==4.9.0
numpy==1.22.4
pandas==1.4.2
python-dateutil==2.8.2
pytz==2022.1
PyYAML==6.0
six==1.16.0
soupsieve==2.3.2.post1

if this is in a reqs.txt, then `pip install -r reqs.txt` should make this environment.


## Making input files from LRMix files
1. The Alele naming is different for vWA. in sibulator it wants vWA vs VWA in the LRMix files
2. The file is organized as was the original in that funky two same named column format that we converted for LRMix. I wrote a small python script to convert the csv file. 


## Running.
All arguments are required. should be run as follows
python sibulator.py --known_sample_path <our sampleFile> --allele_freq_data [nist|strider] --subpop

### --known_sample_path
*REQUIRED*
This will be our converted csv file. not really a "path" but a CSV file. (or .txt that will be run with "eval")

### --allele_freq_data [nist|strider]
*REQURIED*
Bad code alert it must be 'nist' or 'strider' if not it will just fail with a KeyError as it doesn't check the else case.
strider doesn't seem to work without some additional help "bs4.FeatureNotFound: Couldn't find a tree builder with the features you requested: xml. Do you need to install a parser library?" (This error means you have to pip install lxml.)


### --subpop
*REQUIRED*
From the NIST selection this can be AfAm, Asian, Cauc, or Hisp.
Failure to specify one of those specific values results in a KeyError with the bad value you specified.

What's more is that I don't think this works with the STRidER data because those keys aren't in the xml file.
IF you use strider. then this must be one of these:
['AUSTRIA', 'BELGIUM', 'BOSNIA AND HERZEGOWINA', 'CZECH REPUBLIC', 'DENMARK', 'FINLAND', 'FRANCE', 'GERMANY', 'GREECE', 'HUNGARY', 'IRELAND', 'MONTENEGRO', 'NORWAY', 'POLAND', 'SLOVAKIA', 'SLOVENIA', 'SPAIN', 'SWEDEN', 'SWITZERLAND', 'THAILAND', 'Asia', 'Europe', 'Entire Database']


### --num_siblings
*REQUIRED*
The number of siblings to generate


### --dest_path
*REQUIRED*
This is actually should be csvFileOut as that is what it is for the sibling files.




## Programatic warnings
1. Do not feed this files that are *.txt it runs eval(f.read()) on those files. 
2. see --allele_freq_data section, it MUST be either nist or strider or ir will throw a KeyError since no dictionary gets loaded.
3. strider seems to need more help to get running.
4. There 

## Thoughts
I think it could be worthwhile rewriting this to make it work for us instead of using it as is. It would be worthwile making this an object vs the procedural code that it is.

If you want to run this directly.
import sibulator

  known_sample = sibulator.load_known_sample(known_sample_path)
    sibling_samples = sibulator.generate_sibling_samples(known_sample, allele_freq_data, subpop, num_siblings)

