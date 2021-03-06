# eterna_test
test SPOT-RNA on eterna testset (R69) with 1088 RNA sequences.

Overview:
----
* `1088_reactivity` folder contains reactivity labels extracted from `/eternagame-EternaBench-0010251e737b69f65fd9ef259025b36a90abfa2e/ChemMapping/data/datasets_with_predictions/temp/Users/hwayment/das/github/EternaBench/ChemMapping/data/full_datasets/Round69_11May2020.json.zip` pandas dataframe saved under the index name `reactivity`.
* `1088_sequences` folder contains 1088 RNA sequences obtained from above dataframe under the index name `sequence`.
* `RNAfold_prob` folder contains base-pair probability obtained from RNAfold ([ViennaRNA v2.4.14](https://www.tbi.univie.ac.at/RNA/#download)) using `run_RNAfold.sh` script for the full sequence length 107.
* `SPOT-RNA_prob` folder contains base-pair probability obtained from SPOT-RNA using standalone program available at [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA) for the full sequence length 107. 
* `corr.py` is the program to calculate the correlation between predicted probability and reactivity labels.

Results:
----
Need to install numpy, pandas, argparse, and scipy using either pip or conda.
 
```
$ python3 corr.py  --predictor SPOT-RNA

SPOT-RNA
mean Pearson Correlation Coefficent = 0.548
mean Spearman Correlation Coefficent = 0.518

```

```
$ python3 corr.py  --predictor RNAfold

RNAfold
mean Pearson Correlation Coefficent = 0.672
mean Spearman Correlation Coefficent = 0.667

```
Contact
====
jaswinder.singh3@griffithuni.edu.au, yaoqi.zhou@griffith.edu.au
