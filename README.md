# eterna_test
test SPOT-RNA on eterna testset (R69) with 1088 RNA sequences.

Overview:
----
* `1088_reactivity` folder contains reactivity labels extracted form `/eternagame-EternaBench-0010251e737b69f65fd9ef259025b36a90abfa2e/ChemMapping/data/datasets_with_predictions/temp/Users/hwayment/das/github/EternaBench/ChemMapping/data/full_datasets/Round69_11May2020.json.zip` pandas dataframe saved under the index name `reactivity`.
* `1088_sequences` folder contains 1088 RNA sequences obtained from above dataframe under the index name `sequence`.
* `RNAfold_prob` folder contains base-pair probability obtained from RNAfold (ViennaRNA v2.4.14) using `run_RNAfold.sh` script for the full sequence length 107.
* `RNAfold_prob_1_79` folder contains base-pair probability obtained from RNAfold (ViennaRNA v2.4.14) using `run_RNAfold.sh` script for the truncated sequence length 1-79 nucleotides.
* `SPOT-RNA_prob` folder contains base-pair probability obtained from SPOT-RNA using standalone program available at [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA) for the full sequence length 107. 
* `SPOT-RNA_prob_1_79` folder contains base-pair probability obtained from SPOT-RNA using standalone program available at [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA) for the truncated sequence length 1-79 nucleotides.
* `corr.py` is the program to calculate the correlation between predicted probability and reactivity labels.

Results:
----
Need to install numpy, pandas, argparse, and scipy using either pip or conda.
 
Results from the truncated sequence base-pair probability for SPOT-RNA and RNAfold (ViennaRNA v2.4.14, default settings):
```
$ python3 corr.py  --predictor SPOT-RNA --truncation 1

SPOT-RNA
mean Pearson Correlation Coefficent from individual RNA pcc = 0.480
single Pearson Correlation Coefficent by concatenating all nts = 0.489
single Spearman Correlation Coefficent by concatenating all nts = 0.501

```

```
$ python3 corr.py  --predictor RNAfold --truncation 1

RNAfold
mean Pearson Correlation Coefficent from individual RNA pcc = 0.494
single Pearson Correlation Coefficent by concatenating all nts = 0.467
single Spearman Correlation Coefficent by concatenating all nts = 0.337

```

Results from the full sequence base-pair probability for SPOT-RNA and RNAfold (ViennaRNA v2.4.14, default settings):
```
$ python3 corr.py  --predictor SPOT-RNA --truncation 0

SPOT-RNA
mean Pearson Correlation Coefficent from individual RNA pcc = 0.379
single Pearson Correlation Coefficent by concatenating all nts = 0.427
single Spearman Correlation Coefficent by concatenating all nts = 0.419

```

```
$ python3 corr.py  --predictor RNAfold --truncation 0

RNAfold
mean Pearson Correlation Coefficent from individual RNA pcc = 0.498
single Pearson Correlation Coefficent by concatenating all nts = 0.485
single Spearman Correlation Coefficent by concatenating all nts = 0.340

```
Contact
====
jaswinder.singh3@griffithuni.edu.au, yaoqi.zhou@griffith.edu.au
