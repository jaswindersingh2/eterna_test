# eterna_test
test SPOT-RNA on eterna testset (R69) with 1088 RNA sequences.

Overview:
----
* `1088_reactivity` folder contains reactivity labels extracted form `/eternagame-EternaBench-0010251e737b69f65fd9ef259025b36a90abfa2e/ChemMapping/data/datasets_with_predictions/temp/Users/hwayment/das/github/EternaBench/ChemMapping/data/full_datasets/Round69_11May2020.json.zip` pandas dataframe saved under the index name `reactivity`.
* `1088_sequences` folder contains 1088 RNA sequences obtained from above dataframe under the index name `sequence`.
* `RNAfold_prob` folder contains base-pair probability obtained from RNAfold using `run_RNAfold.sh` script.
* `SPOT-RNA_prob` folder contains base-pair probability obtained from SPOT-RNA using standalone program available at [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA).  
* `corr.py` is the program to calculate the correlation between predicted probability and reactivity labels.

Results:
----
```
$ python3 corr.py  --predictor SPOT-RNA

SPOT-RNA
mean pcc from individual RNA pcc = 0.379
single pcc by concatenating all nts = 0.410

```

```
$ python3 corr.py  --predictor RNAfold

RNAfold
mean pcc from individual RNA pcc = 0.496
single pcc by concatenating all nts = 0.416

```
