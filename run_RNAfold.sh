#!/bin/bash

########## for base-pair probability prediction ####################
RNAfold -p /home/jaswinder/Documents/eternagame-EternaBench-0010251e737b69f65fd9ef259025b36a90abfa2e/ChemMapping/data/preprocessing/1088_sequences/$1 | awk '{print $1}' > /home/jaswinder/Documents/eternagame-EternaBench-0010251e737b69f65fd9ef259025b36a90abfa2e/ChemMapping/data/preprocessing/RNAfold_prob/$1.dbn

sed -ne '/%start of base pair probability data/,$ p' /home/jaswinder/Documents/eternagame-EternaBench-0010251e737b69f65fd9ef259025b36a90abfa2e/ChemMapping/data/preprocessing/RNAfold_prob/$1_dp.ps | head -n -3 | tail -n +2 | awk '{print $1, $2, $3}' > /home/jaswinder/Documents/eternagame-EternaBench-0010251e737b69f65fd9ef259025b36a90abfa2e/ChemMapping/data/preprocessing/RNAfold_prob/$1.prob
