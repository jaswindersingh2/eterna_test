#!/bin/bash

########## for base-pair probability prediction. Give path to RNAfold directory if not installed globally. ####################
RNAfold -p ./1088_sequences/$1 | awk '{print $1}' > ./RNAfold_prob/$1.dbn


########## extract pair probability only from .ps file ####################
sed -ne '/%start of base pair probability data/,$ p' ./RNAfold_prob/$1_dp.ps | head -n -3 | tail -n +2 | awk '{print $1, $2, $3, $4}' | grep 'ubox' | awk '{print $1, $2, $3}' > ./RNAfold_prob/$1.prob
