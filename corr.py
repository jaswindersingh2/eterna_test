import numpy as np
import pandas as pd
import pickle as pkl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--predictor', default='SPOT-RNA', type=str, help='provide either SPOT-RNA or RNAfold', metavar='')
args = parser.parse_args()


#########  read file consists of name of each rna file ################
with open('ids') as f:
    ids = f.read().splitlines()

###### read all the reactivities and concatenate together to find out reactivity values above 95% percentile ##############
all_reactivites = []
for id in ids[0:]:
	with open('1088_reactivity/' + str(id)) as f:
		temp = pd.read_csv(f, comment='#', header=None, skiprows=[0]).values

	reactivity = [float(i) for i in temp[1][0].split(' ') if i != ''] # make list of reactivity values in float
	all_reactivites.append(reactivity)

	assert len(reactivity) == 79     # check no. of reactivity values. should be 79

all_reactivities_concat = np.concatenate([i for i in all_reactivites])   # concatenate all the values and make 1D array
thres_remove_above_95 = np.percentile(all_reactivities_concat, 95)       # evaluate theshold for reactivities values above and below the 95% cut-off

######## --------------------- parse base-pair probability RNAfold output ---------------------------- #########################
def RNAfold_bp_prob(id, seq):
    with open('RNAfold_prob/' + str(id) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None).values
    #print(temp.shape)
    output_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:,0]:
        a = i.split(' ')
        #print(a)
        output_pred[int(a[0]) - 1, int(a[1]) - 1] = float(a[2])
    #print(output_pred)
    return output_pred

all_pred_prob = []
all_true_react = []
save_pcc = []

for id in ids[0:]:

###### read sequence and reactivity of query sequence ##############
	with open('1088_reactivity/' + str(id)) as f:
		temp = pd.read_csv(f, comment='#', header=None, skiprows=[0]).values
	seq = [i for i in temp[0][0]]
	assert len(seq) == 107             # check its length. should be 107 

	reactivity = [float(i) for i in temp[1][0].split(' ') if i != ''] # make list of reactivity values in float

	assert len(reactivity) == 79     # check no. of reactivity values. should be 79


########### read either SPOT-RNA or RNAfold probabilites #######################
	if args.predictor == 'RNAfold':
		y_pred = RNAfold_bp_prob(id, seq)                                      # load RNAfold bps probabilties   107 x 107  2D array
	else:
		y_pred = np.loadtxt('SPOT-RNA_prob/' + id + '.prob', delimiter='\t')    # load SPOT-RNA bps probabilties  107 x 107  2D array


	prob = np.sum(y_pred[0:79,0:79] + np.transpose(y_pred[0:79,0:79]), axis=0)  # convert upper triangular matrix to symmetric metric of size 79 x 79 and sum across one axis

	npair_prob = [1-i for i in prob]         # convert pair probability to non-pair probability


########## ignore index of reactivity values less than 1e-5 and above 95% ###############
	ignore_index = []
	for i,I in enumerate(reactivity):
		if I < 0.00001 or I > thres_remove_above_95:
			ignore_index.append(i)	
	npair_prob = [I for i,I in enumerate(npair_prob) if i not in ignore_index]
	reactivity = [I for i,I in enumerate(reactivity) if i not in ignore_index]

######## append non-pair prob and reactivity to calculate single pcc value  ##########
	all_pred_prob.append(npair_prob)
	all_true_react.append(reactivity)

#########  calculate pcc of individual rna  ################
	pcc = np.corrcoef(np.stack((np.array(npair_prob), np.array(reactivity)), axis=0))[0][1]
	save_pcc.append(pcc)
	

####### concatenate all 1-dimensional un-paired prob. and reactivtites ###########
all_probabilities = np.concatenate([i for i in all_pred_prob])
all_reactivities = np.concatenate([i for i in all_true_react])
#print(len(all_probabilities), len(all_reactivities))

###### single pcc value ###########
pcc_all = np.corrcoef(np.stack((np.array(all_probabilities), np.array(all_reactivities)), axis=0))[0][1]

print('\n'+args.predictor)
print('mean pcc from individual RNA pcc = {:.3f}'.format(np.nanmean(save_pcc)))
print('single pcc by concatenating all nts = {:.3f}'.format(pcc_all)) 
print()

