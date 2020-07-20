import numpy as np
import pandas as pd
import pickle as pkl


#########  read file consists of name of each rna file ################
with open('ids') as f:
    ids = f.read().splitlines()

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


########### read either SPOT-RNA of RNAfold probabilites. uncomment either 1 #######################
	y_pred = np.loadtxt('SPOT-RNA_prob/' + id + '.prob', delimiter='\t')    # load SPOT-RNA bps probabilties  107 x 107
	#y_pred = RNAfold_bp_prob(id, seq)                                      # load RNAfold bps probabilties   107 x 107


	prob = np.sum(y_pred[0:79,0:79] + np.transpose(y_pred[0:79,0:79]), axis=0)  # convert upper triangular matrix to symmetric metric of size 79 x 79 and sum across one axis

	npair_prob = [1-i for i in prob]         # convert pair prob. to non-pair prob


########## ignore index of reactivity values less than 1e-5 (given in readme of eternabench) ###############
	ignore_index = []
	for i,I in enumerate(reactivity):
		if I < 0.00001:
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
temp2 = np.concatenate([i for i in all_pred_prob])
temp3 = np.concatenate([i for i in all_true_react])

###### single pcc value ###########
pcc_all = pcc = np.corrcoef(np.stack((np.array(temp2), np.array(temp3)), axis=0))[0][1]

print('mean pcc from individual RNA pcc = ', np.nanmean(save_pcc))
print('single pcc by concatenating all nts = ', pcc_all) 

