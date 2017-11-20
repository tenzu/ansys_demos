import random as rdm
import numpy as np
NUM_F = 3   #number of fibres (number of CSYSs)
i,j = 0,0
tmp_list,tmp_CSYS,move_CSYS,rotation_CSYS = np.zeros((6)),np.zeros((NUM_F,6)),np.zeros((NUM_F,3)),np.zeros((NUM_F,3))
f1 = open('move_CSYS.txt','w')
f2 = open('rotation_CSYS.txt','w')
while i < NUM_F:
    while j < NUM_F*2:
        if j < NUM_F:
            tmp_list[j] = rdm.uniform(-100,0)
        else:
            tmp_list[j] = rdm.uniform(0,360)
        j+=1
    tmp_CSYS[i] = tmp_list
    tmp_list = np.zeros((6))
    j=0
    i+=1
i = 0
while i < NUM_F:
    move_CSYS[i] = tmp_CSYS[i][0:NUM_F]
    rotation_CSYS[i] = tmp_CSYS[i][NUM_F:NUM_F*2]
    for j in move_CSYS[i]:
        f1.write('%3.3f' % j + '\t')
    f1.write('\n')
    for k in rotation_CSYS[i]:
        f2.write('%3.3f' % k + '\t')
    f2.write('\n')
    i += 1
f1.close()
f2.close()
