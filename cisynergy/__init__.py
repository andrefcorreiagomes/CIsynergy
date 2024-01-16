import dit
import numpy as np
import matplotlib.pyplot as plt
from .cisynergy import cisynergy

# This script reproduces the graph of Figure 3 and the synergy values of Table 9

pY00 = 0.25
pY01 = 0.25
pY10 = 0.25
pY11 = 0.25

xx = []
yy = []
yy_kolchinsky_synergy = [] 

for conditional_p in np.linspace(0,1,101):
    d = dit.Distribution(['000', '001', '011', '101', '110'], [pY00*conditional_p, pY00*(1-conditional_p), pY01, pY10, pY11])
    synergy = cisynergy(d)
    d.set_rv_names('XYT')
    xx.append(conditional_p)
    yy.append(synergy)
    yy_kolchinsky_synergy.append(dit.shannon.mutual_information(d, ['T'], ['X', 'Y']) - dit.shannon.mutual_information(d, ['T'], ['Y']))

plt.plot(xx,yy)
plt.plot(xx, yy_kolchinsky_synergy)
plt.legend(["$\mathregular{S^{CI}(Y_1, Y_2)}$", "$\mathregular{S^{d}(Y_1, Y_2)}$"], loc ="lower right")
plt.xlabel("r")
plt.ylabel("Synergy")
plt.show()

#####

distributions = []

distributions.append(dit.Distribution(['000', '011', '101', '110'], [0.25, 0.25, 0.25, 0.25])) #XOR

distributions.append(dit.Distribution(['000', '010', '100', '111'], [0.25, 0.25, 0.25, 0.25])) #AND

distributions.append(dit.Distribution(['000', '011', '102', '113'], [0.25, 0.25, 0.25, 0.25])) #COPY

distributions.append(dit.uniform(['000', '011', '101', '110', '222', '233', '323', '332'])) #RDNXOR

distributions.append(dit.uniform(['000', '011', '101', '110', '022', '033', '123', '132', '204', '215', '305', '314', '226', '237', '327', '336', '448', '459', '549', '558', '46a', '47b', '56b', '57a', '64c', '65d', '74d', '75c', '66e', '67f', '76f', '77e'])) #RDNUNQXOR

distributions.append(dit.uniform(['0000', '0101', '1011', '1110'])) #XORDUPLICATE

distributions.append(dit.uniform(['000', '100', '010', '111'])) #ANDDUPLICATE (Note that we removed the last source variable, Y_3, that satisfies Y_3=Y_1)

distributions.append(dit.uniform(['0000', '0111', '1011', '1100'])) #XORLOSES

distributions.append(dit.uniform(['0000', '1110', '2220', '3330', '2101', '3011', '0321', '1231'])) #XORMULTICOAL

for d in distributions:
    print('Synergy = ', cisynergy(d))