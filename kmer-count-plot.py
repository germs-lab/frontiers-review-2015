# Only slightly modified from Figure 7 in QP's notebook

import sys
import matplotlib.pyplot as plt
import pylab

file = open(sys.argv[1], 'r')

list1 = []
x1 = []

for line in file:
    line = line.rstrip()
    fields = line.split()
    if fields[1] != '0':
        x1.append(float(fields[0]))
        list1.append(float(fields[1]))


fig = plt.figure()

ax = fig.add_subplot(111)
fig.set_size_inches(10,7)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlabel("Starting position of k-mer in read",fontsize=12)
ax.set_ylabel("Number of abund=1 k-mers at that position",fontsize=12)


ax.plot(x1,list1,linestyle='-', label='k=17')
ax.legend(loc='upper left',prop={'size':12})
ax.set_xlim(0,100)

plt.savefig("perc_unique_pos.pdf")
