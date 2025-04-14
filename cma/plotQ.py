import numpy as np
import sys
import matplotlib

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (15.0, 10.0)  #specify size of images produced by matplotlib
plt.rcParams.update({'font.size': 25})         #font size of graphs labels

Qfile    = sys.argv[1]
outfname = sys.argv[2]

print("Data will be loaded from the following file: %s" % Qfile)

Nt = 0
NR = 0

Rgrid = []
tgrid = []
Qdata = []

with open(Qfile,"r") as Qstream:
    for line in Qstream:
        if "#" not in line: #do not read comment lines
            if len(line.strip())!=0:  #skip empty lines
                data = line.split()
                time = float(data[0])
                R    = float(data[1])
                Q    = float(data[2])

                Qdata.append(Q)
                if Nt==0:
                    Rgrid.append(R)
                    NR += 1
            else:
                tgrid.append(time)
                Nt += 1

Qdata = np.array(Qdata)
Qdata = np.transpose(Qdata.reshape(Nt,NR))

#array of levels along z axis of plot
Qmin=np.amin(Qdata)
Qmax=np.amax(Qdata)
levels = np.linspace(Qmin,Qmax,101)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Time [fs]')
ax.set_ylabel('Molecular axis')
ax.set_title("Hole density",y=1.01,fontsize=25)

ax.contourf(tgrid,Rgrid,Qdata)

plt.savefig(outfname,dpi=200)
