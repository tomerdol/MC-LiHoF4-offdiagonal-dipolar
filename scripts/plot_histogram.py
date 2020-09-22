import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import glob

def main():
    path='./data/lattice_output/res_test/*'
    fig, ax = plt.subplots()
    allBx=[]
    for fname in glob.glob(path):
        data = np.genfromtxt(fname,comments='#',skip_header=10)
        allBx.append(data[:,3].tolist())
    ax.hist(np.array(allBx).flatten())
    fig.savefig('./figures/foo.png')

if __name__ == "__main__":
    main()