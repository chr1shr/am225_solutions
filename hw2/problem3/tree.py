import numpy as np
import pylab as plt
import sys

class node:

    def __init__(self,crd,th,nk,lvl):
        self.coord=crd
        self.theta=th
        self.nkids=nk

        n = 5
        low = np.pi/10
        high = np.pi/2
        tpms = np.linspace(np.pi/10,np.pi/2,n)
        if self.nkids <= 2:
            self.tpm=low
        elif self.nkids < n+2:
            self.tpm=tpms[self.nkids-3]
        else:
            self.tpm=tpms[-2]

        self.kids=[]
        self.lvl=lvl

    def from_parent(parent,kid_index,kids):
        theta = parent.theta+(2*((kid_index+1-0.5)/parent.nkids)-1)*parent.tpm
        nkids = kids
        lvl = parent.lvl+1
        coord = (parent.coord[0]+(0.95**lvl)*np.cos(theta),parent.coord[1]+(0.95**lvl)*np.sin(theta))
        return node(coord,theta,nkids,lvl)

    def put(self,nodelabel,nk):
        if nodelabel==0:
            self.kids.append(node.from_parent(self,len(self.kids),nk))
            return -1;
        else:
            for k in self.kids:
                nodelabel = k.put(nodelabel-1,nk)
                if nodelabel<0:
                    return -1
        return nodelabel

    def plot(self,n):
        if self.lvl==0:
            ms=15
        else:
            ms=10
        plt.plot(self.coord[0],self.coord[1],'o',ms=ms,color='black',zorder=10)
#        plt.text(self.coord[0]+0.02,self.coord[1],str(n))
        for k in self.kids:
            n=k.plot(n+1)
            plt.plot([self.coord[0],k.coord[0]],[self.coord[1],k.coord[1]],'-',zorder=0,
                    color=plt.rcParams['axes.prop_cycle'].by_key()['color'][self.lvl],lw=5)
        return n

def plot_tree(tree_list,tbase='',fname=''):

    plt.clf()
    root = node((0,0),-np.pi/2,tree_list[0],0)

    nq=[tree_list[0]]
    lq=[0]
    target=0

    for i in range(1,len(tree_list)):

        root.put(target,tree_list[i])
        if tree_list[i]>0:
            nq.append(tree_list[i])
            lq.append(i)
        else:
            nq[-1] -= 1
            while len(nq)>0 and nq[-1]==0:
               nq.pop()
               lq.pop()
               if len(nq)>0:
                   nq[-1] -= 1
                
        if len(nq)>0:
            target=lq[-1]

    root.plot(0)
    plt.axes().set_aspect('equal','datalim')
    plt.gcf().set_size_inches(5,5)
    plt.axis('off')
#    plt.title(tbase+str(tree_list))
    if (len(fname)==0):
        plt.show()
    else:
        plt.savefig(fname)

#plot_tree([5, 2, 1,1,0, 0,3,0,0,0,0,0,0])
trees = [np.array([[0]])]
max_order = 10
for o in range(2,max_order+1):
    trees.append(np.loadtxt('trees{:d}'.format(o),dtype=int))

# plot order 7 trees
order_to_plot=7
for i in range(len(trees[order_to_plot-1])):
    plot_tree(trees[order_to_plot-1][i],\
#            fname='pics{:d}/{:d}.png'.format(order_to_plot,i)\
            fname='pics{:d}/{:d}_{:d}.pdf'.format(order_to_plot,i%8,i//8)\
            )
