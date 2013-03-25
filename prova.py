import sys
import numpy as np
import pylab as pl
from binaryio import load_binary

def check_same(a, nSame):            #function that checks if the number of particles is the same in an nSame number of intervals 
    if len(a) < nSame: return False
    v0 = a[-1]
    j=1
    while j <= nSame and v0 == a[-j]:
        j += 1
    return j == nSame

def find_borders(step, array, stripe):        
    """
    function which finds borders of stripes. It starts from the center of each stripe and 
    it considers progressively larger intervals around the center. It counts the number of
    particles in each interval. When the number of counts remains the same for nSame intervals
    (attached), then it means the stripe has finished. It does it separately for left and right
    borders
    """
    n_particles_l=np.zeros(step-1)          
    n_particles_r=np.zeros(step-1)          
    for i in range(1,step):                 
        n_particles_r[i-1]=np.sum((stripe < array) * (array < stripe+0.001*i))
        if check_same(n_particles_r[:i], nSame):
            right_side=stripe+0.001*(i-nSame)
            break
    for i in range(1,step):
        n_particles_l[i-1]=np.sum((stripe-0.001*i < array) * (array < stripe))
        if check_same(n_particles_l[:i], nSame):
            left_side=stripe-0.001*(i-nSame)
            return left_side, right_side
            break
def find_eps(epsilon, array_x, array_vel, left_side, right_side):    #function which finds what should be the epsilon to assume in order
    """
    """
    for q in epsilon:                        #to have each stripe populated by at least 5*sigma(Nparticles)=
        e=abs(array_vel) < q                    #5*sqrt(Nparticles). This gives a measure of resolution. 
        array_x_2=array_x[e]                    #The function returns None if the stripes can't be resolved.
        N_str=np.sum((left_side< array_x_2)*(array_x_2 <right_side))
        if N_str > 5*np.sqrt(N_str):
            return dict(eps=q, left_side=left_side, right_side=right_side) 
    
def load(fname, eps=0.01):
    D=load_binary(fname)
    X = D['particles']['x'][:, 0]
    V = D['particles']['v'][:, 0]
    w = abs(V) < eps
    X_1 = X[w]
    V_1 = V[w]
    return dict(X=X,V=V,Xeps=X_1,Veps=V_1)

if __name__ == '__main__':
    data_params = [dict(fname='output/rvlr-00096.binary'), 
               dict(fname='output/rvmr-00096.binary')] 

    N=70                                 
    stripe_base = load('output/zvhr-00096.binary')

    #pl.figure(1)
    #pl.hist(X1, bins=2*N)

    X1 = stripe_base['Xeps']
    hist, bin_edges = np.histogram(X1, bins=2*N)    #calculate histogram of particles along X, 140 bins
    i=hist>80                    #take into account only highest peaks
    #print bin_edges[i]
    #print hist[i], bin_edges[i]
    #pl.ylim(0, 7)
    pl.xlim(-1.5,1.5)
    binsize=(np.max(X1)-np.min(X1))/(2*N)
    s=0
    str_cent=np.zeros(len(bin_edges[i]))
    for s in range(len(bin_edges[i])):        #determine the center of each peak or near group of peaks
        str_cent[s]=np.mean(X1[(bin_edges[i][s]+binsize>X1)*(X1>bin_edges[i][s])])
    #pl.show()
    #fn = fname.split('.')[0]
    #pl.savefig()


    right_side=0
    left_side=0

    step=150    #step assumed when considering intervals around stripe center
    nSame = 7    #number of intervals which have to contain the same number of particles
    eps=np.arange(0, 4, 0.001, dtype=np.float)    #list of possible epsilon values to test 


    border_base = load('output/rvhr-00096.binary')

    #for all the stripe centers (derived from zvhr file), it calculates
    #the borders of the stripes as they are in file rvhr. Then it overlaps these
    #borders to file rvlr (low resolution), in order to calculate the number of
    #particles. Then it calulates the epilon needed to have the stripe resolved,
    #that is Npart > 5*sqrt(Npart)    

    for param in data_params:
        sim = load(param['fname'])
        for snum, stripe in enumerate(str_cent,1):                
            left_side, right_side= find_borders(step, border_base['Xeps'], stripe)        
            r = find_eps(eps, sim['X'], sim['V'], left_side, right_side)        
            if r is None: 
                print 'fname: %s  No result.' % param['fname']
            else:
                print 'fname: %s  epsilon: %f  left: %f  right: %f' % (param['fname'], r['eps'], r['left_side'], r['right_side'])
