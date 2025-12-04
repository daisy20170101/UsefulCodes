import seissolxdmf
import numpy as np
import matplotlib.tri as tri


# Rjb calculation and GME loading code

def calc_rjb(x,y,p1=(1719988.,5407437.69),p2=(1783672.647,5452039.)):

    s = (x,y)
    t = ((s[0]-p1[0])*(p2[0]-p1[0])+ (s[1]-p1[1])*(p2[1]-p1[1]))/ ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

    if t < 0 :
        d = np.sqrt((s[0]-p1[0])**2 + (s[1]-p1[1])**2)
    else:
        if t>1 : 
            d = np.sqrt((s[0]-p2[0])**2 + (s[1]-p2[1])**2)
        else:
            c = (p1[0] + t*(p2[0]-p1[0]), p1[1] + t*(p2[1]-p1[1]))
            d = np.sqrt((s[0]-c[0])**2 + (s[1]-c[1])**2)

    return d

p1 = (1800650.75,5504513.0)
p2 = (1888556.30,5415316.55)
# p3 = (1750194.15,5330949.14)
p3 = (1767281.22,5312777.18)
p4 = (1693832.38,5386412.84)


def calc_rjb3D(x,y,p1=(1,1),p2=(1,2),p3=(2,2),p4=(2,1)):
    '''compute the shortest distance to the surfface projection of a fault plane'''

    s = (x,y)
    rjb = np.zeros(4)

    pxy = np.array([p1,p2,p3,p4]) # pxy is (4,2) in shape

    # check if s is inside the rectangle:

    i = 0

    t = ((s[0]-pxy[i,0])*(pxy[i+1,0]-pxy[i,0])+ (s[1]-pxy[i,1])*(pxy[i+1,1]-pxy[i,1]) )/ ((pxy[i,0]-pxy[i+1,0])**2 +(pxy[i+1,1]-pxy[i,1])**2)
    
    idd = 3

    q = ((s[0]-pxy[idd,0])*(pxy[0,0]-pxy[idd,0])+ (s[1]-pxy[idd,1])*(pxy[0,1]-pxy[idd,1]) )/ ((pxy[idd,0]-pxy[0,0])**2 +(pxy[idd,1]-pxy[0,1])**2)

    if (t>1 or t<0 or q<0 or q>1):

    # if not:

        for i in range(3):
            # print('line:',i)

            t = ((s[0]-pxy[i,0])*(pxy[i+1,0]-pxy[i,0])+ (s[1]-pxy[i,1])*(pxy[i+1,1]-pxy[i,1]) )/ ((pxy[i,0]-pxy[i+1,0])**2 +(pxy[i+1,1]-pxy[i,1])**2)

            if t < 0 :
                d = np.sqrt((s[0]-pxy[i,0])**2 + (s[1]-pxy[i,1])**2)
            else:
                if t>1 : 
                    d = np.sqrt((s[0]-pxy[i+1,0])**2 + (s[1]-pxy[i+1,1])**2)
                else:
                    c = (pxy[i,0] + t*(pxy[i+1,0]-pxy[i,0]), pxy[i,1] + t*(pxy[i+1,1]-pxy[i,1]))
                    d = np.sqrt((s[0]-c[0])**2 + (s[1]-c[1])**2)

            rjb[i] = d

        # print('line:',i+1)
        idd = i + 1

        t = ((s[0]-pxy[idd,0])*(pxy[0,0]-pxy[idd,0])+ (s[1]-pxy[idd,1])*(pxy[0,1]-pxy[idd,1]) )/ ((pxy[idd,0]-pxy[0,0])**2 +(pxy[idd,1]-pxy[0,1])**2)

        if t < 0 :
            d = np.sqrt((s[0]-pxy[idd,0])**2 + (s[1]-pxy[idd,1])**2)
        else:
            if t>1 : 
                d = np.sqrt((s[0]-pxy[0,0])**2 + (s[1]-pxy[0,1])**2)
            else:
                c = (pxy[idd,0] + t*(pxy[0,0]-pxy[idd,0]), pxy[idd,1] + t*(pxy[0,1]-pxy[idd,1]))
                d = np.sqrt((s[0]-c[0])**2 + (s[1]-c[1])**2)
        rjb[i+1] = d

        rjb_min = np.min(rjb)
        
    else: 
        rjb_min = 0

    return rjb_min

def load_gm_data(xdmfFilename,origin=(0,0)):
    '''
    the Rjb is calculated assuming the fault plane Ax+By+D=0

    '''
    p1 = (1800650.75,5504513.0)
    p2 = (1888556.30,5415316.55)
    p3 = (1767281.22,5312777.18)
    p4 = (1693832.38,5386412.84)

    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    
    ndt = sx.ReadNdt()
    surfxyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    
    print(surfxyz.shape,connect.shape)
    
    centers = (surfxyz[connect[:,0]] + surfxyz[connect[:,1]] + surfxyz[connect[:,2]])/3.
    depi = np.sqrt( (centers[:,0]-origin[0])**2 + (centers[:,1]-origin[1])**2)
    
    Rjb = []

    for i in range(len(centers[:,1])):
        Rjb.append(calc_rjb3D(centers[i,0],centers[i,1],p1,p2,p3,p4))

    ##%%
    pga = sx.ReadData('PGA')
    pgv = sx.ReadData('PGV')
    sa1 = sx.ReadData('SA01.000s')
    sa3 = sx.ReadData('SA03.000s')
    
    return pga,pgv,sa1,sa3,depi,Rjb
    

def calc_gm_variability(rjb,sa1,nbins=40):

    '''give the model name, number of bins in distance (km) and ground intensity index (either, pgv, sa1, or sa3)'''
    
    # group into bins of Rjb and compute mean and std

    rjb_bins = np.logspace(0,2,nbins)
    mean_bins = np.zeros(nbins-1)
    std_bins = np.zeros(nbins-1)


    for ik, rjb_bin in enumerate(rjb_bins[0:-1]):

        mean_bins[ik] = np.mean(np.log10(sa1[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))
        std_bins[ik]  = np.std(np.log10(sa1[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))

   

    return rjb_bins, mean_bins, std_bins
        
    
    

    