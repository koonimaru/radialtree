import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib
def plot(Z2,fontsize=8,figsize=[5,5], pallete="gist_rainbow"):
    plt.rcParams['font.family']= 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['svg.fonttype'] = 'none'
    xmax=np.amax(Z2['icoord'])
    ymax=np.amax(Z2['dcoord'])
    
    ucolors=sorted(set(Z2["color_list"]))
    #cmap = cm.gist_rainbow(np.linspace(0, 1, len(ucolors)))
    cmp=cm.get_cmap(pallete, len(ucolors))
    #print(cmp)
    if type(cmp) == matplotlib.colors.LinearSegmentedColormap:
        cmap = cmp(np.linspace(0, 1, len(ucolors)))
    else:
        cmap=cmp.colors
    fig, ax=plt.subplots(figsize=figsize)
    i=0
    label_coords=[]
    for x, y, c in sorted(zip(Z2['icoord'], Z2['dcoord'],Z2["color_list"])):
    #x, y = Z2['icoord'][0], Z2['dcoord'][0]
        r=1-np.array(y)/ymax
        _x=np.cos(2*np.pi*np.array([x[0],x[2]])/xmax)
        _xr0=_x[0]*r[0]
        _xr1=_x[0]*r[1]
        _xr2=_x[1]*r[2]
        _xr3=_x[1]*r[3]
        _y=np.sin(2*np.pi*np.array([x[0],x[2]])/xmax)
        _yr0=_y[0]*r[0]
        _yr1=_y[0]*r[1]
        _yr2=_y[1]*r[2]
        _yr3=_y[1]*r[3]
        #plt.scatter([_xr0, _xr1, _xr2, _xr3],[_yr0, _yr1, _yr2,_yr3], c="b")
        plt.plot([_xr0, _xr1], [_yr0, _yr1], c=cmap[ucolors.index(c)])
        plt.plot([_xr2, _xr3], [_yr2,_yr3], c=cmap[ucolors.index(c)])
        if _yr1> 0 and _yr2>0:
            link=np.sqrt(r[1]**2-np.linspace(_xr1, _xr2, 100)**2)
            plt.plot(np.linspace(_xr1, _xr2, 100), link, c=cmap[ucolors.index(c)])
        elif _yr1 <0 and _yr2 <0:
            link=-np.sqrt(r[1]**2-np.linspace(_xr1, _xr2, 100)**2)
            
            plt.plot(np.linspace(_xr1, _xr2, 100), link, c=cmap[ucolors.index(c)])
        elif _yr1> 0 and _yr2 < 0:
            _r=r[1]
            if _xr1 <0 or _xr2 <0:
                _r=-_r
            link=np.sqrt(r[1]**2-np.linspace(_xr1, _r, 100)**2)
            plt.plot(np.linspace(_xr1, _r, 100), link, c=cmap[ucolors.index(c)])
            link=-np.sqrt(r[1]**2-np.linspace(_r, _xr2, 100)**2)
            plt.plot(np.linspace(_r, _xr2, 100), link, c=cmap[ucolors.index(c)])
        if y[0]==0:
            label_coords.append([1.05*_xr0, 1.05*_yr0,360*x[0]/xmax])
            #plt.text(1.05*_xr0, 1.05*_yr0, Z2['ivl'][i],{'va': 'center'},rotation_mode='anchor', rotation=360*x[0]/xmax)
            i+=1
        if y[3]==0:
            label_coords.append([1.05*_xr3, 1.05*_yr3,360*x[2]/xmax])
            #plt.text(1.05*_xr3, 1.05*_yr3, Z2['ivl'][i],{'va': 'center'},rotation_mode='anchor', rotation=360*x[2]/xmax)
            i+=1
    
    """for i in range(len(label_coords)):
        for j in range(len(label_coords)):
            if i <j:
                _xi,_yi,_roti=label_coords[i]
                _xj,_yj,_rotj=label_coords[j]
                if _xi==_xj and _yi==_yj:
                    print(_xi,_yi,_roti, _xj,_yj,_rotj)"""
    assert len(Z2['ivl'])==len(label_coords), "missmatched label numbers "+str(len(Z2['ivl'])) +" and "+str(len(label_coords))
    
    for (_x, _y,_rot), label in zip(label_coords, Z2['ivl']):
        plt.text(_x, _y, label,{'va': 'center'},rotation_mode='anchor', rotation=_rot,fontsize=fontsize)
            
    #{'ha': 'center', 'va': 'center'}
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.bottom.set_visible(False)
    plt.xticks([])
    plt.yticks([])
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,1.1)
    
    plt.show()
    
if __name__=="__main__":
    # Generate random features and distance matrix.
    np.random.seed(1)
    numleaf=24
    labels=[chr(i)*10 for i in range(97, 97+numleaf)]
    x = np.random.rand(numleaf)
    D = np.zeros([numleaf,numleaf])
    for i in range(numleaf):
        for j in range(numleaf):
            D[i,j] = abs(x[i] - x[j])
    
    #fig = pylab.figure(figsize=(8,8))
    
    # Compute and plot the dendrogram.
    #ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='single')
    Z2 = sch.dendrogram(Y,labels=labels)
    plot(Z2)