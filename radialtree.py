import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib
from matplotlib.lines import Line2D
#import seaborn as sns 
#sns.set_theme()

colormap_list=["nipy_spectral", "terrain","gist_rainbow","CMRmap","coolwarm","gnuplot","gist_stern","brg","rainbow"]

def plot(Z2,fontsize=8,figsize=None, pallete="gist_rainbow", addlabels=True, show=True,sample_classes=None,colorlabels=None,
         colorlabels_legend=None):
    """
    Drawing a radial dendrogram from a scipy dendrogram output.
    Parameters
    ----------
    Z2 : dictionary
        A dictionary returned by scipy.cluster.hierarchy.dendrogram
    addlabels: bool
        A bool to choose if labels are shown.
    fontsize : float
        A float to specify the font size
    figsize : [x, y] array-like
        1D array-like of floats to specify the figure size
    pallete : string
        Matlab colormap name.
    sample_classes : dict
        A dictionary that contains lists of sample subtypes or classes. These classes appear 
        as color labels of each leaf. Colormaps are automatically assigned. Not compatible 
        with options "colorlabels" and "colorlabels_legend".
        e.g., {"color1":["Class1","Class2","Class1","Class3", ....]} 
    colorlabels : dict
        A dictionary to set color labels to leaves. The key is the name of the color label. 
        The value is the list of RGB color codes, each corresponds to the color of a leaf. 
        e.g., {"color1":[[1,0,0,1], ....]}   
    colorlabels_legend : dict
        A nested dictionary to generate the legends of color labels. The key is the name of 
        the color label. The value is a dictionary that has two keys "colors" and "labels". 
        The value of "colors" is the list of RGB color codes, each corresponds to the class of a leaf. 
        e.g., {"color1":{"colors":[[1,0,0,1], ....], "labels":["label1","label2",...]}}   
    
    Returns
    -------
    Raises
    ------
    Notes
    -----
    References
    ----------
    See Also
    --------
    Examples
    --------
    """
    if figsize==None and colorlabels != None:
        figsize=[10,5]
    elif figsize==None and sample_classes != None:
        figsize=[10,5]
    elif figsize==None :
        figsize=[5,5]
    linewidth=0.5
    R=1
    width=R*0.1
    space=R*0.05
    if colorlabels != None:
        offset=width*len(colorlabels)/R+space*(len(colorlabels)-1)/R+0.05
        print(offset)
    elif sample_classes != None:
        offset=width*len(sample_classes)/R+space*(len(sample_classes)-1)/R+0.05
        print(offset)
    else:
        offset=0
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
    labels=[]
    for x, y, c in sorted(zip(Z2['icoord'], Z2['dcoord'],Z2["color_list"]) ):
        
    #x, y = Z2['icoord'][0], Z2['dcoord'][0]
        _color=cmap[ucolors.index(c)]
        if c=="C0": #np.abs(_xr1)<0.000000001 and np.abs(_yr1) <0.000000001:
            _color="black"
        
        # transforming original x coordinates into relative circumference positions and y into radius
        # the rightmost leaf is going to [1, 0] 
        r=R*(1-np.array(y)/ymax)
        _x=np.cos(2*np.pi*np.array([x[0],x[2]])/xmax) # transforming original x coordinates into x circumference positions
        _xr0=_x[0]*r[0]
        _xr1=_x[0]*r[1]
        _xr2=_x[1]*r[2]
        _xr3=_x[1]*r[3]
        _y=np.sin(2*np.pi*np.array([x[0],x[2]])/xmax) # transforming original x coordinates into y circumference positions
        _yr0=_y[0]*r[0]
        _yr1=_y[0]*r[1]
        _yr2=_y[1]*r[2]
        _yr3=_y[1]*r[3]
        #plt.scatter([_xr0, _xr1, _xr2, _xr3],[_yr0, _yr1, _yr2,_yr3], c="b")
        
        
        #if y[0]>0 and y[3]>0:
            #_color="black"
        #plotting radial lines
        plt.plot([_xr0, _xr1], [_yr0, _yr1], c=_color,linewidth=linewidth)
        plt.plot([_xr2, _xr3], [_yr2,_yr3], c=_color,linewidth=linewidth)
        
        #plotting circular links between nodes
        if _yr1> 0 and _yr2>0:
            link=np.sqrt(r[1]**2-np.linspace(_xr1, _xr2, 100)**2)
            plt.plot(np.linspace(_xr1, _xr2, 100), link, c=_color,linewidth=linewidth)
        elif _yr1 <0 and _yr2 <0:
            link=-np.sqrt(r[1]**2-np.linspace(_xr1, _xr2, 100)**2)
            
            plt.plot(np.linspace(_xr1, _xr2, 100), link, c=_color,linewidth=linewidth)
        elif _yr1> 0 and _yr2 < 0:
            _r=r[1]
            if _xr1 <0 or _xr2 <0:
                _r=-_r
            link=np.sqrt(r[1]**2-np.linspace(_xr1, _r, 100)**2)
            plt.plot(np.linspace(_xr1, _r, 100), link, c=_color,linewidth=linewidth)
            link=-np.sqrt(r[1]**2-np.linspace(_r, _xr2, 100)**2)
            plt.plot(np.linspace(_r, _xr2, 100), link, c=_color,linewidth=linewidth)
        
        #Calculating the x, y coordinates and rotation angles of labels
        
        if y[0]==0:
            label_coords.append([(1.05+offset)*_xr0, (1.05+offset)*_yr0,360*x[0]/xmax])
            #plt.text(1.05*_xr0, 1.05*_yr0, Z2['ivl'][i],{'va': 'center'},rotation_mode='anchor', rotation=360*x[0]/xmax)
            #labels.append(label)
            i+=1
        if y[3]==0:
            label_coords.append([(1.05+offset)*_xr3, (1.05+offset)*_yr3,360*x[2]/xmax])
            #plt.text(1.05*_xr3, 1.05*_yr3, Z2['ivl'][i],{'va': 'center'},rotation_mode='anchor', rotation=360*x[2]/xmax)
            #labels.append(label)
            i+=1
    

    if addlabels==True:
        if len(Z2['ivl'])!=len(label_coords):
            print("Warning several labels (samples) may be missing in the tree. This may be due to the data structure you have \
                  (e.g., too many similar samples) or a bug from scipy.")
        
        #Adding labels
        for (_x, _y,_rot), label in zip(label_coords, Z2['ivl']):
            plt.text(_x, _y, label,{'va': 'center'},rotation_mode='anchor', rotation=_rot,fontsize=fontsize)
    
    
    
    if colorlabels != None:
        if len(Z2['ivl'])!=len(label_coords):
            print("Warning several labels (samples) may be missing in the tree. This may be due to the data structure you have \
                  (e.g., too many similar samples) or a bug from scipy.")
        
        j=0
        outerrad=R*1.05+width*len(colorlabels)+space*(len(colorlabels)-1)
        print(outerrad)
        #sort_index=np.argsort(Z2['icoord'])
        #print(sort_index)
        intervals=[]
        for i in range(len(label_coords)):
            _xl,_yl,_rotl =label_coords[i-1]
            _x,_y,_rot =label_coords[i]
            if i==len(label_coords)-1:
                _xr,_yr,_rotr =label_coords[0]
            else:
                _xr,_yr,_rotr =label_coords[i+1]
            d=((_xr-_xl)**2+(_yr-_yl)**2)**0.5
            intervals.append(d)
        colorpos=intervals#np.ones([len(label_coords)])
        labelnames=[]
        for labelname, colorlist in colorlabels.items():
            colorlist=np.array(colorlist)[Z2['leaves']]
            if j!=0:
                outerrad=outerrad-width-space
            innerrad=outerrad-width
            patches, texts =plt.pie(colorpos, colors=colorlist,
                    radius=outerrad,
                    counterclock=True,
                    startangle=label_coords[0][2]*0.5)
            circle=plt.Circle((0,0),innerrad, fc='whitesmoke')
            plt.gca().add_patch(circle)
            labelnames.append(labelname)
            j+=1
        
        if colorlabels_legend!=None:
            for i, labelname in enumerate(labelnames):
                print(colorlabels_legend[labelname]["colors"])
                colorlines=[]
                for c in colorlabels_legend[labelname]["colors"]:
                    colorlines.append(Line2D([0], [0], color=c, lw=4))
                leg=plt.legend(colorlines,
                           colorlabels_legend[labelname]["labels"],
                       bbox_to_anchor=(1.5+0.3*i, 1.0),
                       title=labelname)
                plt.gca().add_artist(leg)   
    elif sample_classes!=None:
        if len(Z2['ivl'])!=len(label_coords):
            print("Warning several labels (samples) may be missing in the tree. This may be due to the data structure you have \
                  (e.g., too many similar samples) or a bug from scipy.")
        j=0
        outerrad=R*1.05+width*len(sample_classes)+space*(len(sample_classes)-1)
        print(outerrad)
        #sort_index=np.argsort(Z2['icoord'])
        #print(sort_index)
        intervals=[]
        for i in range(len(label_coords)):
            _xl,_yl,_rotl =label_coords[i-1]
            _x,_y,_rot =label_coords[i]
            if i==len(label_coords)-1:
                _xr,_yr,_rotr =label_coords[0]
            else:
                _xr,_yr,_rotr =label_coords[i+1]
            d=((_xr-_xl)**2+(_yr-_yl)**2)**0.5
            intervals.append(d)
        colorpos=intervals#np.ones([len(label_coords)])
        labelnames=[]
        colorlabels_legend={}
        for labelname, colorlist in sample_classes.items():
            ucolors=sorted(list(np.unique(colorlist)))
            type_num=len(ucolors)
            _cmp=cm.get_cmap(colormap_list[j], type_num)
            _colorlist=[_cmp(ucolors.index(c)/(type_num-1)) for c in colorlist]
            _colorlist=np.array(_colorlist)[Z2['leaves']]
            if j!=0:
                outerrad=outerrad-width-space
            innerrad=outerrad-width
            patches, texts =plt.pie(colorpos, colors=_colorlist,
                    radius=outerrad,
                    counterclock=True,
                    startangle=label_coords[0][2]*0.5)
            circle=plt.Circle((0,0),innerrad, fc='whitesmoke')
            plt.gca().add_patch(circle)
            labelnames.append(labelname)
            colorlabels_legend[labelname]={}
            colorlabels_legend[labelname]["colors"]=_cmp(np.linspace(0, 1, type_num))
            colorlabels_legend[labelname]["labels"]=ucolors
            j+=1
        
        if colorlabels_legend!=None:
            for i, labelname in enumerate(labelnames):
                print(colorlabels_legend[labelname]["colors"])
                colorlines=[]
                for c in colorlabels_legend[labelname]["colors"]:
                    colorlines.append(Line2D([0], [0], color=c, lw=4))
                leg=plt.legend(colorlines,
                           colorlabels_legend[labelname]["labels"],
                       bbox_to_anchor=(1.5+0.3*i, 1.0),
                       title=labelname)
                plt.gca().add_artist(leg)
            #break
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.bottom.set_visible(False)
    plt.xticks([])
    plt.yticks([])
    if colorlabels!=None:
        maxr=R*1.05+width*len(colorlabels)+space*(len(colorlabels)-1)
    elif sample_classes !=None:
        maxr=R*1.05+width*len(sample_classes)+space*(len(sample_classes)-1)
    else:
        maxr=R*1.05
    plt.xlim(-maxr,maxr)
    plt.ylim(-maxr,maxr)
    if show==True:
        plt.show()
    else:
        return ax

def mat_plot(mat):
    #Take a matrix data instead of a dendrogram data, calculate dendrogram and draw a circular dendrogram
    pass 

def pandas_plot(df):
    
    pass


if __name__=="__main__":
    # Generate random features and distance matrix.
    test=2
    np.random.seed(1)
    numleaf=200
    _alphabets=[chr(i) for i in range(97, 97+24)]
    labels=sorted(["".join(list(np.random.choice(_alphabets, 10))) for i in range(numleaf)])
    x = np.random.rand(numleaf)
    D = np.zeros([numleaf,numleaf])
    for i in range(numleaf):
        for j in range(numleaf):
            D[i,j] = abs(x[i] - x[j])
    if test==1:
        
        #optionally leaves can be labeled by colors
        type_num=12
        _cmp=cm.get_cmap("bwr", type_num)
        _cmp2=cm.get_cmap("hot", type_num)
        colors_dict={"example_color":_cmp(np.random.rand(numleaf)),
                     "example_color2":_cmp2(np.random.rand(numleaf))}
        colors_legends={"example_color":{"colors":_cmp(np.linspace(0, 1, type_num)),
                                         "labels": ["ex1_"+str(i+1) for i in range(type_num)]},
                        "example_color2":{"colors":_cmp2(np.linspace(0, 1, type_num)),
                                          "labels": ["ex2_"+str(i+1) for i in range(type_num)]}}
        #fig = pylab.figure(figsize=(8,8))
        
        # Compute and plot the dendrogram.
        #ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
        Y = sch.linkage(D, method='single')
        Z2 = sch.dendrogram(Y,labels=labels,no_plot=True)
       
        plot(Z2, colorlabels=colors_dict,colorlabels_legend=colors_legends)
    elif test==2:
        Y = sch.linkage(D, method='single')
        Z2 = sch.dendrogram(Y,labels=labels,no_plot=True)
        type_num=6
        type_list=["ex"+str(i) for i in range(type_num)]
        sample_classes={"example_color": [np.random.choice(type_list) for i in range(numleaf)]}
        plot(Z2, sample_classes=sample_classes)
    elif test==0:
        Y = sch.linkage(D, method='single')
        Z2 = sch.dendrogram(Y,labels=labels,no_plot=True)
       
        plot(Z2)
    