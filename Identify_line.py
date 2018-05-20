#Sample of code from finding lines in microtubule images project.
#The main program, Find_lines would be used several times on several different windows of the image.

def smooth0(y, box_pts): #convolve y by box_pts length array of equal length
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode = 'same')
    return y_smooth
def smooth(y, box_pts): # Smooths twice
    box = np.ones(30)/30.0
    y = y - np.convolve(y, box, mode = 'same')
    box = np.ones(box_pts)/box_pts
    return np.convolve(y, box, mode = 'same')

#Create arrays find_lines_ind and find_lines_sum_ind that are used in main program Find_line(a)

n_angle = 80 #number of angles
n = 60 #size of window

m = n #Keeping window square
c=[]
y=[]
middle=int(m/np.sqrt(2))-12
y = .5*np.arange(-middle-5,middle+6)
len_y = middle*2+11 
X = np.tile(np.repeat([np.arange(m)-middle],n,axis=0),[n_angle,len_y,1,1])
Y = np.tile(np.repeat([np.arange(n)-middle],m,axis=0).transpose(),[n_angle,len_y,1,1])

M = np.swapaxes(np.tile(np.tan((np.arange(n_angle)*1.0/n_angle-0.5)*np.pi)+0.01,[n,len_y,n,1]),0,3)
i = np.swapaxes(np.tile(np.arange(-middle-5,middle+6),[n_angle,n,n,1]),1,3)
b = np.sqrt(1+M*M)*i*0.5*np.sign(M)
ind = (np.abs(b+M*X-Y)/np.sqrt(1+M*M)<1)
b = np.sqrt(225*(1+1/(M*M)))
ind2 = (-1/M*X+b>Y)*(-1/M*X-b<Y)
find_lines_ind = ind*ind2 #binary array, indicates where the square is for each angle 
find_lines_sum_ind = np.sum(find_lines_ind,axis=(2,3))


def Find_line(a):  #Identifies location and slope of line in image (end point identification is
#     not included as it still needs to be optimized)

    gauss_ratio_0 = np.sqrt(np.log(0.5)/np.log(0.75))
    gauss_ratio_1 = np.sqrt(np.log(0.5)/np.log(0.8))
    def Analyze_peak(c,y,disp): 
        if(len(c)>7):
            Max_ind = np.argmax(c)
            thresh = np.max(c)*0.5
            Var = np.sum(c>thresh)
            if(c[0] < thresh and c[-1] < thresh):
                return [Var,Max_ind+disp]
            else:
                thresh = np.max(c)*0.75
                Var = np.sum(c>thresh)*gauss_ratio_0
                if(c[0] < thresh and c[-1] < thresh):
                    return [Var,Max_ind+disp]
                else:
                    thresh = np.max(c)*0.8
                    Var = np.sum(c>thresh)*gauss_ratio_1
                    if(c[0] < thresh and c[-1] < thresh):
                        return [Var,Max_ind+disp]
                    else:
                        return [np.inf,Max_ind+disp]
        else:
            return [np.inf,disp]

    M = (np.arange(80)/80.0-0.5)*np.pi+0.001
    Var_0,Max_ind_0,Max_ind_lst,Slope,Peak_slope,Peak_var,Peak_max_ind,Cnt_0 = [],[],[],[],[],[],[],[]
    C = np.sum(find_lines_ind*a,axis=(2,3))/find_lines_sum_ind
    y = .5*np.arange(-middle-5,middle+6)
    for m_i,m in enumerate(M):
        c =C[m_i]
        #Analyze plots
        Analy = []
        peak_ind = signal.argrelmin(smooth(c,11)[5:-5])[0]+5
        l = len(peak_ind)
        if(l==0):
            Analy.append(Analyze_peak(c[5:-5],y[5:-5],5))
        else:
            Analy.append(Analyze_peak(c[5:peak_ind[0]],y[5:peak_ind[0]],5))
            Analy.append(Analyze_peak(c[peak_ind[-1]:-5],y[peak_ind[-1]:-5],peak_ind[-1]))
            for i in xrange(l-1):
                Analy.append(Analyze_peak(c[peak_ind[i]:peak_ind[i+1]],y[peak_ind[i]:peak_ind[i+1]],peak_ind[i]))
        Var,Max_ind = np.transpose(Analy)
        ind = np.where(Var<50)[0]
        Var = Var[ind]
        Max_ind = Max_ind[ind]
        peaks = []
        i = 0
        Max_ind = list(Max_ind)
        Var = list(Var)
        Max_ind_1 = Max_ind_0
        for i2 in xrange(len(Max_ind_0)):
            j = Max_ind_0[i]
            if(j in Max_ind):
                ind = Max_ind.index(j)
            elif(j+1 in Max_ind):
                ind = Max_ind.index(j+1)
            elif(j-1 in Max_ind):
                ind = Max_ind.index(j-1)
            elif(j+2 in Max_ind):
                ind = Max_ind.index(j+2)
            elif(j-2 in Max_ind):
                ind = Max_ind.index(j-2)
            else:
                ind = -1
            if(ind != -1):
                Cnt_0[i] = 0
                Var_0[i].append(Var.pop(ind))
                Max_ind_0[i] = Max_ind[ind]
                Max_ind_lst[i].append(Max_ind.pop(ind))
                Slope[i].append(m)
                i = i+1
            else:
                if(Cnt_0[i] < 3):# and np.abs(np.tan(Slope[i][-1])-np.tan(Slope[i][0]))<2):
                    Cnt_0[i] = Cnt_0[i]+1
                    i = i+1
                else:
                    Peak_var.append(Var_0.pop(i))
                    Peak_max_ind.append(Max_ind_lst.pop(i))
                    Peak_slope.append(Slope.pop(i))
                    Cnt_0.pop(i)
                    Max_ind_0.pop(i)
        
        for i in xrange(len(Max_ind)):
            Cnt_0.append(0)
            Max_ind_0.append(Max_ind[i])
            Var_0.append([Var[i]])
            Max_ind_lst.append([Max_ind[i]])
            Slope.append([m])
    #Merge chains
    
    Merged_ind = []
    Merged_dir = []
    for i,s in enumerate(Peak_slope):
     for t in xrange(3):
        if(M[t] in s):
            j = Peak_max_ind[i][0]
            if(j in Max_ind_0):
                ind = Max_ind_0.index(j)
            elif(j+1 in Max_ind_0):
                ind = Max_ind_0.index(j+1)
            elif(j-1 in Max_ind_0):
                ind = Max_ind_0.index(j-1)
            elif(j+2 in Max_ind_0):
                ind = Max_ind_0.index(j+2)
            elif(j-2 in Max_ind_0):
                ind = Max_ind_0.index(j-2)
            else:
                ind = -1
            if(ind != -1 and Cnt_0[ind]+t<3):
                Merged_ind.append(i)
                Merged_dir.append(len(Peak_slope[i])>len(Slope[ind]))
                Peak_slope[i] = Slope.pop(ind)+Peak_slope[i]
                Peak_max_ind[i] = Max_ind_lst.pop(ind)+Peak_max_ind[i]
                Peak_var[i] = Var_0.pop(ind)+Peak_var[i]
                Max_ind_0.pop(ind)
                Cnt_0.pop(ind)             
    for i in xrange(len(Var_0)):
        Peak_var.append(Var_0.pop(0))
        Peak_max_ind.append(Max_ind_lst.pop(0))
        Peak_slope.append(Slope.pop(0))
    Theta = []
    Max_ind = []
    certinity = []
    i2=0
    
    
    for i,m in enumerate(Peak_slope):
        var = np.array(Peak_var[i])
        Min = np.min(var)
        l = len(var)
        if(l>3): # and l-np.argmin(var)>2 and np.argmin(var)>2):
            clr = Color[i2]
            i2=i2+1
            var = smooth0(var,3)[1:-1]
            x = m[1:-1]
            x = np.array(x)
            if(i in Merged_ind):
                j = Merged_ind.index(i)
                if(Merged_dir[j]):
                    k = np.argmax(x)
                    x[:k+1]= x[:k+1]-np.pi
                else:
                    k = np.argmin(x)
                    x[k:]= x[k:]+np.pi
            var = np.exp((Min*1.1-var)**3)
            if(np.sum(var) != 0):
                Theta.append(np.sum(var*x)/np.sum(var))
                Max_ind.append(np.sum(var*np.array(Peak_max_ind[i])[1:-1])/np.sum(var))
