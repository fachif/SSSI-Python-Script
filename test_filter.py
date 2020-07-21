def test_filter(sig,t,var_thr,step_min):
    '''
    Filter signal Y based on variance. The signals with threshold lower than Thr are removed. 
    The time vector T is required, which has variable step size. The output of the function are:
    Yh- Output signal after filter
    IDX- indexes of the most relevant signals
    VAR-Variance of signals Y.

    If no oscillation is detected, the function gives empty matrices.
    
    '''
    
    import numpy as np
    import pandas as pd
    
    nL = len(sig[0])
    ts = t
    trshold = var_thr
    tsmin = step_min
    
    # Identify disturbance on time-series with variable step
    tk = np.zeros((len(ts),len(ts[0]))) #Generate a zero vector that has the same dimmensions as ts.
    tq = np.zeros((len(ts),len(ts[0]))) #Generate a zero vector that has the same dimmensions as ts.
    tk[1::] = ts[:(len(ts)-1)] #Repositioned all of the time steps by one position.
    tm = ts-tk #This is the time step
    tm[0] = 0
    tq[2:] = tm[1:(len(tm)-1)]
    h = tm-tq
    n = len(tm)
    tr = np.arange(1,(n+1),1) #This looks like a vector that will be used to point the location
    hmax = max(h)
    
    jj = 0
    ff = 0
    ta = np.zeros([100,len(tm),len(tm[0])])
    
    
    for iter in range(0,n-1):
        ta[ff,jj,:] = tm[iter]
        if tm[iter] > float(hmax):
            #qq = qq + 1
            ff = ff + 1
        jj = jj + 1
        
    for iter in range(0,n-1):   
        if int(max(ta[iter,:,:])) == 0:
            t_a = ta[0:iter,:,:]
            break
        
        
    p = len(t_a)
    ll = len(t_a[0])
    jj = 0
    tw = pd.DataFrame()
    tu = pd.DataFrame()
    nm = pd.DataFrame()
    jo = pd.DataFrame()
    yy = 0
    pp = 0
    ww = 0
    gg = 0
    
    for i in range(0,p):
        for k in range(0,ll):
            if t_a[i,k,:] >= 0.0001:
                tw.loc[jj, 0] = t_a[i,k,:]
                tw.loc[jj, 1] = k + 1
                jj = jj + 1
            
            if k == (ll-1):
                nx = len(tw)
                nm.loc[i,0] = nx - gg
                gg = float(nm.sum(0))
    
        for a in range(yy,nx):
            if tw.iloc[a,0] <= 0.2:
                tu.loc[pp,0] = tw.iloc[a,0]
                tu.loc[pp,1] = tw.iloc[a,1]
                pp = pp+1
            
            
            if a == (nx-1):
                jo.loc[i,0] = len(tu) - ww
                ww = float(jo.sum(0))
                
        yy = a + 1
        
    tw = tw.to_numpy()
    tu = tu.to_numpy()
    jo = jo.to_numpy()
    Var = np.zeros([p,len(sig[0])])
    Var_t = np.zeros([p,1])
    
    osc_signal = pd.DataFrame()
    var_s = pd.DataFrame()
    
    for i in range(0,p):
        if float(jo[i]) > 0.001:
            Var[i,:] = np.var(sig[int(tu[int(np.sum(jo[0:i])),1] - 1):(int(tu[int(np.sum(jo[0:i-3])-1),1])+1),:],axis = 0)
            if i == 3:
                Var[i,:] = np.var(sig[int(tu[int(np.sum(jo[0:i])),1]):(int(tu[int(np.sum(jo[0:i+1])-1),1])+1),:],axis = 0)
            Var_t[i,:] = np.sum(Var[i,:])
            if Var_t[i,0] > 0.01:
                time = ts[int(tu[int(np.sum(jo[0:i])),1])-1:int(tu[int(np.sum(jo[0:i-3])-1),1]),:]
                signal = sig[int(tu[int(np.sum(jo[0:i])),1])-1:int(tu[int(np.sum(jo[0:i-3])-1),1]),:]
                
    if time.size == 0:
        osc_signal=pd.DataFrame()
        var_s=[];
        yh=[];
        ydh=[];
        print('Signals not suitable for SSS analysis!')
    else:
        j0=0
        time0=np.zeros([len(time), len(time[0])])
        time0[0]=time[0]
        time0[1:]=time[0:-1]
        dt0 = time - time0
        gg=time.size
        timex = pd.DataFrame()
        signalx = pd.DataFrame()
        dt0x = pd.DataFrame()
        
        for k in range(gg):
            if dt0[k] >= tsmin:
                dt0x.loc[j0,0] = dt0[k]
                timex.loc[j0,0] = time[k]
                signalx[k] = (signal[k,:])
                j0 = j0 + 1
        dt0x = dt0x.to_numpy()
        timex = timex.to_numpy()
        signalx = np.transpose(signalx.to_numpy())       
                
        signalx0 = pd.DataFrame()
        
        if signalx.size == 0:
            print("")
        else:
            for k in range(nL):
                signalx0[k] = (signalx[:,k] - np.mean(signalx[:,k])).tolist()
            
            signalx0 = signalx0.to_numpy()
            ydh = np.concatenate((timex, signalx0), axis=1)
            sig_var = np.transpose(np.var(signalx0,0))
            sig_std = np.transpose(np.std(signalx0,0))
            
            var_sort = np.sort(sig_var)
            var_sort_desc = var_sort[::-1] # THIS IS THE ONE YOU WANT!
            s_idx = np.argsort(sig_var)
            s_idx_desc = s_idx[::-1] # THIS IS THE ONE YOU WANT!
            
            std_sort = np.sort(sig_std)
            std_sort_desc = std_sort[::-1] # THIS IS THE ONE YOU WANT!
            std_idx = np.argsort(sig_std)
            std_idx_desc = std_idx[::-1] # THIS IS THE ONE YOU WANT!
            
            nL0 = len(var_sort_desc)
            vA = np.max(var_sort_desc)
            var_sort0 = var_sort_desc/vA
            
            for i in range(nL):
                if var_sort0[i] >= trshold:
                    osc_signal.loc[i,0] = s_idx_desc[i]
                    var_s.loc[i,0] = var_sort_desc[i]
           
            osc_signal = osc_signal.to_numpy()
            var_s = var_s.to_numpy()
            
    yh = sig[:,np.transpose(osc_signal.astype(int))].reshape(616,10)
    
    return yh, osc_signal, var_s, ydh