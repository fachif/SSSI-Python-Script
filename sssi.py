def sssi(signal,ts,stepsize_min,variance_min,Freq,Damp,Nm):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import pandas as pd
    from prony import pronyitesla
    from test_filter import test_filter
    
    # %Small-signal-stability index (SSI) function 
    # %
    # %[ INDEX ] = sssi( signal,ts,step_min,var_min,F,D,Nm)
    # %
    # % OUTPUTS
    # %
    # %    INDEX - Three layer sss index with smi, ami and gmi.
    # %
    # % INPUTS
    # %   signal - Active power flow of relevant lines 
    # %       ts - Time vector with varable step size 
    # % step_min - Minimum varable step for signal analysis 
    # %  var_min - Filter signals with varianza lower than var_min
    # %       F  - [fmin,fmax] 1x2 vector with frequencies of interest in Hz
    # %       D  - [d0,d1,d2] 1x3 vector with damping ratios to compute index in
    # %                       percent
    # %       Nm - Number of modes used in Prony to reconstruc input signals
    # %
    # %
    # %    * Empty matrices will be delivered if input signals used are not
    # %      suitable for sss analysis
    
    yh,x,l,yhd = test_filter(signal,ts,variance_min,stepsize_min)
    
    if yhd.size == 0:
        out = np.array([])
        sss_smi = np.array([])
        sss_ami = np.array([])
        sss_gmi = np.array([])
        modelb = np.array([])
        detail_l = np.array([])
        detail_f = np.array([])
        detail_d = np.array([])
    else:
        t = yhd[:,0:1]
        y = yhd[:, (x.astype(int))+1]
        y = np.reshape(y,(len(y),len(y[0])))
        
        plt.figure(1)
        plt.plot(t,y)
        
        out_t = t
        out_y = y # WILL NOT USE STRUCT AS MATLAB! WILL NOT USE STRUCT AS MATLAB! WILL NOT USE STRUCT AS MATLAB! WILL NOT USE STRUCT AS MATLAB!
        
    
    # fmin= Freq(1);  % Remove frequencies less than fmin Hz
    # fmax= Freq(2);  % Remove frequencies greater than fmax Hz
    # dmax=0.25; % Remove modes with damping greater than dmax
    # Nm;   % Number of modes used to reconstruct signals using Prony 
    
    # sig=y(:,:);
    # nL=size(sig,2);
    
    # [lamda,modelb]=pronyiTesla(t,sig,Nm,t(1)*ones(1,nL),t(end)*ones(1,nL),0.1,1,t(1),t(end),1);
    
    # Poles.sys=diag(lamda);
    
    # n_modes1 = size(lamda,1); %number of modes
    # fmodes0  = abs(imag(lamda))/2/pi; % frequency of modes in Hz
    # dmodes0  = -cos(atan2(imag(lamda),real(lamda))); %damping of modes
              
    # [fmodes1,mode_idx1]=sort(fmodes0,'ascend'); %sorted frequencies
    # dmodes1=dmodes0(mode_idx1); %sorted dampings
    
    # jj=1; kk=1;    
        
        
        fmin = Freq[0]
        fmax = Freq[1]
        dmax = 0.25
        
        sig = y[:,:]
        nL = len(sig[0])
        
        lamda, model_Poles, model_Res, model_K, model_that, model_yhat = pronyitesla(t, sig, Nm, (t[0])*np.ones([1,nL]), t[len(t)-1]*np.ones([1,nL]), 0.1, 1, t[0], float(t[len(t)-1]), 1)
        

    # n_modes1 = size(lamda,1); %number of modes
    # fmodes0  = abs(imag(lamda))/2/pi; % frequency of modes in Hz
    # dmodes0  = -cos(atan2(imag(lamda),real(lamda))); %damping of modes
              
    # [fmodes1,mode_idx1]=sort(fmodes0,'ascend'); %sorted frequencies
    # dmodes1=dmodes0(mode_idx1); %sorted dampings
    
    # jj=1; kk=1;
    # for i=1:n_modes1,
        
    #     if (fmodes1(i)>=fmin)&(fmodes1(i)<=fmax) %modes within fmin and fmax
    #         fmodes2(jj,1)=fmodes1(i);
    #         dmodes2(jj,1)=dmodes1(i);
    #         mode_idx2(jj,1)=mode_idx1(i);
    #        if dmodes2(jj,1)<=dmax   % discard modes with large damping ratio
    #            fmodes3(kk,1)=fmodes2(jj);
    #            dmodes3(kk,1)=dmodes2(jj);
    #            mode_idx3(kk,1)=mode_idx2(jj);
    #            kk=kk+1;           
    #        end
    #         jj=jj+1;
    #     end
        
    # end
        
        
        dmodes0 = pd.DataFrame()
        
        Poles_sys = np.diag(lamda)
        n_modes1 = len(lamda)
        fmodes0 = np.abs(np.imag(lamda))/(2*np.pi)
        #dmodes0 = -np.cos(math.atan2(np.imag(lamda),np.real(lamda)))
        
        for i in range(n_modes1):
            dmodes0.loc[i,0] = -np.cos(math.atan2(np.imag(lamda[i,0]),np.real(lamda[i,0])))
        
        dmodes0 = dmodes0.to_numpy()
        
        fmodes1 = np.array(sorted(fmodes0))     # OK
        mode_idx1 = np.argsort(np.transpose(fmodes0))   # OK
        
        dmodes1 = dmodes0[mode_idx1]
        
        jj = 0; kk = 0;
        
        fmodes2 = pd.DataFrame()
        dmodes2 = pd.DataFrame()
        mode_idx2 = pd.DataFrame()
        
        fmodes3 = pd.DataFrame()
        dmodes3 = pd.DataFrame()
        mode_idx3 = pd.DataFrame()
        
        for i in range(n_modes1):
            if (fmodes1[i,0] >= fmin) and (fmodes1[i,0] <= fmax):
                fmodes2.loc[jj,0] = fmodes1[i,0]
                dmodes2.loc[jj,0] = dmodes1[0,i,0]
                mode_idx2.loc[jj,0] = mode_idx1[0,i]
                if dmodes2.iloc[jj,0] <= dmax:
                    fmodes3.loc[kk,0] = fmodes2.iloc[jj,0]
                    dmodes3.loc[kk,0] = dmodes2.iloc[jj,0]
                    mode_idx3.loc[kk,0] = mode_idx2.iloc[jj,0]
                    kk = kk + 1
                jj = jj + 1
                
        fmodes2 = fmodes2.to_numpy()
        dmodes2 = dmodes2.to_numpy()
        mode_idx2 = mode_idx2.to_numpy()
        
        fmodes3 = fmodes3.to_numpy()
        dmodes3 = dmodes3.to_numpy()
        mode_idx3 = mode_idx3.astype(int).to_numpy()
        
        
        
        
        
        
        n_modes2 = len(fmodes3)
        jj = np.array(range(0,n_modes2,2)).reshape(-1,1)
        fmodes = fmodes3[jj]
        dmodes = dmodes3[jj,0]
        
        print('\n  f(Hz)        d(%)')
        print('---------------------')
        for i in range(int(n_modes2/2)):
            print('\n {:.4f}       {:.4f}'.format(float(fmodes[i]), float(dmodes[i]*100)))
        print('\n\n')
        
        detail_l = lamda[mode_idx3]
        detail_f = fmodes0[mode_idx3]
        detail_d = dmodes0[mode_idx3]
        
        
        th = pd.DataFrame()
        
        for i in range(len(Damp)):
            th.loc[0,i] = math.acos(Damp[i]/100)
            
        th = th.to_numpy()
        ths = np.pi - th
        
        thl = pd.DataFrame()
        for i in range(len(dmodes)):
            thl.loc[i,0] = math.acos(dmodes[i,0])
        thl = thl.to_numpy()
        thls = np.pi - thl
        
        jj = 0;  kk = len(Damp);
        
        
        SMI = pd.DataFrame()
        
        for i in range(int(n_modes2/2)):
            for m in range(kk):
                SMI.loc[jj,m] = thls[0,i] - ths[0,m]
            jj = jj + 1
         
        SMI = SMI.to_numpy()
            
            
        AMI = pd.DataFrame()
        
        for m in range(kk):
            AMI.loc[0,m] = np.min(SMI[:,m])
        
        AMI = AMI.to_numpy()
        
        GMI = np.min(AMI)
        
        sss_smi = SMI
        sss_ami = AMI
        sss_gmi = GMI
        
        
        
        
        
        
    return  sss_smi, sss_ami, sss_gmi, out_t, out_y, model_Poles, model_Res, model_K, model_that, model_yhat, detail_d, detail_f, detail_l