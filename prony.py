def pronyitesla(t,y,XX,tstart,tend,T,shift,tstplot,tedplot,plotFlag):
    
    # %[model]=pronyanalysis(t,y,n,tstart,tend,T,shift,tstplot,tedplot,plotFlag);
    # %Prony analysis program for fitting to a ringdown
    # % Inputs:
    # %   t = time vector (does not need to be equally spaced) -- column vector
    # %       t(1) is assumed to be 0; if not, all time variables are shifted by t(1).
    # %   y = ringdown matrix of order N by NChan corresponding to t, each column is a different signal
    # %   n = order of the estimated model
    # %   tstart = row vector of order 1 by NChan.  tstart(k) is the starting
    # %            time for analysis for y(:,k).
    # %   tend = ending times for analysis (same dimension as tstart)
    # %   shift = flag; if shift = 1, residues are shifted to t = 0.
    # %           If the data is noisy, reccomend shift = 0.
    # %   T = sample period for analysis
    # %   tsttplot = starting time for model simulation
    # %   tedplot = ending time for model simulation
    # %   plotFlag = if = 1, plot results
    # % Outputs (structured array):
    # %   model.Poles = ringdown pole matrix -- column k is for column k of y
    # %   model.Res = ringdown residue matrix
    # %   model.K = ringdown offset row vector
    # %   model.yhat = model signal matrix
    # %   model.that = time vector for yhat (starts at tstart)
    # %
    # % NOTE: It is reccomended that the N/(1+NChan) < 200, where
    # %       N = total number of data points in y to be analyzed, and
    # %       NChan = number of columns in y. 
    
    # % Written by D. Trudnowski, 1999.
    # % Last update, D. Trudnowski, 2005.
    
    # % 1.0 Setup
    
    # %Basic error checks 
    
 
    
    import sys
    import numpy as np
    import math
    import pandas as pd
    import numpy.linalg as lin
    import matplotlib.pyplot as plt
    
    #from scipy.interpolate import CubicSpline
    from scipy import interpolate
    # from scipy.interpolate import UnivariateSpline
    
    # -----------------------------------------------------------------------------------------------------
    # ------------------------------------ BASIC ERROR CHECKS ---------------------------------------------
    # -----------------------------------------------------------------------------------------------------
    
    if ((len(y[0])) != (len(tstart[0]))) or ((len(y[0])) != (len(tend[0]))):
        sys.exit('Dimension error in y, tstart, tend.')
    if ((len(tstart)) != 1) or ((len(tend)) != 1):
        sys.exit('Dimension error in tstart or tend.')
    if (T <= 0) or (tstplot >= tedplot) or (np.max(t) < np.max(np.transpose(tend))) or (np.min(t) > np.min(np.transpose(tstart))):
        sys.exit('Data Error.')
    if (len(t) != len(y)) or (len(t[0]) != 1):
        sys.exit('Dimension error in y or t')
    
    # -----------------------------------------------------------------------------------------------------
    # ----------------------------------- DATA PARAMETERS -------------------------------------------------
    # -----------------------------------------------------------------------------------------------------
        
    NChannels = len(y[0])
    
    # -----------------------------------------------------------------------------------------------------
    # ----------------------------------- SHIFT TIME PARAMETERS -------------------------------------------
    # -----------------------------------------------------------------------------------------------------
    
    tstart = tstart-t[0]                                    # OK!    
    tend = tend-t[0]                                        # OK!
    tstplot = tstplot-t[0]                                  # OK!
    tedplot = tedplot-t[0]                                  # OK!
    t = t-t[0]                                              # OK!
    
    # -----------------------------------------------------------------------------------------------------
    # ---------------------------- Set up analysis data and calculate offset ------------------------------
    # -----------------------------------------------------------------------------------------------------
        
    tanalysis = T*np.transpose(np.array(range(math.ceil(np.max(np.transpose(tend))/T)+2)))     # OK!
    Nstart = np.zeros([1,NChannels])                                                           # OK!
    #Nstart = Nstart.reshape(-1,1)
    Nend = np.zeros([1,NChannels])                                                             # OK!
    #Nend = Nend.reshape(-1,1)
    yanal = np.zeros([len(tanalysis),NChannels])                                               # OK!
    #K = np.zeros([1,NChannels])
    
    # for k=1:NChannels
    #     Nstart(k) = floor(tstart(k)/T)+1;
    #     Nend(k) = ceil(tend(k)/T)+1;
    #     yanal(:,k) = spline(t,y(:,k),tanalysis);
    #     K(1,k) = mean(yanal(Nstart(k):Nend(k),k));
    #     yanal(:,k) = yanal(:,k)-K(k);
    # end
    
    
    K = np.zeros([1,NChannels])                                                                  # OK!
    
    for i in range(NChannels):
        Nstart[0,i] = math.floor(tstart[:,i]/T)
        Nend[0,i] = math.ceil(tend[:,i]/T)+1
        #yanal[:,i] = UnivariateSpline(t.flatten(), y[:,i].flatten())
        #spl = CubicSpline(t.flatten(), y[:,i].flatten())
        spl = interpolate.splrep(t.flatten(), y[:,i].flatten())
        #yanal[:,i] = spl(tanalysis)
        yanal[:,i] = interpolate.splev(tanalysis, spl)
        K[0,i] = np.mean(yanal[int(Nstart[0,i]):int((Nend[0,i])),i])
        yanal[:,i] = yanal[:,i] - K[0,i]
        
        
    # -----------------------------------------------------------------------------------------------------
    # --------------------------------------- FIND MODEL ORDER --------------------------------------------
    # ----------------------------------------------------------------------------------------------------- 
    
    # NdataPoints = Nend-Nstart+1; %Number of data points used for analysis on each channel
    # Ntotal = sum(NdataPoints'); %Total number of data points used for Prony analysis
    # nOrder=XX;
    
    NdataPoints = Nend - Nstart
    Ntotal = np.sum(np.transpose(NdataPoints))
    nOrder = XX
    
    
    
    
    # -----------------------------------------------------------------------------------------------------
    # ------------------------------------- BUILD MATRIX AND VECTOR ---------------------------------------
    # -----------------------------------------------------------------------------------------------------
    
    
    # for k=1:NChannels
    #     Ym = zeros(NdataPoints(k)-nOrder,nOrder);
    #     for kk=1:nOrder
    #         Ym(:,kk) = yanal(Nstart(k)+nOrder-kk:Nstart(k)-kk+NdataPoints(k)-1,k);
    #     end
    #     yv = yanal(Nstart(k)+nOrder:Nstart(k)+NdataPoints(k)-1,k);
    #     if k==1;
    #         Ymatrix = Ym;
    #         yvector = yv;
    #     else
    #         Ymatrix = [Ymatrix;Ym]; %Cancatinate the channels
    #         yvector = [yvector;yv];
    #     end
    # end
    
    for k in range(NChannels):
        Ym = np.zeros([int(NdataPoints[0,k])-nOrder,nOrder])
        for kk in range(nOrder):
            Ym[:,kk] = yanal[(int(Nstart[0,k]) + nOrder-kk-1):(int(Nstart[0,k])-kk+int(NdataPoints[0,k]))-1,k]
        yv = yanal[int(Nstart[0,k])+nOrder:int(Nstart[0,k])+int(NdataPoints[0,k]),k] # + 1
        if k == 0:
            Ymatrix = Ym
            yvector = yv
        else:
            Ymatrix = np.r_[Ymatrix,Ym]
            yvector = np.r_[yvector,yv]
            
    # del Ym, yv, k, kk
    
    yvector = yvector.reshape(-1,1)
    # acoef = np.multiply(lin.pinv(Ymatrix),(yvector))        # ERROR!!!!!!!!!
    # acoef = np.matmul(yvector, lin.pinv(Ymatrix))
   # acoef = np.linalg.lstsq(Ymatrix, yvector, rcond = None)
    acoef = lin.pinv(Ymatrix).dot(yvector)
    
    # -----------------------------------------------------------------------------------------------------
    # ----------------------------------------- FIND POLES ------------------------------------------------
    # -----------------------------------------------------------------------------------------------------
    
    # zPoles = roots([1;-acoef]);
    # sPoles = log(zPoles)/T;        
    
    zPoles = np.roots(np.insert(-acoef,0,1)).reshape(-1,1)
    sPoles = np.log(zPoles).reshape(-1,1)/T
    
    # -----------------------------------------------------------------------------------------------------
    # ------------------------------------ SOLVE FOR RESIDUALS --------------------------------------------
    # -----------------------------------------------------------------------------------------------------
    
    # % 3.0 Solve for residues
    # Res = zeros(nOrder,NChannels);
    # for k=1:NChannels
    #     ZMatrix = zeros(NdataPoints(k),nOrder);
    #     for kk=1:NdataPoints(k)
    #         ZMatrix(kk,:) = (zPoles.').^(kk-1);
    #     end
    #     Res(:,k) = ZMatrix\yanal(Nstart(k):Nend(k),k);
    #     if shift==1; 
    #         Res(:,k) = Res(:,k).*(zPoles.^(-Nstart(k)+1)); %Shift residues to time 0
    #     end
    # end
    # clear k kk
    
    Res = np.zeros([nOrder,NChannels], dtype=complex)
    for k in range(NChannels):
        ZMatrix = np.zeros([int(NdataPoints[0,k]), nOrder], dtype=complex)
        for kk in range(int(NdataPoints[0,k])):
            ZMatrix[kk,:] = np.transpose(zPoles)**kk  #.conjugate()
        # Res[:,k] = np.linalg.lstsq(ZMatrix,yanal[int(Nstart[k]):int(Nend[k]),k])
        Res[:,k] = np.matmul(lin.pinv(ZMatrix),(yanal[int(Nstart[0,k]):int(Nend[0,k]),k]))
        if shift == 1:
            Res[:,k] = (Res[:,k:k+1]*(zPoles**(-Nstart[0,k]))).flatten()
            
    del k, kk
    
    
    # -----------------------------------------------------------------------------------------------------
    # ------------------------------------ REORDER USING ENERGY -------------------------------------------
    # -----------------------------------------------------------------------------------------------------       
            
    # P = zeros(nOrder,NChannels);
    # R = zeros(size(Res));
    # for k=1:NChannels
    #     clear E
    #     for kk=1:nOrder
    #         if abs(real(sPoles(kk)))<1e-8
    #             E(kk)=Res(kk,k)^2*(tend(k)-tstart(k));
    #         else
    #             E(kk)=(Res(kk,k)^2/(2*real(sPoles(kk))))*(exp(2*real(sPoles(kk))*(tend(k)-tstart(k)))-1);
    #         end
    #     end
    #     E=E(:);
    #     [x,ii]=sort(E);
    #     R(:,k)=Res(ii,k); 
    #     P(:,k)=sPoles(ii);
    #     M=[length(ii):-1:1]';
    #     R(:,k)=R(M,k); 
    #     P(:,k)=P(M,k);
    # end
    # clear M k x ii E
    
    P = np.zeros([nOrder,NChannels], dtype = complex)
    R = np.zeros([len(Res), len(Res[0])], dtype = complex)
    for k in range(NChannels):
        E = pd.DataFrame()
        for kk in range(nOrder):
            if np.abs(np.real(sPoles[kk])) < 1e-8:
                E.loc[kk,0] = (Res[kk,k]**2)*(tend[0,k]-tstart[0,k])
            else:
                E.loc[kk,0] = complex((Res[kk,k]**2/(2*np.real(sPoles[kk,0])))*(np.exp(2*np.real(sPoles[kk])*(tend[0,k]-tstart[0,k]))-1))
        
        E = E.to_numpy()
        x = sorted(E, key = abs)
        ii = np.argsort(np.transpose(abs(E)))
        #ii = [i[0] for i in sorted(enumerate(E), key=lambda x:x[1].imag)]
        
        R[:,k] = Res[ii,k]
        P[:,k]=sPoles[ii,0]
        M = np.transpose(np.array(range(len(ii[0]),0,-1)))
        R[:,k] = R[M-1,k]
        P[:,k] = P[M-1,k]
        
    # del M, k, x, ii, E
    
    
    # -----------------------------------------------------------------------------------------------------
    # ----------------------------------------- SIMULATE --------------------------------------------------
    # -----------------------------------------------------------------------------------------------------      
        
    # that = [tstplot:T:tedplot]';
    # yhat = zeros(length(that),NChannels);
    # for k=1:NChannels
    #     yhat(:,k) = K(k).*ones(size(that));
    #     for kk=1:length(that);
    #         for n=1:nOrder;
    #             yhat(kk,k) = yhat(kk,k) + R(n,k)*exp(P(n,k)*that(kk));
    #         end
    #     end;
    # end;
    # yhat=real(yhat);    
        
    that = (np.transpose(np.array(np.arange(tstplot,tedplot,T))).reshape(-1,1)).reshape(-1,1)
    yhat = np.zeros([len(that), NChannels], dtype = complex)
    for k in range(NChannels):
        yhat[:,k] = K[0,k]*np.ones([len(that)])
        for kk in range(len(that)):
            for n in range(nOrder):
                yhat[kk,k] = yhat[kk,k] + R[n,k]*np.exp(P[n,k]*that[kk,0])
    yhat = np.real(yhat)
    
    
    # -----------------------------------------------------------------------------------------------------
    # --------------------------------- Place output in structured array ----------------------------------
    # -----------------------------------------------------------------------------------------------------
    
    # model.Poles = P;
    # model.Res = R;
    # model.K = K;
    # model.that = that;
    # model.yhat = yhat;
    
    # lamda=sPoles;
    
    model_Poles = P
    model_Res = R
    model_K = K
    model_that = that
    model_yhat = yhat
    
    lamda=sPoles
    
    
    # -----------------------------------------------------------------------------------------------------
    # ---------------------------------------- PLOT RESULTS -----------------------------------------------
    # -----------------------------------------------------------------------------------------------------
    
    # % 6.0 Plot results
    # if plotFlag==1
    #        figure
    #        hold on
    #         h1 = plot(t,y);
    #         h2 = plot(that,yhat,'--*');axis tight
    #         hold off
    #         xlabel('Time (sec)')
    #         %legend('Actual','Prony')
    # end
    
    
    if plotFlag == 1:
        plt.figure(2)
        # plt.hold(True)
        plt.plot(t,y)
        #plt.figure(3)
        plt.plot(that, yhat,'--*')
        # plt.hold(False)
        plt.xlabel('Time (sec)')
        # plt.legend('Actual', 'Prony')
        

    
    return lamda, model_Poles, model_Res, model_K, model_that, model_yhat
