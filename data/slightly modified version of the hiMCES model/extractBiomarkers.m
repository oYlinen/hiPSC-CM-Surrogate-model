function [featureVectorAP, featureVectorCaT, featureVectorContr] = extractBiomarkers(window, stimFlag)
% Adapted from Forouzandehmehr2021 Model (https://doi.org/10.14814/phy2.15124)

try
    %% 
    t       = window.Time;                % ms 
    i_tot   = window.i_tot;%              % A/F 
    Cai     = window.Cai;                 % mM
    Nai     = window.Nai;                 % mM
    Istim   = window.Istim;               % A/F
    Frc     = window.Frc;                 % mN/mm^2
    Lsrc    = window.Lsrc;                % um
    svl     = window.svl;
    
    V_ode = 1000*window.V_ode;            % V to mV
    dVm = diff(V_ode)./diff(t);
    dVm(end+1) = dVm(end);
    %%  A) dV_dt max, B) mdp
    mdp   = [];
    ap    = [];
    
    for i=2:length(V_ode)-1
        if(V_ode(i)<V_ode(i-1) && V_ode(i)<V_ode(i+1))
            mdp(end+1) = i; 

            idx = find(V_ode(mdp)>0);
            mdp(idx)=[];
        elseif(V_ode(i)>V_ode(i-1) && V_ode(i)>V_ode(i+1))
            ap(end+1) = i;  
        end
    end
    
    %% 
    cl = [];
    for i=1:length(mdp)-1
        cl(end+1) = t(mdp(i+1))-t(mdp(i)); 
    end
    
    %% 

    top = [];
    
    if stimFlag == 0
        try
            % 
            numLastAPForBiomarkers = length(mdp);
            tSoglieSing = [];
            VSoglieSing = [];
            idxSoglieSingole = [];
            alfa = 1;
            for ppp = 1:1:numLastAPForBiomarkers
                [tTOP, TOP, idxTOP] = myDDfeatures(t, V_ode, t(mdp(ppp)), V_ode(mdp(ppp)));
                tSoglieSing(alfa) = tTOP;
                VSoglieSing(alfa) = TOP/1000;
                idxSoglieSingole(alfa) = idxTOP;
                alfa = alfa+1;
            end
            %     
            for ppp = 1: numel(tSoglieSing)
                top(ppp) = find(t == tSoglieSing(ppp));
            end
            
        catch
            % 
            %
            %
            for i=2:length(i_tot)
                if(i_tot(i)<=-0.008 & i_tot(i-1)>=-0.008)
                    top(end+1) = i;
                end
            end
        end
    elseif stimFlag == 1
        for i=2:length(Istim)
            if(Istim(i)>=0.5 & Istim(i-1)<=0.5)
                top(end+1) = i;
            end
        end
    end
    
    figure(5001), hold on, plot(t, V_ode, t(top), V_ode(top), 'or')
    %% APD10 APD30 APD90 RappAPD
    % 
    apd10 = [];
    apd20 = [];
    apd30 = [];
    apd40 = [];
    apd50 = [];
    apd70 = [];
    apd80 = [];
    apd90 = [];
    
    for i=1:length(top)
        if(i<length(mdp) && mdp(i)>top(i))
            % 
            %        
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-10/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-10/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd10(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-20/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-20/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd20(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-30/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-30/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd30(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-40/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-40/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd40(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-50/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-50/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd50(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-70/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-70/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd70(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-80/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-80/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd80(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-90/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i))>=V_ode(top(i))+(1-90/100)*(max(V_ode(top(i):mdp(i)))-V_ode(top(i))),1,'first');
            apd90(end+1) = t(j+top(i))-t(jj+top(i));
        elseif(i+1<=length(mdp) && mdp(i)<top(i))

            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-10/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-10/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd10(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-20/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-20/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd20(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-30/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-30/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd30(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-40/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-40/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd40(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-50/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-50/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd50(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-70/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-70/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd70(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-80/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-80/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd80(end+1) = t(j+top(i))-t(jj+top(i));
            j=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-90/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'last');
            jj=find(V_ode(top(i):mdp(i+1))>=V_ode(top(i))+(1-90/100)*(max(V_ode(top(i):mdp(i+1)))-V_ode(top(i))),1,'first');
            apd90(end+1) = t(j+top(i))-t(jj+top(i));
        end
    end
    RappAPD=(apd30-apd40)./(apd70-apd80);
    
    %% dV/dt max
    [~,LOCS] =  findpeaks(dVm,'MinPeakDistance',sum(apd90)/length(apd90));
    dv_dt = dVm(LOCS);
    
    %% AP features
    APD10    = sum(apd10)/length(apd10);                   % ms
    APD20    = sum(apd20)/length(apd20);                   % ms
    APD30    = sum(apd30)/length(apd30);                   % ms
    APD50    = sum(apd50)/length(apd50);                   % ms
    APD70    = sum(apd70)/length(apd70);                   % ms
    APD90    = sum(apd90)/length(apd90);                   % ms
    RAPP_APD = sum(RappAPD)/length(RappAPD);
    
    dV_dt_max  = sum(dv_dt)/length(dv_dt);                  % V/s
    CL         = sum(cl)/length(cl);                        % ms
    Rate_AP    = 60/(CL/1000);                              % bpm
    MDP        = sum(V_ode(mdp))/length(mdp);               % mV
    APA        = max(V_ode)-MDP;                            % mV
    Peak       = sum(V_ode(ap))/length(ap);
    
    % featureVectorAP = [APA, MDP, CL, dV_dt_max, APD10, APD30, APD90, RAPP_APD, Rate_AP];
    featureVectorAP = [MDP, dV_dt_max, APA, Peak,  APD10, APD20, APD30, APD50, APD70, APD90, Rate_AP, RAPP_APD, CL];
catch
    %     featureVectorAP = NaN(1,9);
    featureVectorAP = NaN(1,13);
end

%% Features Calculation
try
    %% Max min standard
    [max_cai, indmax_cai]=findpeaks(Cai);  %300
    [min_cai, indmin_cai]=findpeaks(-Cai);
    
    %% CL
    CLCa = [];   %ms
    for i=1:length(indmax_cai)-1
        CLCa(end+1) = t(indmax_cai(i+1))-t(indmax_cai(i)); % CL = distanza tra 2 picchi cal
    end
    Freq=1000./CLCa;  %Hz
    
    %% Cai
    % % %         if stimFlag == 0
    bbb = mean(diff(indmax_cai));
    soglia=[1];
    h=1;
    shift=20;%5;
    
    for i=(shift+1):(length(Cai)-shift-1)
        if  (Cai(i+shift)>=(Cai(i)+Cai(i)*0.20))
            %                
            if h==1
                bbb_den = 3;
            else
                bbb_den = 2;
            end
            if i>=soglia(h)+(bbb/bbb_den) 
                soglia(end+1)=i;
                h=h+1;
            end
        end
    end
    figure(3547), plot(t, Cai,'b',t(soglia),Cai(soglia),'m*')
    soglia=soglia(2:end);
    if soglia(1)> indmax_cai(1)
        indmax_cai(1) = [];
        indmin_cai(1) = [];
    end
    
    %% Rise slope/time 10 to 90% and decay slope/time 90 to 10 %
    risetime1090=[];
    riseslope1090=[];
    dectime9010=[];
    decslope9010=[];
    risetime1050 = [];
    Tpeak=[];
    
    dectime9010_mio = [];
    
    for i=1:length(indmax_cai)
        if indmax_cai(i)>soglia(i) && i<=length(soglia) 
            
            j90=find(Cai(soglia(i):indmax_cai(i))>=Cai(soglia(i))+0.9*(Cai(indmax_cai(i))-Cai(soglia(i))),1,'first');
            j10=find(Cai(soglia(i):indmax_cai(i))<=Cai(soglia(i))+0.1*(Cai(indmax_cai(i))-Cai(soglia(i))),1,'last');
            j50=find(Cai(soglia(i):indmax_cai(i))>=Cai(soglia(i))+0.5*(Cai(indmax_cai(i))-Cai(soglia(i))),1,'first');
            risetime1090(end+1)=t(soglia(i)+j90)-t(soglia(i)+j10);
            riseslope1090(end+1)=((Cai(soglia(i)+j90)-Cai(soglia(i)+j10)))/risetime1090(end)*1000*1000;   %nM/s
            risetime1050(end+1)=t(soglia(i)+j50)-t(soglia(i)+j10);
            
        elseif indmax_cai(i)<soglia(i) && i+1<=length(indmax_cai) 
            j90=find(Cai(soglia(i):indmax_cai(i+1))>=Cai(soglia(i))+0.9*(Cai(indmax_cai(i+1))-Cai(soglia(i))),1,'first');
            j10=find(Cai(soglia(i):indmax_cai(i+1))<=Cai(soglia(i))+0.1*(Cai(indmax_cai(i+1))-Cai(soglia(i))),1,'last');
            j50=find(Cai(soglia(i):indmax_cai(i+1))>=Cai(soglia(i))+0.5*(Cai(indmax_cai(i+1))-Cai(soglia(i))),1,'first');
            risetime1090(end+1)=t(soglia(i)+j90)-t(soglia(i)+j10);
            riseslope1090(end+1)=((Cai(soglia(i)+j90)-Cai(soglia(i)+j10)))/risetime1090(end)*1000*1000;
            risetime1050(end+1)=t(soglia(i)+j50)-t(soglia(i)+j10);
        end
        t10_start(i) = t(soglia(i)+j10);
        t50_start(i) = t(soglia(i)+j50);
        Tpeak(i)=t(indmax_cai(i))-t(soglia(i)+j10);
    end
    
    
    
    for i=1:length(indmax_cai)
        if (indmax_cai(i)>indmin_cai(i) && i+1<=length(indmin_cai))
            % 
            j90=find(Cai(indmax_cai(i):indmin_cai(i+1))>=Cai(indmin_cai(i+1))+0.9*(Cai(indmax_cai(i))-Cai(indmin_cai(i+1))),1,'last');
            j10=find(Cai(indmax_cai(i):indmin_cai(i+1))<=Cai(indmin_cai(i+1))+0.1*(Cai(indmax_cai(i))-Cai(indmin_cai(i+1))),1,'first');
            dectime9010(end+1)=t(indmax_cai(i)+j10)-t(indmax_cai(i)+j90);
            decslope9010(end+1)=((Cai(indmax_cai(i)+j10)-Cai(indmax_cai(i)+j90)))/dectime9010(end)*1000*1000;
            
            % 
            j90_mio=find(Cai(indmax_cai(i):soglia(i+1))>=Cai(soglia(i+1))+0.9*(Cai(indmax_cai(i))-Cai(soglia(i+1))),1,'last');
            j10_mio=find(Cai(indmax_cai(i):soglia(i+1))<=Cai(soglia(i+1))+0.1*(Cai(indmax_cai(i))-Cai(soglia(i+1))),1,'first');
            dectime9010_mio(end+1)=t(indmax_cai(i)+j10_mio)-t(indmax_cai(i)+j90_mio);
            
            % 
            jSogliaOtherSide(i) = indmax_cai(i) + find(Cai(indmax_cai(i):indmin_cai(i+1))>=Cai(soglia(i)), 1,'last');
            
        elseif (indmax_cai(i)<indmin_cai(i) && i<=length(indmin_cai))
            % 
            j90=find(Cai(indmax_cai(i):indmin_cai(i))>=Cai(indmin_cai(i))+0.9*(Cai(indmax_cai(i))-Cai(indmin_cai(i))),1,'last');
            j10=find(Cai(indmax_cai(i):indmin_cai(i))<=Cai(indmin_cai(i))+0.1*(Cai(indmax_cai(i))-Cai(indmin_cai(i))),1,'first');
            dectime9010(end+1)=t(indmax_cai(i)+j10)-t(indmax_cai(i)+j90);
            decslope9010(end+1)=((Cai(indmax_cai(i)+j10)-Cai(indmax_cai(i)+j90)))/dectime9010(end)*1000*1000;
            %
            j90_mio=find(Cai(indmax_cai(i):soglia(i))>=Cai(soglia(i))+0.9*(Cai(indmax_cai(i))-Cai(soglia(i))),1,'last');
            j10_mio=find(Cai(indmax_cai(i):soglia(i))<=Cai(soglia(i))+0.1*(Cai(indmax_cai(i))-Cai(soglia(i))),1,'first');
            dectime9010_mio(end+1)=t(indmax_cai(i)+j10_mio)-t(indmax_cai(i)+j90_mio);
            
            % 
            jSogliaOtherSide(i) = indmax_cai(i) + find(Cai(indmax_cai(i):indmin_cai(i+1))>=Cai(soglia(i)), 1,'last');
        end
        %     
        %     
    end
    
    %% new DURATION from Paci2020 paper: https://doi.org/10.1016/j.bpj.2020.03.018
    tSoglia = t(soglia); tSoglia(end) = [];
    CaiSoglia = Cai(soglia); CaiSoglia(end) = [];
    tSogliaOtherSide = t(jSogliaOtherSide);
    CaiSogliaOtherSide = Cai(jSogliaOtherSide);
    DURATION = tSogliaOtherSide-tSoglia;
    
    figure(5000), plot(t, Cai,'b',t(soglia),Cai(soglia),'m*',tSoglia,CaiSoglia,'ko',tSogliaOtherSide,CaiSogliaOtherSide,'k*')
    
    %% CTD 
    CaT_mdp = indmin_cai;
    ctd30 = [];
    ctd50 = [];
    ctd90 = [];
    triangle = [];
    for i=1:length(soglia)
        if(i<length(CaT_mdp) && CaT_mdp(i)>soglia(i))

            j=find(Cai(soglia(i):CaT_mdp(i))>=Cai(soglia(i))+(1-30/100)*(max(Cai(soglia(i):CaT_mdp(i)))-Cai(soglia(i))),1,'last');
            jj=find(Cai(soglia(i):CaT_mdp(i))>=Cai(soglia(i))+(1-30/100)*(max(Cai(soglia(i):CaT_mdp(i)))-Cai(soglia(i))),1,'first');
            ctd30(end+1) = t(j+soglia(i))-t(jj+soglia(i));
            j=find(Cai(soglia(i):CaT_mdp(i))>=Cai(soglia(i))+(1-50/100)*(max(Cai(soglia(i):CaT_mdp(i)))-Cai(soglia(i))),1,'last');
            jj=find(Cai(soglia(i):CaT_mdp(i))>=Cai(soglia(i))+(1-50/100)*(max(Cai(soglia(i):CaT_mdp(i)))-Cai(soglia(i))),1,'first');
            ctd50(end+1) = t(j+soglia(i))-t(jj+soglia(i));
            j=find(Cai(soglia(i):CaT_mdp(i))>=Cai(soglia(i))+(1-90/100)*(max(Cai(soglia(i):CaT_mdp(i)))-Cai(soglia(i))),1,'last');
            jj=find(Cai(soglia(i):CaT_mdp(i))>=Cai(soglia(i))+(1-90/100)*(max(Cai(soglia(i):CaT_mdp(i)))-Cai(soglia(i))),1,'first');
            ctd90(end+1) = t(j+soglia(i))-t(jj+soglia(i));
            %                
            triangle(end+1) = ctd90(i)-ctd30(i);
        elseif(i+1<=length(CaT_mdp) && CaT_mdp(i)<soglia(i))

            j=find(Cai(soglia(i):CaT_mdp(i+1))>=Cai(soglia(i))+(1-30/100)*(max(Cai(soglia(i):CaT_mdp(i+1)))-Cai(soglia(i))),1,'last');
            jj=find(Cai(soglia(i):CaT_mdp(i+1))>=Cai(soglia(i))+(1-30/100)*(max(Cai(soglia(i):CaT_mdp(i+1)))-Cai(soglia(i))),1,'first');
            ctd30(end+1) = t(j+soglia(i))-t(jj+soglia(i));
            %  
            j=find(Cai(soglia(i):CaT_mdp(i+1))>=Cai(soglia(i))+(1-50/100)*(max(Cai(soglia(i):CaT_mdp(i+1)))-Cai(soglia(i))),1,'last');
            jj=find(Cai(soglia(i):CaT_mdp(i+1))>=Cai(soglia(i))+(1-50/100)*(max(Cai(soglia(i):CaT_mdp(i+1)))-Cai(soglia(i))),1,'first');
            ctd50(end+1) = t(j+soglia(i))-t(jj+soglia(i));
            %       
            j=find(Cai(soglia(i):CaT_mdp(i+1))>=Cai(soglia(i))+(1-90/100)*(max(Cai(soglia(i):CaT_mdp(i+1)))-Cai(soglia(i))),1,'last');
            jj=find(Cai(soglia(i):CaT_mdp(i+1))>=Cai(soglia(i))+(1-90/100)*(max(Cai(soglia(i):CaT_mdp(i+1)))-Cai(soglia(i))),1,'first');
            ctd90(end+1) = t(j+soglia(i))-t(jj+soglia(i));
            %       
            triangle(end+1) = ctd90(i)-ctd30(i);
        end
    end
    
    %% Cai features
    DURATION     = mean(DURATION);                         %ms
    RT10Peak     = mean(Tpeak);                            %ms
    RT1050       = mean(risetime1050);                     %ms
    RT1090       = mean(risetime1090);                     %ms
    DT9010       = mean(dectime9010);                      %ms
    Rate_Cai     = mean(Freq);                             %Hz
    Max_Cai      = mean(max_cai);                          %mM
    Min_Cai      = mean(-min_cai);                         %mM
    CTD30        = mean(ctd30);                            %ms
    CTD50        = mean(ctd50);                            %ms
    CTD90        = mean(ctd90);                            %ms
    
    featureVectorCaT = [DURATION, RT10Peak, RT1050, RT1090, DT9010, Rate_Cai, Max_Cai, Min_Cai, CTD30, CTD50, CTD90];
catch
    featureVectorCaT = NaN(1,11);
end

%% features Contractility
try
    idxMaxTension = [];
    idxMinLsrc = [];
    relaxTime50Buf = [];
    
    %% Peak Tension
    for i=2:length(Frc)-1
        if(Frc(i)>Frc(i-1) && Frc(i)>Frc(i+1))
            idxMaxTension(end+1) = i;  
        end
    end
    figure(6000), plot(t, Frc,'b',t(idxMaxTension),Frc(idxMaxTension),'r*')
    
    %% Sarcomere shortening
    for i=2:length(Lsrc)-1
        if(Lsrc(i)<Lsrc(i-1) && Lsrc(i)<Lsrc(i+1))
            idxMinLsrc(end+1) = i; % individuo il minimo di Lsrc
        end
    end
    figure(6001), plot(t, Lsrc,'b',t(idxMinLsrc),Lsrc(idxMinLsrc),'r*')
    
    %% Relaxation time 
    for i=1:length(top)-1
        tSec = t*1e-3;
        %     tSecStartRT50
        %     tSecStopRT50
        relaxTime50Buf(i) = computeContrRT50(tSec, tSec(top(i)), tSec(top(i+1)), Frc);
    end
    
    %% Contractility features
    peakTension     = sum(Frc(idxMaxTension))/length(idxMaxTension); % mN/mm^2
    cellShort       = sum(Lsrc(idxMinLsrc))/length(idxMinLsrc);      % uM
    cellShortPerc   = cellShort*100/max(Lsrc);                       % -
    relaxTime50     = sum(relaxTime50Buf)/length(relaxTime50Buf);    % ms
    
    featureVectorContr = [peakTension cellShortPerc relaxTime50];
catch
    featureVectorContr = NaN(1,3);
end
return;
