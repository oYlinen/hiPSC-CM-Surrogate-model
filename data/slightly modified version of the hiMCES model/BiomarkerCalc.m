%% Final Biomarker Complete Test
clear; close all; clc
%load hiMCES_cnt_1Hz

Frc = AT; Lsrc = Lsarc; svl = Velo;
window = extractWindowForBiomarkers(Vm,I_tot, Istim ,t,Cai,Nai,Frc,Lsrc,svl, 30);
[featureVectorAP, featureVectorCaT, featureVectorContr] = extractBiomarkers(window, stimFlag) 

% Array guide:
VectorAP = '[MDP, dV_dt_max, APA, Peak,  APD10, APD20, APD30, APD50, APD70, APD90, Rate_AP, RAPP_APD, CL]'
VectorCaT = '[DURATION, RT10Peak, RT1050, RT1090, DT9010, Rate_Cai, Max_Cai, Min_Cai, CTD30, CTD50, CTD90]'
VectorContr = '[peakTension cellShortPerc relaxTime50]'

%close all