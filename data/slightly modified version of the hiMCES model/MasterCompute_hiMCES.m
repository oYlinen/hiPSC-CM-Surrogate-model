
%===============================================================================
% Version Disease Models & Mechanisms (https://doi.org/10.1242/dmm.050365)
% Author: Mohamadamin Forouzandehmehr
%===============================================================================

%%
clear; close all; clc
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);
tic
%% Myofilament original initial state
Nxb = 0.9997;  % Non-permissive fraction in overlap region
XBpreR = 0.0001;   % Fraction of XBs in pre-rotated state (Stiffness generation)
XBpostR = 0.0001;  % Fraction of XBS in post-rotated state (force generation)
x_XBpreR = 0;    % Strain of XBs in the pre-rotated state
x_XBpostR = 0.007;  % Strain of XBs in the post-rotated state
TropCaL = 0.0144; % Ca bound to the low affinity troponin site
TropCaH = 0.1276; % Ca bound to the high affinity troponin sit
IntegF = 0;
SL = 1.9; % Sarcomere length

init = [Nxb XBpreR XBpostR x_XBpreR x_XBpostR TropCaL TropCaH IntegF SL];

% Y=[-0.070  0.32    0.0002  0    0    1     1     1      0      1      0   0.75  0.75  0   0.1    1    0    9.2    0     0.75    0.3     0.9     0.1    0.97    0.01    0.01    0.01    1.9     0.0070  0   0   0   0   0   0];
% Yn = Y(1:23);
% Yn(24:32) = init;
% Yn(33) = [0.133];

load ('Yss.mat');
%Yn(33) = [0.1323];

%% Levosimendan simulation
tDrugApplication = 1400;        % it is 400 for Levo ischemia-reperfusion simulations
dose_uM = 0.3;
fkon = 1.1;                     %fkon for 0.3 uM Levo: 1.1, 1 uM Levo: 1.21, 2 uM Levo: 1.32, 3 uM Levo: 1.39, 5 uM Levo: 1.57, 10 uM: 1.88, 15 uM: 2.2, 20 uM: 2.51.
flv =  1.002385;                %fKATP for different concentrations of Levosimendan. 10 uM Levo = 1.843903, 2 uM = 1.118538, 0.3 uM = 1.002385. 

INaFIC50 = 85.714; % uM
INaFHill = 1;
INaFRedMed = 1./(((dose_uM)./(INaFIC50)).^INaFHill+1);

ICaLIC50 = 28.45; % uM
ICaLHill = 1;   
ICaLRedMed = 1./(((dose_uM)./(ICaLIC50)).^ICaLHill+1);

IKrIC50 = 22; % uM
IKrHill = 0.68;
IKrRedMed = 1./(((dose_uM)./(IKrIC50)).^IKrHill+1);

INaLRedMed = 1;
IKsRedMed  = 1;
INaCaRedMed = 1;
IfRedMed = 1;

%% To provide SET version of the model
% load('params_and_states_7_EAD_models.mat')
% clearvars initialStates
% latinHypercubeVector = parameters(5,:);
 myLatinHypercubeVector = ones(1,23);
% myLatinHypercubeVector(1:12) = latinHypercubeVector; %used for hiMCES_SET

%% 0: Temperature = 37°C , 1: Temperature = 21°C
classicVSoptic = 0;

%% 0: spontaneous beating, 1: paced.
stimFlag = 0;
% Set the twitch protocol
contrflag = 1;     % 0 = isometric  1 = active contraction

%% ischemia setting
% For hiMCES_SET version (arrhythmic ischemia-reperfusion simulations) ischemiaFlag must be 1.

ischemiaFlag = 0; % 0 = no ischemia, 1 = SEV1, 2 = SEV2 
myFlags = zeros(1,7); % 1: Hyperkalemia OFF:0 / ON:1, 2: Acidosis OFF:0 / ON:1, 3: Hypoxia IKatp OFF:0 / ON:1, 4: Hypoxia on INaK OFF:0 / ON:1, 5: Hypoxia on IpCa OFF:0 / ON:1, 6: Hypoxia on Jup OFF:0 / ON:1, 7: Hypoxia on CE OFF:0 / ON:1.
myFlags = [1  1  1  1  1  1  1];
IKatpSelectionString = 'Kazbanov'; % it should take a string: Ledezma or Kazbanov

syms x
% z takes the extent of inhibition of INaK & IpCa and calculates
% corresponding normalized O2s
% fnp takes the level of oxygantion (0 to 1) and calculates the corresponding
% inhibition for INaK & IpCa

z = 1; fnp = 0;                   % Physoxia

if ischemiaFlag == 1
z = solve (1-rho(x) == 0.2, x);
z = double(z);
fnp = 1 - rho(z);

else if ischemiaFlag == 2
z = solve (1-rho(x) == 0.31, x);
z = double(z);
fnp = 1 - rho(z);

end
end

%% Integration
% @hiMCES is the standard model
% @hiMCES_SET is the SET version used for arrhythmic ischemia-reperfusion simulations where all ischemic changes (table 1 manus)
% occur following ouabain-like curve, i.e., it considers all changes based on rho
% (Eq. 15).

[t,Yc] = ode15s(@hiMCES,[0 800],Yn, options, contrflag, myLatinHypercubeVector(2:end), stimFlag, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, fkon, flv);
Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];
caSR = Yc(:,2);
Cai  = Yc(:,3);
Nai  = Yc(:,18);


for i= 1:size(Yc,1)
[~, dati]    = hiMCES(t(i), Yc(i,:), contrflag, myLatinHypercubeVector(2:end), stimFlag, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, fkon, flv);
    INa(i)   = dati(1);
    If(i)    = dati(2);
    ICaL(i)   = dati(3);
    Ito(i)   = dati(4);
    IKs(i)   = dati(5);
    IKr(i)   = dati(6);
    IK1(i)   = dati(7);
    INaCa(i) = dati(8);
    INaK(i)  = dati(9);
    IpCa(i)  = dati(10);
    IbNa(i)  = dati(11);
    IbCa(i)  = dati(12);
    Irel(i)  = dati(13);
    Iup(i)   = dati(14);
    Ileak(i) = dati(15); 
    Istim(i) = dati(16);
    E_K(i)   = dati(17);
    E_Na(i)  = dati(18);
    INaL(i)  = dati(19);
    Fps(i)  = dati(20);
    AT(i)  = dati(21);
    ATPase(i)  = dati(22);
    Lsarc(i) = dati(23);
    Velo(i) = dati(24);
    JCB(i) = dati(25);
    PSP(i) = dati(26);          % SERCA phosphorylation
    IKATP(i) = dati(27);
    ronak(i) = dati(28);         
    O2s(i) = dati(29);          % Bath/source O2 concentration   
    O2e(i) = dati(30);          % Extracellular O2 concentration
    vcy(i) = dati(31);          % SERCA cycle rate
    do2dt(i) = dati(32);
    fhib(i) = dati(33);         % Returns finhib
    Koo(i) = dati(34);
    ro1(i) = dati(35);
    ro2(i) = dati(36);
    fkatp(i) = dati(37);
    kon(i) = dati(38);
end
result       = [INa; If; ICaL; Ito; IKs; IKr; IK1; INaCa; INaK; IpCa; IbNa; IbCa; Irel; Iup; Ileak; Istim; E_K; E_Na; INaL];
mat_currents = [INa; If; ICaL; Ito; IKs; IKr; IK1; INaCa; INaK; IpCa; IbNa; IbCa; Irel; Iup; Ileak; Istim; INaL];
I_tot=sum(mat_currents);

toc
% End
