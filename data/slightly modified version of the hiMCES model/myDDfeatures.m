function [tTOP, TOP, idxTOP] = myDDfeatures(t_short, Vm_short, tMin, VmMin)
VmMin = VmMin*1000;
Vm_short = Vm_short*1000;
% figure(1000), hold on,plot(t_short,Vm_short)
% hold on
% [MDP,MDPind]=min(a(:,2))
MDP = VmMin;
MDPind=find(t_short==tMin);
a = [t_short Vm_short];
MDPtime = tMin;
% MDPtime=a(MDPind,1)
DDendind=min(find(a(:,2)>-10 & a(:,1)>MDPtime));
DDendtime=a(DDendind,1);
linfit=polyfit(a(MDPind:round((DDendind-MDPind)*2/3+MDPind),1),a(MDPind:round((DDendind-MDPind)*2/3+MDPind),2),1);
MDDS=linfit(1);
afit=a(MDPind:DDendind,2);
tfit=a(MDPind:DDendind,1);
NDCC0=20;
t0=DDendtime;
tau0=0.017;
C=linfit(2);
P0=[NDCC0 tau0];
% parott=fminsearch('fcosto',P0);
parott=fminsearch(@(par) fcosto(par, afit, tfit, t0, MDDS, C), P0);
DD= parott(1)*exp(-(t0-tfit)/parott(2))+ MDDS*tfit +C;

%length(DD)
DDduration=DDendtime-MDPtime;
parott(2);
% plot(tfit,afit,tfit,DD,'g--')
TOP=linfit(1)*DDendtime+linfit(2);
LINAMP=TOP-MDP;
NDDC5time=parott(2)*log(0.05)+t0;
NDDCd=NDDC5time-MDPtime;

idxPertTOP = find(t_short>tMin & t_short<DDendtime);
t_shorter = t_short(idxPertTOP);
Vm_shorter = Vm_short(idxPertTOP);
idxTOP = max(find(Vm_shorter<=TOP));
tTOP = t_shorter(idxTOP);

end