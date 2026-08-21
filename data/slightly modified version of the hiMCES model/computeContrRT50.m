function [RT50] = computeContrRT50(t, tStartRT50, tStopRT50, Frc)
minTol = 1e-4;
maxTol = 5e-4;
stepTol = 0.5e-4;
count = 0;
tt= find(tStartRT50<=t & t<=tStopRT50);
ts = tt(1); te = tt(end);
tlc = t(ts:te);                 %final cycle time (1 second at 1Hz)
frct = Frc(ts:te);              % Tension in the cycle
pkt = max(frct);                % peak tension
ppp = find(frct == pkt);
pkt_time = tlc(ppp);            % time at peak tension


while minTol <= maxTol
    %         count;
    minTol = minTol + stepTol*count;
    rel50 = find(abs(frct(ppp:end)-0.5*pkt)<=minTol);
    if numel(rel50)==1
        break
    elseif numel(rel50)>1
        rel50 = rel50(1);
        break
    else
        count = count + 1;
    end
end
if isempty(rel50)
    error('Error in RT50');
end

halfreltime = tlc(ppp+rel50-1);     % time at 50% of relaxation time
RT50 = (halfreltime-pkt_time)*1000; % RT50

end