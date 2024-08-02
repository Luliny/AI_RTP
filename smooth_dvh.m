function out_dvh = smooth_dvh(in_dvh,xaxis_dvh)

%% smooth DVHs.
%% L.Yuan 09/14/2011

work_dvh = in_dvh;

%%xaxis_dvh = [0: 0.05:4.8]/4;

ind_end = numel(xaxis_dvh);

ind = 1;
while ind <= ind_end-1
    ind2 = ind+1;
    if in_dvh(ind)<in_dvh(ind2)-0.0001
    while ind2 <= ind_end-1 && in_dvh(ind)<in_dvh(ind2)-0.0001
       	  ind2 = ind2 +1;
    end

yi = interp1([xaxis_dvh(1:ind) xaxis_dvh(ind2:end)],[work_dvh(1:ind) work_dvh(ind2:end)],xaxis_dvh(ind+1:ind2),'pchip');
work_dvh(ind+1:ind2) = yi;
end

ind = ind2;
end

out_dvh = max(work_dvh,0.);

return
	    
