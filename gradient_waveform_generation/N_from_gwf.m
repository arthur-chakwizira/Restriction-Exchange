function [N, g] = N_from_gwf(gwf)
%return number of points in gwf (start to end) excluding possibly appended
%zeros
%also returns the cropped waveform

wasrow = false;
if size(gwf,1) < size(gwf,2)
    gwf = gwf';
    wasrow = true;
end


start = (gwf(1) == 0 & gwf(2)==0);
en_d = (gwf(end)==0 & gwf(end-1) == 0);

ind_start = 0;
ind_end = 0;

g = gwf;

if start
    ind_start = find(gwf ~= 0, 1, 'first')-2;
    g(1:ind_start) = [];
end

if en_d
    tmp_gwf = flip(gwf);
    ind_end = find(tmp_gwf ~= 0, 1, 'first') -2;
    g((end-ind_end+1):end) = [];
end

N = numel(gwf) - ind_start - ind_end;

if wasrow; g = g'; end


end