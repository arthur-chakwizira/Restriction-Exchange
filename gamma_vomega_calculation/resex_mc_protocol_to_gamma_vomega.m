function [gams, voms] = resex_mc_protocol_to_gamma_vomega(gwf_fn, dt)
%Returns Gamma and Vomega given protocol

if nargin < 2
    load(gwf_fn, 'gwf', 'dt')
else
    gwf = gwf_fn;
end

Ngwf = size(gwf, 1);

gams = zeros(Ngwf, 1);
voms = zeros(Ngwf, 1);

for c_wf = 1:Ngwf
    g = gwf(c_wf, :);
    q = msf_const_gamma*cumsum(g)*dt;
    b = sum(q.^2)*dt;
    
    N = numel(g);
    t = 0:dt:(N-1)*dt;
    
    if ~isequal(size(t), size(g)); t = t'; end
    
    vom = (1/b)*msf_const_gamma()^2*trapz(t, g.^2);
    
    q4 = (1/b^2)*resex_mc_correlate(q.^2, q.^2, dt);
    gam = 2*trapz(t, t.*q4);
    
    gams(c_wf) = gam;
    voms(c_wf) = vom;
    
end
end