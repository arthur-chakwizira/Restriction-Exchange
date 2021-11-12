function xps = resex_mc_xps_build(gwf_fn, b_s, opt, r_fn)
%Rescales waveforms according to the target b-values specified in b_s
%Does this while obeying amplitude constraint in waveform_params
cla

if nargin < 4
    if opt.exp_opt.using_sim_data
        warning("No trajectory file specified. Creating xps based on supplied opt...")
    end
else
    try
        using_sim_data = opt.exp_opt.using_sim_data;
    catch
        using_sim_data = true;
    end
    exp_opt = opt.exp_opt;
    clear opt
    warning off
    tmp = load(r_fn, 'opt');
    try opt = tmp.opt;
    catch
        load(r_fn, 'positions')
        try
            opt = positions.opt;
            clear positions
        catch
            error("Could not find options structure in trajectory file. Aborting...")
        end
    end
    exp_opt.using_sim_data = using_sim_data;
    opt.exp_opt = exp_opt;
end


if isempty(gwf_fn)
    gwf_fn = 'gwf.mat';
end
try
    load(gwf_fn, 'gwf', 'dt');
catch
    errordlg(strcat('Could not find the varialbes ''gwf and dt'' in :', gwf_fn))
end


waveforms = gwf;
[N_wf, ~] = size(waveforms); %number of waveforms, number of time points in each waveform
n =  N_wf*numel(b_s); %total number of acquisitions
b_values = zeros(n, 1); %number of waveforms x how many times we scale it
q_vectors = cell(n,1);%zeros(n, N_tp); %# waveforms x # time points
scaled_wfs = cell(n,1);%zeros(size(q_vectors)); %#waveforms x # time points # number scalings
Gamma = zeros(n,1);
Vomega = zeros(n,1);
gam_ma = msf_const_gamma();
avg_dt = opt.sim_opt.avg_dt; %resolution of saved positions
g_max = opt.exp_opt.g_max;

c_n = 0;

for c_wf = 1:N_wf %for each waveform
    tmp_wf = waveforms(c_wf, :);
    [~, tmp_wf] = N_from_gwf(tmp_wf);
    tmp_wf = g_max*tmp_wf./(max(tmp_wf(:)));
    tmp_q = gam_ma*cumsum(tmp_wf)*dt;
    max_b = sum(tmp_q.^2)*dt;
    if max_b < max(b_s)
        tmp_bs = linspace(min(b_s), max_b, numel(b_s));
        max_b = max(tmp_bs);
    else
        tmp_bs = b_s;
        tmp_wf = sqrt(max(b_s)/max_b)*tmp_wf;
        max_b = max(b_s);
    end
    
    for c_b = 1:numel(b_s)
        c_n = c_n + 1;
        
        b = tmp_bs(c_b);
        
        if opt.exp_opt.using_sim_data
            g_fact = sqrt((b/max_b)/opt.sim_opt.n_dim); %expect waveform will be repeated in each direction
        else
            g_fact = sqrt((b/max_b)/1); %for simulation-independent signal prediction, we do not need to do this scaling
        end
        
        tmp_tmp_wf = g_fact*tmp_wf;
        q_vector = gam_ma*cumsum(tmp_tmp_wf)*dt;
        b_value = sum(q_vector.^2)*dt;
        %Get gamma and vomega
        [gam, vom] = resex_mc_gwf_to_gamma_vomega(tmp_tmp_wf, q_vector, b_value, dt);
        
        %reduce resolution of waveforms to quicken analysis later
        [gwf_new, rf_new] = gwf_subsample_1d(tmp_tmp_wf, dt, avg_dt);
        tmp_tmp_wf = gwf_new.*rf_new;
        
        %         plot(tmp_tmp_wf, 'k', 'LineWidth', 2); hold on; plot(gwf_new.*rf_new, 'r--',  'LineWidth', 2)
        
        scaled_wfs{c_n} = tmp_tmp_wf;
        q_vector = gam_ma*cumsum(tmp_tmp_wf)*avg_dt;
        %         b_value = sum(q_vector.^2)*avg_dt;
        
        b_values(c_n) = b;%b_value;
        q_vectors{c_n} = q_vector;
        Gamma(c_n) = gam;
        Vomega(c_n) = vom;
    end
    
    
end

xps.n = n;
xps.gwf = scaled_wfs; %this is nxN_tp
xps.b = b_values; %
xps.q = q_vectors; %
xps.dt = avg_dt;
xps.Gamma = Gamma;
xps.Vomega = Vomega;
xps.bvec = b_s;

end