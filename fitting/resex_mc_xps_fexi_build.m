function xps = resex_mc_xps_fexi_build(gwf_fn, bs_f, bs_d, opt)
%Rescales gwf according to the target b-values specified in b_f and
%b_d (filtering and detection)
%Does this while obeying amplitude constraint in fexi_opt

load(gwf_fn, 'gwf', 'fexi_opt');

tms = fexi_opt.tm;
[N_wf, ~] = size(gwf); %number of gwf, number of time points in each waveform
n =  N_wf*numel(bs_f); %total number of acquisitions
b_values = zeros(n, 1); %number of gwf x how many times we scale it
b_1 = zeros(n,1);
b_2 = zeros(n,1);
tm = zeros(n,1); %mixing times
q_vectors = cell(n,1);%zeros(n, N_tp); %# gwf x # time points # number scalings
scaled_wfs = cell(n,1);%zeros(size(q_vectors)); %#gwf x # time points # number scalings
gam_ma = msf_const_gamma();
dt = fexi_opt.dt;
avg_dt = opt.sim_opt.avg_dt; %resolution of saved positions

g_max = opt.exp_opt.g_max;

c_n = 0;

Gamma = zeros(n,1);
Vomega = zeros(n,1);

for c_wf = 1:N_wf %for each waveform
    
    tmp_tm = tms(c_wf);
    
    tmp_wf = gwf(c_wf,:);
    [~, tmp_wf] = N_from_gwf(tmp_wf);
    
    %filtering block
    gwf_f = tmp_wf(1:round(end/2));
%     [~, gwf_f] = N_from_gwf(gwf_f);
    gwf_f = g_max*gwf_f./(max(gwf_f(:)));
    q_f = gam_ma*cumsum(gwf_f)*dt;
    max_b_f = sum(q_f.^2)*dt;
    if max_b_f < max(bs_f)
        tmp_bs_f = linspace(min(bs_f), max_b_f, numel(bs_f));
        max_b_f = max(tmp_bs_f);
    else
        tmp_bs_f = bs_f;
        gwf_f = sqrt(max(bs_f)/max_b_f)*gwf_f;
        max_b_f = max(bs_f);
    end
    
    %detection block
    gwf_d = tmp_wf(round(end/2)+1:end);
%     [~, gwf_d] = N_from_gwf(gwf_d);
    gwf_d = g_max*gwf_d./(max(gwf_d(:)));
    q_d = gam_ma*cumsum(gwf_d)*dt;
    max_b_d = sum(q_d.^2)*dt;
    if max_b_d < max(bs_d)
        tmp_bs_d = linspace(min(bs_d), max_b_d, numel(bs_d));
        max_b_d = max(tmp_bs_d);
    else
        tmp_bs_d = bs_d;
        gwf_d = sqrt(max(bs_d)/max_b_d)*gwf_d;
        max_b_d = max(bs_d);
    end
    
    
    
    
    for c_b = 1:numel(bs_f)
        c_n = c_n + 1;
        
        b_d = tmp_bs_d(c_b);
        b_f = tmp_bs_f(c_b);
        
        
         if opt.exp_opt.using_sim_data
            g_fact_d = sqrt((b_d/max_b_d)/opt.sim_opt.n_dim); %expect waveform will be repeated in each direction
            g_fact_f = sqrt((b_f/max_b_f)/opt.sim_opt.n_dim); %expect waveform will be repeated in each direction
         else
            g_fact_d = sqrt((b_d/max_b_d)/1); %expect waveform will be repeated in each direction
            g_fact_f = sqrt((b_f/max_b_f)/1); %expect waveform will be repeated in each direction
         end
        
        %Adding this line below because its needed. AC. 2021-aug-17
        if ~isfinite(g_fact_f); g_fact_f = 0; end
        tmp_gwf_f = g_fact_f*gwf_f;
        tmp_gwf_d = g_fact_d*gwf_d;
        
        tmp_q_f = gam_ma*cumsum(tmp_gwf_f)*dt;
        tmp_q_d = gam_ma*cumsum(tmp_gwf_d)*dt;
        
        tmp_b_1 = sum(tmp_q_f.^2)*dt;
        tmp_b_2 = sum(tmp_q_d.^2)*dt;
        
        
        tmp_tmp_wf = [tmp_gwf_f, tmp_gwf_d];
        q_vector = gam_ma*cumsum(tmp_tmp_wf)*dt;
        b_value = sum(q_vector.^2)*dt;
        
       %Get gamma and vomega
         [gam, vom] = resex_mc_gwf_to_gamma_vomega(tmp_tmp_wf, q_vector, b_value, dt);
        
        %reduce resolution of waveforms to quicken analysis later
        [gwf_new, rf_new] = gwf_subsample_1d(tmp_tmp_wf, dt, avg_dt);
        tmp_tmp_wf = gwf_new.*rf_new;
        
        
        
        scaled_wfs{c_n} = tmp_tmp_wf;
        q_vector = gam_ma*cumsum(tmp_tmp_wf)*avg_dt;
        b_value = sum(q_vector.^2)*avg_dt;
        
        b_values(c_n) = b_value;
        q_vectors{c_n} = q_vector;
        
        
        b_1(c_n) = tmp_b_1;
        b_2(c_n) = tmp_b_2;
        
        tm(c_n) = tmp_tm;
        Gamma(c_n) = gam;
        Vomega(c_n) = vom;
    end
    
    
    
%     attained_max_b = (max_b_d==max(bs_d));
end


xps.gwf = scaled_wfs; %this is N_wf x N_tp x N_b
xps.b = b_values; %N_wf x N_b
xps.q = q_vectors; %N_wf x N_p x N_b
xps.DELTA_1 = fexi_opt.dde.DELTA_1;
xps.b_1 = b_1;
xps.b_2 = b_2;
xps.tm = tm;
xps.n = n;
xps.dt = avg_dt;
xps.Gamma = Gamma;
xps.Vomega = Vomega;
xps.bvec = bs_d;
% if attained_max_b; msg = 'Sir, all gwf attained the maximum b-value';
% else; msg = 'Sir, you need to increase g_max'; end
%     msgbox(msg)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

