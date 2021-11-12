function [signal, xps] = resex_mc_signal_from_r(r_fn, xps)
%generates simulated signal using saved positions and waveforms
%returns
% dock_figures(1)
try 
    tmp = load(r_fn, 'positions', 'opt');
    opt = tmp.opt;
    positions = tmp.positions;
    clear tmp
catch
    load(r_fn, 'positions')
    try
        opt = positions.opt;
    catch
        error("Could not find options structure anywhere. Aborting...") 
    end
end


R = positions.R;

%retrieve relevant parameters
gam_ma = msf_const_gamma();
dt = xps.dt;

if ~isfield(opt.sim_opt, 'do_samp'); opt.sim_opt.do_samp = false; end

%get the correct simulation time-step; gwf must match this
if opt.sim_opt.do_avg||opt.sim_opt.do_samp
    true_dt = opt.sim_opt.avg_dt;
else
    true_dt = opt.sim_opt.sim_dt;
end

if dt ~= true_dt; warning("Gradient waveform resolution does not match simulation. We had to resample waveforms from " + num2str(dt*1e3) + " ms to " + num2str(true_dt*1e3) + " ms."); end

gwf = xps.gwf;
b = xps.b;
bvec = xps.bvec;
signal = zeros(xps.n, 1);

wf_counter = 0;
period = numel(bvec);
axs = myfig(1,2);
for c_n = 1:xps.n %for each waveform
        tmp_wf= gwf{c_n};
        %resample if needed
        tmp_wf = gwf_resample_1d(tmp_wf, dt, true_dt);
        %recompute gamma and vomega because they might change a bit from
        %the resampling
        z_q = msf_const_gamma*cumsum(tmp_wf)*true_dt;
        z_b = sum(z_q.^2)*true_dt;
        [xps.Gamma(c_n), xps.Vomega(c_n)] = resex_mc_gwf_to_gamma_vomega(tmp_wf, z_q, z_b, true_dt);

        N_tp = numel(tmp_wf);
        waveform = repmat(tmp_wf, opt.sim_opt.n_dim, 1); %[tmp_wf; tmp_wf];
        N_part = size(R, 1); %number of particles
        N_steps = N_tp; %number of time steps; equal to number of time points in waveform
        sum_g_times_r = zeros(1, N_part);   %initialise phase integral
        for step = 1:N_steps
            sum_g_times_r = sum_g_times_r + waveform(:,step)'*squeeze(R(:,:, step))';
        end
        phase = gam_ma*sum_g_times_r*true_dt; %evaluate phase integral
        tmp_sig = (1/numel(phase))*sum(cos(phase));
        
        signal(c_n) = tmp_sig;
        if mod(c_n, period) == 0
            wf_counter = wf_counter+1;
            where = (c_n-period+1:c_n);
            tmp_b = b(where);
            tmp_sig = signal(where);
            plt = plot(axs{1}, tmp_b/1e9, log(tmp_sig), 'o');
            set(plt, 'MarkerFaceColor', plt.Color)
            tmp_b_fine = linspace(min(tmp_b), max(tmp_b));
            tmp_sig_fine = interp1(tmp_b, log(tmp_sig), tmp_b_fine, 'pchip');
            plot(axs{1}, tmp_b_fine/1e9, tmp_sig_fine, '--', 'Color', plt.Color,...
                'HandleVisibility', 'off')
            %also show gwf
            t = 0:numel(tmp_wf)-1;
            fill(axs{2}, t*1e3, tmp_wf*1e3, plt.Color, 'FaceAlpha', 0.3)
        end
end

legend(axs{1}, array2strings(1:wf_counter))
xlabel(axs{1}, 'b [ms/\mum^2]')
ylabel(axs{1}, 'log(Signal)')
title(axs{1}, 'Simulated signals ')

legend(axs{2}, array2strings(1:wf_counter))
xlabel(axs{2}, 'time [ms]')
ylabel(axs{2}, 'g [mT/m]')
title(axs{2}, 'Gradient waveforms ')

end