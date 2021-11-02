function s_noise = resex_mc_add_noise(s_simulated, snr, N_samp)
%adds Rice distributed noise to signal
%at snr = snr
%generates N_samp realisations of noise
n_sig = numel(s_simulated);
s_noise = zeros(N_samp, n_sig);
for n = 1:N_samp
    noise_real = randn(size(s_simulated))*(1/snr);
    noise_imag = randn(size(s_simulated))*(1/snr);
    tmp_s_noise = sqrt( (s_simulated + noise_real).^2 + noise_imag.^2);
    s_noise(n,:) = tmp_s_noise;
end
end