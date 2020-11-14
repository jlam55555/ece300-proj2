% simulate_transmission: simulate transmitting symbols with noise, and
% returns the transmitted symbols and the estimates of the transmitted
% symbols after adding noise
%
% params:
% base_con  = base constellation (M x 1)
% sym       = noise-free transmitted symbols (from base_cons, N x 1)
% N0        = noise power
% SNR       = desired SNR per bit (dB)
%
% returns:
% true_sym  = transmitted (scaled) symbols
% est_sym   = estimated symbols
function [true_sym, est_sym, scaling_factor] = simulate_transmission(base_con, sym, N0, SNR)

    N = length(sym);
    M = length(base_con);

    % find desired average symbol energy
    E_avg = mean(abs(base_con).^2);
    E_bav = E_avg / ceil(log2(M));
    
    % using SNR = E_av / (N0 / 2)
    % Eav = SNR * N0 / 2
    % Ebav = SNR * N0 / 2 / log2(M)
    E_bav_des = 10^(SNR/20) * N0 / 2 / ceil(log2(M));
    scaling_factor = sqrt(E_bav_des/E_bav);
    
    % scale base constellation and symbols to true constellation
    scaled_con = base_con * scaling_factor;
    true_sym = sym * scaling_factor;
    
    % produce noise and add to transmitted vectors
    variance = N0 / 2;
    noise_proc = sqrt(variance/2) * (randn([N, 1]) + 1j*randn([N, 1]));
    noisy_transmitted = true_sym + noise_proc;
    
    % determine estimates of transmitted signal
    est_sym = l2_nearest(scaled_con, noisy_transmitted);
end