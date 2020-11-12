% simulate_transmission_diff: simulate transmitting symbols with noise, and
% returns the transmitted symbols and the estimates of the transmitted
% symbols after adding noise; assumes a differential scheme constellation
% (i.e., evenly spaced throughout a circle) where one constellation point
% lies at theta=0
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
function [true_sym, est_sym, scaling_factor] = simulate_transmission_diff(base_con, sym, N0, SNR)

    N = length(sym);
    M = length(base_con);

    % find desired average symbol energy
    E_avg = mean(abs(base_con).^2);
    E_bav = E_avg / ceil(log2(M));
    E_bav_des = 10^(SNR/20) * N0;
    scaling_factor = sqrt(E_bav_des/E_bav);
    
    % scale base constellation and symbols to true constellation
    scaled_con = base_con * scaling_factor;
    true_sym = sym * scaling_factor;
    
    % produce noise and add to transmitted vectors
    variance = N0/2;
    noise_proc = sqrt(variance/2) * (randn([N, 1]) + 1j*randn([N, 1]));
    noisy_transmitted = true_sym + noise_proc;
    
    % randomly rotate symbols to constant channel phase shift
    phase_shift = 2*pi*rand();
    noisy_transmitted = noisy_transmitted * exp(1j*phase_shift);
    
    % unshift by first element phase to get differentials
    true_sym_phase = true_sym(1)/abs(true_sym(1));
    noisy_sym_phase = noisy_transmitted(1)/abs(noisy_transmitted(1));
    true_sym = true_sym * true_sym_phase';
    noisy_transmitted = noisy_transmitted * noisy_sym_phase';
    
    % determine estimates of transmitted signal
    est_sym = l2_nearest(scaled_con, noisy_transmitted);
end