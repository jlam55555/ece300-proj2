% simulate_transmission: simulate transmitting symbols with noise, and
% returns the transmitted symbols and the estimates of the transmitted
% symbols after adding noise
%
% params:
% base_con  = base constellation
% sym       = noise-free transmitted symbols (from base_cons, N x 1)
% N0        = noise power
% SNR       = desired SNR per bit (dB)
%
% returns:
% true_sym  = transmitted (scaled) symbols
% est_sym   = estimated symbols
function [true_sym, est_sym] = simulate_transmission(base_con, sym, N0, SNR)

    % find desired average symbol energy
    % ?????
    
    % scale base constellation and symbols to true constellation
    scaled_con = base_con * ave_sym_energy;
    true_sym = sym * ave_sym_energy;
    
    % produce noise and add to transmitted vectors
    variance = N0/2;
    noise_proc = sqrt(variance/2) * (randn([N, 1]) + 1j*randn([N, 1]));
    noisy_transmitted = true_sym + noise_proc;
    
    % determine estimates of transmitted signal
    est_sym = l2_nearest(scaled_con, noisy_transmitted);
end