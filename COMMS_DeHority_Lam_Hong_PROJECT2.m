N = 100000;
N0 = 3;
SNR = -4;

randnumb = round(rand(N,1));
bin_ap_m = (randnumb*2)-1;
bin_ap_const = [-1;1];


[true_sym, est_sym] = simulate_transmission(bin_ap_const, bin_ap_m, N0, SNR);
errors = num_errors(true_sym, est_sym);
Perrors_SNR = zeros(1,25);

for i = [1:25]
    SNR = i-5;
    [true_sym, est_sym] = simulate_transmission(bin_ap_const, bin_ap_m, N0, SNR);
    Perror = num_errors(true_sym, est_sym)/N;
    Perrors_SNR(i) = [Perror];
end

figure;
scatter([-4:20],Perrors_SNR);
