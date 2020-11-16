%% showing the simulation process for PSK/DPSK M=16
M = 8;
base_con = exp(1j*(0:(2*pi/M):(2*pi-0.001)).');

N = 1000;
sym = base_con(ceil(M*rand(N, 1)));

N = length(sym);
M = length(base_con);
SNR = 20;
N0 = 1;

% find desired average symbol energy
E_avg = mean(abs(base_con).^2);
E_bav = E_avg / ceil(log2(M));

% using SNR = E_bav / (N0 / 2)
% Ebav = SNR * N0 / 2
E_bav_des = 10^(SNR/20) * N0 / 2;
scaling_factor = sqrt(E_bav_des/E_bav);

% scale base constellation and symbols to true constellation
scaled_con = base_con * scaling_factor;
true_sym = sym * scaling_factor;

variance = N0 / 2;
noise_proc = sqrt(variance/2) * (randn([N, 1]) + 1j*randn([N, 1]));
noisy_transmitted = true_sym + noise_proc;

phase_shift = 2*pi*rand();
noisy_transmitted = noisy_transmitted * exp(1j*phase_shift);

noisy_transmitted = noisy_transmitted(2:N) ./ ...
    (noisy_transmitted(1:N-1) ./ abs(noisy_transmitted(1:N-1)));
true_sym = true_sym(2:N) ./ ...
    (true_sym(1:N-1) ./ abs(true_sym(1:N-1)));

est_sym = l2_nearest(scaled_con, noisy_transmitted);

figure();
tiledlayout(1, 1, 'TileSpacing', 'Compact');
hold on;
scatter(real(scaled_con), imag(scaled_con), 'x');
scatter(real(noisy_transmitted), imag(noisy_transmitted), 'o');
hold off;
xlabel('I');
ylabel('Q');
title('Calculate differentials');
axis equal;
grid on;

exportgraphics(gcf(), 'process_4b.png');

num_errors(est_sym, true_sym) / N