% number of samples
N = 100000;

% noise power
N0 = 1;

tic();

%% binary antipodal

% binary antipodal constellation and random signal sequence
bin_ap_const = [-1; 1];
bin_ap_sym = bin_ap_const(2 * ceil(rand(N, 1)));

% loop through SNRs from -4 to 20
Perrors_ap = zeros(1, 25);
Perrors_theo_ap = zeros(1, 25);
for i = 1:25
    SNR = i-5;
    
    [true_sym_ap, est_sym_ap, root_Eb] ...
        = simulate_transmission(bin_ap_const, bin_ap_sym, N0, SNR);
    Perror = num_errors(true_sym_ap, est_sym_ap) / N;
    Perrors_ap(i) = Perror;
    
    % see note about fudge factor
    Perrors_theo_ap(i) = qfunc(root_Eb*sqrt(4/N0));
    
end

fprintf('Elapsed for binary antipodal for N=%d: %.2fs\n', N, toc());
tic();

figure('visible', 'off', 'position', [0 0 1000 500]);
tiledlayout(1, 1, 'TileSpacing', 'Compact');
nexttile();
hold on;
scatter(-4:20, Perrors_ap, 'b');
plot(-4:20, Perrors_theo_ap, 'b');

% binary orthogonal
bi_ortho_cons = [1; 1j];
bi_ortho_m = bi_ortho_cons(2 * ceil(rand(N, 1)));

Perrors_orth = zeros(1,25);
Perrors_theo_orth = zeros(1,25);
for i = 1:25
    SNR = i-5;
    [true_sym_ortho, est_sym_ortho, root_Eb] ...
        = simulate_transmission(bi_ortho_cons, bi_ortho_m, N0, SNR);
    Perror = num_errors(true_sym_ortho, est_sym_ortho)/N;
    Perrors_orth(i) = Perror;
    
    % see note about fudge factor
    Perrors_theo_orth(i) = qfunc(root_Eb*sqrt(2)/sqrt(N0));
end

scatter(-4:20, Perrors_orth, 'r');
plot(-4:20, Perrors_theo_orth, 'r');
hold off;

title('Binary PAM bit error rate vs. SNR');
ylabel('P_{error}');
xlabel('SNR (dB)');
legend([
    "binary antipodal empirical", "binary antipodal theoretical", ...
    "binary orthogonal empirical", "binary orthogonal theoretical" ...
]);
exportgraphics(gcf(), 'binary.eps');

fprintf('Elapsed for binary orthogonal for N=%d: %.2fs\n', N, toc());
tic();

%% PSK

Ms = [4 8 16 32];
colors = ["b", "r", "g", "k"];
figure('visible', 'off', 'position', [0 0 1000 500]);
tiledlayout(1, 1, 'TileSpacing', 'Compact');
nexttile();
hold on;
for i = 1:length(Ms)
    M = Ms(i);
    
    % generate basis
    theta = linspace(0, 2*pi, M+1);
    theta = theta(1:M).';
    PSK_cons = cos(theta) + 1j*sin(theta);

    % generate signal sequence
    PSK_sym = PSK_cons(ceil(M * rand(N, 1)));

    Perrors_PSK = zeros(1,25);
    Perrors_theo_PSK = zeros(1,25);
    for j = 1:25
        SNR = j-5;
        [true_sym_psk, est_sym_psk, root_Eb] ...
            = simulate_transmission(PSK_cons, PSK_sym, N0, SNR);
        Perror = num_errors(true_sym_psk, est_sym_psk)/N;
        Perrors_PSK(j) = Perror;
        
        % see note about fudge factor
        Perrors_theo_PSK(j) = 2*qfunc(root_Eb*sqrt(4/N0)*sin(pi/M));
    end

    scatter(-4:20, Perrors_PSK, colors(i));
    plot(-4:20, Perrors_theo_PSK, colors(i));
end

hold off;
title('PSK bit error rate vs. SNR');
ylabel('P_{error}');
xlabel('SNR (dB)');
legend([
    "M=4 empirical", "M=4 theoretical", ...
    "M=8 empirical", "M=8 theoretical", ...
    "M=16 empirical", "M=16 theoretical", ...
    "M=32 empirical", "M=32 theoretical"
]);
exportgraphics(gcf(), 'psk.eps');

fprintf('Elapsed for PSK for N=%d: %.2fs\n', N, toc());
tic();

%% DPSK

Ms = [4 8 16 32];
colors = ["b", "r", "g", "k"];
figure('visible', 'off', 'position', [0 0 1000 500]);
tiledlayout(1, 1, 'TileSpacing', 'Compact');
nexttile();
hold on;
for i = 1:length(Ms)
    M = Ms(i);
    
    % generate basis
    theta = linspace(0, 2*pi, M+1);
    theta = theta(1:M).';
    DPSK_cons = cos(theta) + 1j*sin(theta);

    % generate signal sequence
    PSK_sym = DPSK_cons(ceil(M * rand(N, 1)));

    Perrors_PSK = zeros(1,25);
    Perrors_theo_DPSK = zeros(1,25);
    for j = 1:25
        SNR = j-5;
        [true_sym_psk, est_sym_psk, root_Eb] ...
            = simulate_transmission_diff(DPSK_cons, PSK_sym, N0, SNR);
        Perror = num_errors(true_sym_psk, est_sym_psk)/N;
        Perrors_PSK(j) = Perror;
        
        % see note about fudge factor
        Perrors_theo_DPSK(j) = 2*qfunc(root_Eb*sqrt(2/N0)*sin(pi/M));
    end

    scatter(-4:20, Perrors_DPSK, colors(i));
    plot(-4:20, Perrors_theo_DPSK, colors(i));
end

hold off;
title('DPSK bit error rate vs. SNR');
ylabel('P_{error}');
xlabel('SNR (dB)');
legend([
    "M=4 empirical", "M=4 theoretical", ...
    "M=8 empirical", "M=8 theoretical", ...
    "M=16 empirical", "M=16 theoretical", ...
    "M=32 empirical", "M=32 theoretical"
]);
exportgraphics(gcf(), 'dpsk.eps');

fprintf('Elapsed for DPSK for N=%d: %.2fs\n', N, toc());
tic();

%% QAM
Ms = [16, 32, 64];
colors = ["b", "r", "g", "k"];

figure('visible', 'off', 'position', [0 0 1000 500]);
tiledlayout(1, 1, 'TileSpacing', 'Compact');
nexttile();
hold on;
for i = 1:length(Ms)
    M = Ms(i);
    
    % generate constellation (assume M is a power of 2)
    if mod(log2(M), 2) == 0
        x = (1:sqrt(M)) - (sqrt(M)+1)/2;
        % use broadcasting to generate sqrt(M)*sqrt(M) square
        qam_cons = x + x.'*1j;
    else
        x = (1:sqrt(M/2)) - (sqrt(M/2)+1)/2;
        y = (1:sqrt(M*2)) - (sqrt(M*2)+1)/2;
        % use broadcasting to generate sqrt(M/2)*sqrt(M*2) rectangular
        qam_cons = x + y.'*1j;
    end
    qam_cons = qam_cons(:);     % flatten
    
    % generate signal sequence
    qam_sym = qam_cons(ceil(M * rand(N, 1)));

    Perrors_QAM = zeros(1, 25);
    Perrors_theo_QAM = zeros(1, 25);
    for j = 1:25
        SNR = j-5;
        [qam_true_sym, qam_est_sym, scale_factor] ...
            = simulate_transmission(qam_cons, qam_sym, N0, SNR);
        Perrors_QAM(j) = num_errors(qam_true_sym, qam_est_sym) / N;
        
        % see note about fudge factor
        E_avg = mean(abs(qam_true_sym).^2);
        if mod(log2(M), 2) == 0
            Prootm = 2*(1-1/sqrt(M))*qfunc(sqrt(2*3*E_avg/((M-1) * N0)));
            Perrors_theo_QAM(j) = 1 - (1 - Prootm)^2;
        else
            Perrors_theo_QAM(j) = 1 - (1 - 2*qfunc(sqrt(2*3*E_avg/((M-1)*N0))))^2;
        end
    end 
    
    scatter(-4:20, Perrors_QAM, colors(i));
    plot(-4:20, Perrors_theo_QAM, colors(i));
end

hold off;
title(sprintf('QAM bit error rate vs. SNR'));
ylabel('P_{error}');
xlabel('SNR (dB)');
legend([
    "M=16 empirical", "M=16 theoretical", ...
    "M=32 empirical", "M=32 theoretical", ...
    "M=64 empirical", "M=64 theoretical" ...
]);
exportgraphics(gcf(), 'qam.eps');

fprintf('Elapsed for QAM for N=%d: %.2fs\n', N, toc());