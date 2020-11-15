%% TODO: make variable names more clear & concise

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
Perrorsap_SNR = zeros(1, 25);
Theoretical_Perror_ap = zeros(1, 25);
for i = 1:25
    SNR = i-5;
    
    [true_sym_ap, est_sym_ap, root_Eb] = simulate_transmission(bin_ap_const, bin_ap_sym, N0, SNR);
    Perror = num_errors(true_sym_ap, est_sym_ap) / N;
    Perrorsap_SNR(i) = Perror;
    Theoretical_Perror_ap(i) = qfunc(root_Eb*sqrt(4/N0));   % Pg. 406 of textbook WITH FUDGE FACTOR
    
end

fprintf('Elapsed for binary antipodal for N=%d: %.2fs\n', N, toc());
tic();

figure('visible', 'off', 'position', [0 0 1000 500]);
tiledlayout(1, 1, 'TileSpacing', 'Compact');
nexttile();
hold on;
scatter(-4:20, Perrorsap_SNR, 'b');
plot(-4:20, Theoretical_Perror_ap, 'b');

% binary orthogonal
bi_ortho_cons = [1; 1j];
bi_ortho_m = bi_ortho_cons(2 * ceil(rand(N, 1)));

Perrorsortho_SNR = zeros(1,25);
Theoretical_Perror_ortho = zeros(1,25);
for i = 1:25
    SNR = i-5;
    [true_sym_ortho, est_sym_ortho, root_Eb] = simulate_transmission(bi_ortho_cons, bi_ortho_m, N0, SNR);
    Perror = num_errors(true_sym_ortho, est_sym_ortho)/N;
    Perrorsortho_SNR(i) = Perror;
    Theoretical_Perror_ortho(i) = qfunc(root_Eb*sqrt(2)/sqrt(N0));  %Pg. 408 of textbook WITH FUDGE FACTOR
end

scatter(-4:20, Perrorsortho_SNR, 'r');
plot(-4:20, Theoretical_Perror_ortho, 'r');
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
    fouraryPSK_cons = cos(theta) + 1j*sin(theta);

    % generate signal sequence
    fouraryPSK_m = fouraryPSK_cons(ceil(M * rand(N, 1)));

    PerrorsfouraryPSK_SNR = zeros(1,25);
    Theoretical_Perror_PSK = zeros(1,25);
    for j = 1:25
        SNR = j-5;
        [true_sym_fourPSK, est_sym_fourPSK, root_Eb] = simulate_transmission(fouraryPSK_cons, fouraryPSK_m, N0, SNR);
        Perror = num_errors(true_sym_fourPSK, est_sym_fourPSK)/N;
        PerrorsfouraryPSK_SNR(j) = Perror;
        Theoretical_Perror_PSK(j) = 2*qfunc(root_Eb*sqrt(4/N0)*sin(pi/M));  %Pg. 416 of textbook WITH FUDGE FACTOR
    end

    scatter(-4:20, PerrorsfouraryPSK_SNR, colors(i));
    plot(-4:20, Theoretical_Perror_PSK, colors(i));
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
    fouraryPSK_cons = cos(theta) + 1j*sin(theta);

    % generate signal sequence
    fouraryPSK_m = fouraryPSK_cons(ceil(M * rand(N, 1)));

    PerrorsfouraryPSK_SNR = zeros(1,25);
    Theoretical_Perror_PSK = zeros(1,25);
    for j = 1:25
        SNR = j-5;
        [true_sym_fourPSK, est_sym_fourPSK, root_Eb] = simulate_transmission_diff(fouraryPSK_cons, fouraryPSK_m, N0, SNR);
        Perror = num_errors(true_sym_fourPSK, est_sym_fourPSK)/N;
        PerrorsfouraryPSK_SNR(j) = Perror;
        Theoretical_Perror_PSK(j) = 2*qfunc(root_Eb*sqrt(2/N0)*sin(pi/M));  %Pg. 416 of textbook WITH FUDGE FACTOR
    end

    scatter(-4:20, PerrorsfouraryPSK_SNR, colors(i));
    plot(-4:20, Theoretical_Perror_PSK, colors(i));
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
        [qam_true_sym, qam_est_sym, scale_factor] = simulate_transmission(qam_cons, qam_sym, N0, SNR);
        Perrors_QAM(j) = num_errors(qam_true_sym, qam_est_sym) / N;
        
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
tic();