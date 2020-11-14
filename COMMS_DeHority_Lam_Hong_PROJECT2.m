%% TODO: make variable names more clear & concise

% number of samples
N = 1000;

% noise power
N0 = 1;

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
    Theoretical_Perror_ap(i) = qfunc(root_Eb*sqrt(2)*sqrt(2/N0));   % Pg. 406 of textbook WITH FUDGE FACTOR
    
end

figure;
hold on;
scatter(-4:20, Perrorsap_SNR);
plot(-4:20, Theoretical_Perror_ap);
hold off;
title('Binary Antipodal');
ylabel('P_{error}');
xlabel('SNR (dB)');

%% binary orthogonal
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

figure;
hold on;
scatter(-4:20, Perrorsortho_SNR);
plot(-4:20, Theoretical_Perror_ortho);
hold off;
title('Binary Orthogonal');
ylabel('P_{error}');
xlabel('SNR (dB)');

%% PSK

Ms = [4 8 16 32];
figure();
tiledlayout(1, length(Ms), 'TileSpacing', 'Compact');
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

    nexttile();
    hold on;
    scatter(-4:20, PerrorsfouraryPSK_SNR);
    plot(-4:20, Theoretical_Perror_PSK);
    hold off;
    title(sprintf('PSK; M=%d', M));
    ylabel('P_{error}');
    xlabel('SNR (dB)');
end

%% DPSK

Ms = [4 8 16 32];
figure();
tiledlayout(1, length(Ms), 'TileSpacing', 'Compact');
for i = 1:length(Ms)
    M = Ms(i);
    
    % generate constellation
    theta = linspace(0, 2*pi, M+1);
    theta = theta(1:M).';
    dpsk_cons = cos(theta) + 1j*sin(theta);

    % generate signal sequence
    dpsk_sig = dpsk_cons(ceil(M * rand(N, 1)));

    PerrorsfouraryPSK_SNR = zeros(1, 25);
    Theoretical_Perror_PSK = zeros(1, 25);
    test = zeros(1, 25);
    for j = 1:25
        SNR = j-5;
        [dpsk_true_sym, dpsk_est_sym, root_Eb] = simulate_transmission_diff(dpsk_cons, dpsk_sig, N0, SNR);
        Perror = num_errors(dpsk_true_sym, dpsk_est_sym)/N;
        PerrorsfouraryPSK_SNR(j) = Perror;
        Theoretical_Perror_PSK(j) = 2*qfunc(root_Eb*sqrt(2/N0)*sin(pi/M));  %Pg. 416 of textbook WITH FUDGE FACTOR
    end

    nexttile();
    hold on;
    scatter(-4:20, PerrorsfouraryPSK_SNR);
    plot(-4:20, Theoretical_Perror_PSK);
    hold off;
    title(sprintf('DPSK; M=%d', M));
    ylabel('P_{error}');
    xlabel('SNR (dB)');
end

%% QAM
Ms = [4, 16, 64];

figure();
tiledlayout(1, length(Ms), 'TileSpacing', 'Compact');
for i = 1:length(Ms)
    M = Ms(i);
    
    % generate constellation
    x = (1:sqrt(M)) - (sqrt(M)+1)/2;
    qam_cons = x + x.'*1j;      % use broadcasting to generate sqrt(M)*sqrt(M) square
    qam_cons = qam_cons(:);     % flatten
    
    % generate signal sequence
    qam_sym = qam_cons(ceil(M * rand(N, 1)));

    Perrors_QAM = zeros(1, 25);
    Perrors_theo_QAM = zeros(1, 25);
    for i = 1:25
        SNR = i-5;
        [qam_true_sym, qam_est_sym, scale_factor] = simulate_transmission(qam_cons, qam_sym, N0, SNR);
        Perrors_QAM(i) = num_errors(qam_true_sym, qam_est_sym) / N;
        
        E_avg = mean(abs(qam_true_sym).^2);
%         Perrors_theo_QAM(i) = 4*qfunc(sqrt(2*3*E_avg/((M-1)*N0)));
        Prootm = 2*(1-1/sqrt(M))*qfunc(sqrt(2*3*E_avg/((M-1) * N0)));
        Perrors_theo_QAM(i) = 1 - (1 - Prootm)^2;
    end
    
    nexttile();
    hold on;
    scatter(-4:20, Perrors_QAM);
    plot(-4:20, Perrors_theo_QAM);
    hold off;
    title(sprintf('QAM; M=%d', M));
    ylabel('P_{error}');
    xlabel('SNR (dB)');
end