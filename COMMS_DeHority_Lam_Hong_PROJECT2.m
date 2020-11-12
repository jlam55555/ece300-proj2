% number of samples
N = 10000;

% noise power
N0 = 3;

%% binary antipodal

% binary antipodal constellation and random signal sequence
bin_ap_const = [-1; 1];
bin_ap_sym = bin_ap_const(2 * ceil(rand(N, 1)));

% loop through SNRs from -4 to 20
Perrorsap_SNR = zeros(1, 25);
Theoretical_Perror_ap = zeros(1, 25);
for i = 1:25
    SNR = i-5;
    
    [true_sym_ap, est_sym_ap, root_Eb] = simulate_transmission(bin_ap_const, bin_ap_m, N0, SNR);
    Perror = num_errors(true_sym_ap, est_sym_ap) / N;
    Perrorsap_SNR(i) = Perror;
    Theoretical_Perror_ap(i) = qfunc(root_Eb*sqrt(4/N0));   % Pg. 406 of textbook WITH FUDGE FACTOR
    
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
for i = [1:25]
    SNR = i-5;
    [true_sym_ortho, est_sym_ortho, root_Eb] = simulate_transmission(bi_ortho_cons, bi_ortho_m, N0, SNR);
    Perror = num_errors(true_sym_ortho, est_sym_ortho)/N;
    Perrorsortho_SNR(i) = Perror;
    Theoretical_Perror_ortho(i) = qfunc(root_Eb*sqrt(2/N0));  %Pg. 408 of textbook WITH FUDGE FACTOR
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
M = 32
theta = linspace(0,2*pi,M+1);
theta = theta(2:M+1).';

fouraryPSK_cons = cos(theta)+j*sin(theta);

rand12 = floor(4*rand(N,1))+1;
fouraryPSK_m = fouraryPSK_cons(round(rand12));

PerrorsfouraryPSK_SNR = zeros(1,25);
Theoretical_Perror_fourary_PSK = zeros(1,25);
for i = [1:25]
    SNR = i-5;
    [true_sym_fourPSK, est_sym_fourPSK, root_Eb] = simulate_transmission(fouraryPSK_cons, fouraryPSK_m, N0, SNR);
    Perror = num_errors(true_sym_fourPSK, est_sym_fourPSK)/N;
    PerrorsfouraryPSK_SNR(i) = Perror;
    Theoretical_Perror_fouraryPSK(i) = 2*qfunc(root_Eb*sqrt(4/N0)*sin(pi/M));  %Pg. 416 of textbook WITH FUDGE FACTOR
end

figure;
scatter([-4:20],PerrorsfouraryPSK_SNR);
title('4-ary PSK')
hold on;
plot([-4:20], Theoretical_Perror_fouraryPSK);

%% DPSK -- needs to be fixed
% M = 32
% thetaD = linspace(0,2*pi,M+1);
% thetaD = thetaD(2:M+1).';
% 
% fouraryDPSK_cons = cos(thetaD)+j*sin(thetaD);
% 
% rand12D = floor(4*rand(N,1))+1;
% fouraryDPSK_m = fouraryDPSK_cons(round(rand12D));
% 
% PerrorsfouraryDPSK_SNR = zeros(1,25);
% Theoretical_Perror_fourary_DPSK = zeros(1,25);
% for i = [1:25]
%     SNR = i-5;
%     [true_sym_fourDPSK, est_sym_fourDPSK, root_Eb] = simulate_transmission(fouraryDPSK_cons, fouraryDPSK_m, N0, SNR);
%     Perror = num_errors_D(true_sym_fourDPSK, est_sym_fourDPSK)/N;
%     PerrorsfouraryDPSK_SNR(i) = Perror;
%     Theoretical_Perror_fouraryDPSK(i) = 2*qfunc(root_Eb*sqrt(4/N0)*sin(pi/M));  %Pg. 416 of textbook WITH FUDGE FACTOR
% end
% 
% figure;
% scatter([-4:20],PerrorsfouraryDPSK_SNR);
% title('4-ary DPSK')
% hold on;
% plot([-4:20], Theoretical_Perror_fouraryDPSK);

%% DPSK
M = 4
theta = linspace(0,2*pi,M+1);
theta = theta(2:M+1).';

fouraryPSK_cons = cos(theta)+j*sin(theta);

rand12 = floor(4*rand(N,1))+1;
fouraryPSK_m = fouraryPSK_cons(round(rand12));

PerrorsfouraryPSK_SNR = zeros(1,25);
Theoretical_Perror_fourary_PSK = zeros(1,25);
for i = [1:25]
    SNR = i-5;
    [true_sym_fourPSK, est_sym_fourPSK, root_Eb] = simulate_transmission_diff(fouraryPSK_cons, fouraryPSK_m, N0, SNR);
    Perror = num_errors(true_sym_fourPSK, est_sym_fourPSK)/N;
    PerrorsfouraryPSK_SNR(i) = Perror;
    Theoretical_Perror_fouraryPSK(i) = 2*qfunc(root_Eb*sqrt(4/N0)*sin(pi/M));  %Pg. 416 of textbook WITH FUDGE FACTOR
end

figure;
scatter([-4:20],PerrorsfouraryPSK_SNR);
title('4-ary PSK')
hold on;
plot([-4:20], Theoretical_Perror_fouraryPSK);

%% QAM
M = 16;
