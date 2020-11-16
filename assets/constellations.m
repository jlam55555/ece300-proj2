%% visualizing the base constellations

figure('visible', 'off', 'position', [0 0 1500 1500]);
tiledlayout(3, 3, 'TileSpacing', 'Compact');

% binary antipodal
pc([-1; 1], 'Binary antipodal');

% binary orthogonal
pc([1; 1j], 'Binary orthogonal');

% (D)PSK
for M = [4 8 16 32]
    pc(exp(1j*(0:(2*pi/M):(2*pi-0.1)).'), sprintf('(D)PSK M=%d', M))
end

% QAM
for M = [16 32 64]
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
    pc(qam_cons(:), sprintf('QAM M=%d', M));
end

exportgraphics(gcf(), 'constellations.eps');

% plot constellation
function pc(const, name)
    nexttile();
    scatter(real(const), imag(const), 'x');
    title(name);
    ylabel('Q');
    xlabel('I');
    axis equal;
    grid on;
end