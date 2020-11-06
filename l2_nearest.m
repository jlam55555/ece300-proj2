% l2_nearest: Finds nearest symbols in a constellation to a set of received
% (noisy) signals
%
% params:
% con   = constellation (M x 1)
% est   = received signals (N x 1)
%
% returns:
% N_hat = estimated signals
function nearest = l2_nearest(con,  est)
    
    % transform complex numbers into a 2D real vectors
    con_2d = [real(con); imag(con)]';
    est_2d = [real(est); imag(est)]';
    
    % find nearest points (query points are the estimated vectors)
    nearest = con(dsearchn(con_2d, est_2d));
end