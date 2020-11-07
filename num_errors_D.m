% num_errors_D: count the number of errors in a sequence of transmitted
% symbols assuming a differential scheme
%
% params:
% true_sym  = sequence of true symbols
% est_sym   = sequence of transmitted symbols
function errorsD = num_errors_D(true_sym, est_sym)
    threshold = 0.0001;
    est1_phase = est_sym(1)/abs(est_sym(1));
    true1_phase = true_sym(1)/abs(true_sym(1));
  
    est_sym_rot = est_sym.*est1_phase';
    true_sym_rot = true_sym.*true1_phase';
    
    errorsD = nnz(abs(true_sym_rot(2:end)-est_sym_rot(2:end))>threshold);
end