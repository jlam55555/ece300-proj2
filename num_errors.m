% num_errors: count the number of errors in a sequence of transmitted
% symbols assuming a nondifferential scheme
%
% params:
% true_sym  = sequence of true symbols
% est_sym   = sequence of transmitted symbols
%
function errors = num_errors(true_sym, est_sym)
    threshold = 0.0001;
    errors = nnz(abs(true_sym-est_sym)>threshold);
end