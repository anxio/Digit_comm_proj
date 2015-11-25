function phihat = phase_estimation(r, b_train)
% phihat = phase_estimation(r, b_train)
%
% Phase estimator using the training sequence. The phase estimate is
% obtained by minimizing the norm of the difference between the known
% transmitted QPSK-modulated training sequence and the received training
% part. NB! There are other ways of estimating the phase, this is just
% one example.
%
% Input:
%   r       = received baseband signal
%   b_train = the training sequence bits
%
% Output:
%   phihat     = estimated phase

% Setting of constants
nr_training_bits=length(b_train);
Q=8;

%Length of the QPSK training sequence
L = nr_training_bits / 2;

%Generated symbols of the QPSK training sequence
c= qpsk(b_train);

k= 1:L;

phihat= sum(angle(r(Q.*k).*conj(c)))./L;



end

