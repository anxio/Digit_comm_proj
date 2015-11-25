function t_samp = sync(mf, b_train, Q, t_start, t_end)
% t_samp = sync(mf, b_train, Q, t_start, t_end)
%
% Determines when to sample the matched filter outputs. The synchronization
% algorithm is based on cross-correlating a replica of the (known)
% transmitted training sequence with the output from the matched filter
% (before sampling). Different starting points in the matched filter output
% are tried and the shift yielding the largest value of the absolute value
% of the cross-correlation is chosen as the sampling instant for the first
% symbol.
%
% Input:
%   mf            = matched filter output, complex baseband, I+jQ, (r)
%   b_train       = the training sequence bits, known ()
%   Q             = number of samples per symbol, known (Q)
%   t_start       = start of search window, known (tstart)
%   t_end         = end of search window, known (tend)
%
% Output:
%   t_samp = sampling instant for first symbol

% Setting of constants
nr_training_bits=length(b_train);

%Length of the QPSK training sequence
L = nr_training_bits / 2;

%Generated symbols of the QPSK training sequence
c= qpsk(b_train);

%Vectors for the time and for the summation
%t = linspace(t_start,t_end,length(b_train/2)); CORR ADRI
t= t_start:t_end;

k= 1:L;

%Creation of the timeshift index matrix (i,j) of size (length(t),length(k))
% Q.k(j)+t(i)
mf_index=repmat(k.*Q,length(t),1) + repmat(t.',1,length(k));

%Computing of mf, the ith raw corresponds to the received sequence
% with a timeshift of  t(i) 
mf_matrix= mf(mf_index);

%Computing of the cross-correlation with the training sequence
corr = abs(mf_matrix * c.');

%Computing of t_samp
[~,i_samp] = max(corr);
t_samp=t(i_samp);

%Displaying of the cross-correlation 
figure(1)
    hold all;
    plot(t,corr,'b'); 
    xlabel('Time shift');
    ylabel('Norm of the correlation function');
end








