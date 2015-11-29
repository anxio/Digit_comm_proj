%%%
%%% Test of sync function
%%%

%% Using the data of simulation
load('data_test_sync_without_noise');
b_double=double(b_train);
max_hyp=b_double*b_double';
t_start=1+Q*nr_guard_bits/2;
t_end=t_start+length(b_train)/2;
t_samp = sync(mf, b_train, Q, t_start, t_end);

%% Homemade data
nr_bits=100;
Q=8;

b = training_sequence(nr_bits);
d = qpsk(b);

%pulse_shape = ones(1, Q);
pulse_shape = root_raised_cosine(Q);
mf_pulse_shape = fliplr(pulse_shape);
tx = upfirdn(d, pulse_shape, Q, 1);
rx=tx; t_samp=1;
mf=conv(mf_pulse_shape,rx);
r = mf(t_samp:Q:t_samp+Q*nr_bits/2-1);
bhat=detect(r);
norm_diff=norm(bhat-b);

t=1:length(b);
corr=(repmat(bhat,length(b),1)+repmat(t',1,length(b))) * b';
corr2=xcorr(b,bhat);
figure
    hold all
    plot(corr2);
    plot(xcorr(b));
    hold off;
