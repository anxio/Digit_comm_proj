% Test of sync function

load('data_test_sync_without_noise');
b_double=double(b_train);
max_hyp=b_double*b_double.';
t_start=1+Q*nr_guard_bits/2;
t_end=t_start+length(b_train)/2;
t_samp = sync(mf, b_train, Q, t_start, t_end);