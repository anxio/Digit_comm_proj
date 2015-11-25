% Test of sync function

load('data_test_sync_without_noise');

t_start=1+Q*nr_guard_bits/2;
t_end=t_start+70;
t_samp = sync(mf, b_train, Q, t_start, t_end);