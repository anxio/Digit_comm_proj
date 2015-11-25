% Test of QPSK and detect function

nr_data_bits=10;
b = random_data(nr_data_bits);
k=qpsk(b);
bhat=detect(k);