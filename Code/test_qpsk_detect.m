% Test of QPSK and detect function

nr_data_bits=10;
b = random_data(nr_data_bits);

r=qpsk(b);
b=double(b);
bhat=detect(r);
n=norm(bhat-b);