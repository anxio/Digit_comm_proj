% Test of DQPSK and detect function

nr_data_bits=10;
b = random_data(nr_data_bits);
r=dqpsk(b);
b=double(b);
bhat=ddetect(r);
n=norm(bhat-b);