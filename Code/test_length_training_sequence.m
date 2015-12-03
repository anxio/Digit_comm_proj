% Skeleton code for simulation chain

% History:
%   2000-06-28  written /Stefan Parkvall
%   2001-10-22  modified /George Jongren


% Initialization
EbN0_db = 0:10:20;                     % Eb/N0 values to simulate (in dB)
nr_bits_per_symbol = 2;             % Corresponds to k in the report
nr_guard_bits = 50;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 1000;                % Size of each data sequence (in nr bits)
nr_training_bits = 10:2:30;             % Size of training sequence (in nr bits)
    nr_blocks = 400;                     % The number of blocks to simulate

Q = 8;                              % Number of samples per symbol in baseband

% Define the pulse-shape used in the transmitter. 
% Pick one of the pulse shapes below or experiemnt
% with a pulse of your own.
pulse_shape = ones(1, Q);
%pulse_shape = root_raised_cosine(Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);


% Loop over different values of Eb/No.
nr_errors = zeros(length(EbN0_db), length(nr_training_bits));   % Error counter
phi_estimation = zeros(length(EbN0_db),length(nr_training_bits)); 
t_samp_estimation=zeros(length(EbN0_db),length(nr_training_bits)); 


for training_point = 1:length(nr_training_bits)
%Loop over different length of training sequence
for snr_point = 1:length(EbN0_db)
  % Loop over several blocks to get sufficient statistics.
  for blk = 1:nr_blocks

    %%%
    %%% Transmitter
    %%%

    % Generate training sequence.
    b_train = training_sequence(nr_training_bits(training_point));
    
    % Generate random source data {0, 1}.
    b_data = random_data(nr_data_bits);

    % Generate guard sequence.
    b_guard = random_data(nr_guard_bits);
 
    % Multiplex training and data into one sequence.
    b = [b_guard b_train b_data b_guard];
    
    % Map bits into complex-valued QPSK symbols.
    d = qpsk(b);

    % Upsample the signal, apply pulse shaping.
    tx = upfirdn(d, pulse_shape, Q, 1);
    theta=pi/2;
    tx=tx*exp(1i*theta);

    timeshift=4;
    tx=[ random_data(timeshift-1),tx(1:length(tx)-(timeshift-1))];
    
    %%%
    %%% AWGN Channel
    %%%
    
    % Compute variance of complex noise according to report.
    sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol / 10.^(EbN0_db(snr_point)/10);

    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(tx))+j*randn(size(tx)));
    %n=zeros(size(tx));    

    % Received signal
    rx = tx + n;
    
    %%%
    %%% Receiver
    %%%
    
    % Matched filtering.
    mf=conv(mf_pulse_shape,rx);
    
    % Synchronization. The position and size of the search window
    % is here set arbitrarily. Note that you might need to change these
    % parameters. Use sensible values (hint: plot the correlation
    % function used for syncing)! 
    t_start=1+Q*nr_guard_bits/2;
    t_end=t_start+nr_training_bits(training_point)/2;
    t_samp = sync(mf, b_train, Q, t_start, t_end);
    t_samp_estimation(snr_point,training_point)=t_samp_estimation(snr_point,training_point)+t_samp;
    

    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r = mf(t_samp:Q:t_samp+Q*(nr_training_bits(training_point)+nr_data_bits)/2-1);

    % Phase estimation and correction.
    %phihat = phase_estimation(r, b_train);
    phihat=phase_estimation(r,b_train);
    phi_estimation(snr_point,training_point)=phi_estimation(snr_point,training_point)+phihat;
    
    r=r * exp(-j*phihat);
    % Make decisions. Note that dhat will include training sequence bits
    % as well.
    bhat = detect(r);
    % Count errors. Note that only the data bits and not the training bits
    % are included in the comparison. The last data bits are missing as well
    % since the whole impulse response due to the last symbol is not
    % included in the simulation program above.
    temp=bhat(1+nr_training_bits(training_point):nr_training_bits(training_point)+nr_data_bits) ~= b_data;
    nr_errors(snr_point,training_point) = nr_errors(snr_point,training_point) + sum(temp);
    % Next block.
  end
  % Next Eb/No value.
end
% next training sequence
end

BER=nr_errors./nr_blocks./nr_data_bits;

phi_estimation=phi_estimation./nr_blocks;
t_samp_estimation=t_samp_estimation./nr_blocks;

%% Displaying
figure
hold all
for snr_point=1:length(EbN0_db)
    p=plot(nr_training_bits,BER(snr_point,:));
end
legend(strcat(num2str(EbN0_db(1)),32,'dB'), strcat(num2str(EbN0_db(2)),32,'dB') , strcat(num2str(EbN0_db(3)),32,'dB'))
set(gca, 'FontSize', 12);
ylabel('BER')
xlabel('Size of the sequence ');
title(strcat('Study of the synchronization in function of the length of the training sequence with \tau =',32,num2str(timeshift-1),32,'and \theta =',32,num2str(theta)));

