function bhat = ddetect(r)
% bhat = detect(r)
%
% DQPSK detection
% Input:
%   r  = sequence of complex-valued QPSK symbols
%
% Output:
%   bhat  = bits {0,1} corresponding to the QPSK symbols

%Computing of the phase of each symbol
r_phase=angle(r)-[pi/4,angle(r(1:end-1))];
r_phase=mod(r_phase,2*pi);

%Computing of the first bit of each symbol
r_phase_0=r_phase;
r_phase_0(((abs(r_phase_0)<=pi/4)|(abs(r_phase_0)>5*pi/4)))=0;
r_phase_0(~((abs(r_phase_0)<=pi/4)|(abs(r_phase_0)>5*pi/4)))=1;

%Computing of the second bit of each symbol
r_phase_1=r_phase;
r_phase_1(((abs(r_phase_1)<3*pi/4)|(abs(r_phase_1)>=7*pi/4)))=0;
r_phase_1(~((abs(r_phase_1)<3*pi/4)|(abs(r_phase_1)>=7*pi/4)))=1;

%Merging of the first and second bit for each symbol
inter=[r_phase_0 ; r_phase_1];
bhat=reshape(inter,2*size(r,2),1);

bhat=bhat.';
end


