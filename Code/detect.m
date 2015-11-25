function bhat = detect(r)
% bhat = detect(r)
%
% Computes the received bits given a received sequence of (phase-corrected)
% QPSK symbols. Gray coding of the individual bits is assumed. Hence, the
% two bits for each symbol can be detected from the real and imaginary
% parts, respectively. The first of the two bits below is output first in
% the bhat-sequence.
%
% Assumed mapping:
%
%  10 x   |   x 00
%         |
%  -------+-------
%         |
%  11 x   |   x 01
%
% Input:
%   r  = sequence of complex-valued QPSK symbols
%
% Output:
%   bhat  = bits {0,1} corresponding to the QPSK symbols

%Computing of the phase of each symbol
r_phase=angle(r);

%Computing of the first bit of each symbol
r_phase_0=r_phase;
r_phase_0((abs(r_phase_0)<=pi/2))=0;
r_phase_0((abs(r_phase_0)>pi/2))=1;

%Computing of the second bit of each symbol
r_phase_1=r_phase;
r_phase_1((r_phase_1>=0))=0;
r_phase_1((r_phase_1<0))=1;

%Merging of the first and second bit for each symbol
inter=[r_phase_0 ; r_phase_1];
bhat=reshape(inter,2*size(r,2),1);

bhat=bhat.';
end


