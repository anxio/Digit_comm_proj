function pulse_shape = root_raised_cosine(Q, b, trunc)
if nargin < 3
if nargin < 2
% The symbol time T is assumed to be unity.
T = 1*Q;
pulse_shape = 4 * b * ...




