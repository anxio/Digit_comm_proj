function d = qpsk(b)
% d = qpsk(b)
%
% Map the bits to be transmitted into QPSK symbols using Gray coding. The
% resulting QPSK symbol is complex-valued, where one of the two bits in each
% QPSK symbol affects the real part (I channel) of the symbol and the other
% bit the imaginary part (Q channel). Each part is subsequently PAM
% modulated to form the complex-valued QPSK symbol. The energy per QPSK
% symbol is normalized to unity.
%
% The mapping resulting from the two PAM branches are:
%
% complex part (Q channel)
%         ^
%         |
%  10 x   |   x 00   (odd bit, even bit)
%         |
%  -------+------->  real part (I channel)
%         |
%  11 x   |   x 01
%         |
%
%
%
% Input:
%   b = bits {0, 1} to be mapped into QPSK symbols
%
% Output:
%   d = complex-valued QPSK symbols

temp=zeros(1,length(b)/2); %2 bits represent 1 symbol
%Definition of the QPSK symbols using Gray Coding
for n=1:length(b)/2
    p=b(2*n);%Pair bits
    imp=b(2*n - 1);%Impair bits
    if (imp==0)&&(p==0)
        temp(n)=exp(1i*pi/4); %45 degres symbol
    end
    if (imp==1)&&(p==0)
        temp(n)=exp(1i*3*pi/4); %135 degres symbol
    end
    if (imp==1)&&(p==1)
        temp(n)=exp(1i*5*pi/4); %225 degres symbol
    end
    if (imp==0)&&(p==1)
        temp(n)=exp(1i*7*pi/4); %315 degres symbol
    end
end
d=temp;
% figure(1);
% plot(d,'o');%plot constellation without noise
% axis([-2 2 -2 2]);
% grid on;
% xlabel('real'); ylabel('imag');
% title('QPSK constellation');
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
