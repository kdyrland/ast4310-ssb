function [ u ] = partfunc_Na(temp)

u = zeros([1 3]);
theta = 5040 ./ temp;
c0 = 0.30955;
c1 = -0.17778;
c2 = 1.10594;
c3 = -2.42847;
c4 = 1.70721;
	
logU1 = ( c0 + c1 .* log10(theta) + c2 .* log10(theta).^2 + c3 .* log10(theta).^3 ...
	+ c4 .* log10(theta).^4 );
u(1) = 10.^logU1;
u(2) = 1;
u(3) = 6;

end
