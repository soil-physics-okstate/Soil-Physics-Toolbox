function [y] = cumnansum(x)
%CUMNANSUM returns the cumulative sums of the columns in array "x" ignoring
%nan's

i = isnan(x); %find nan's
x(i) = 0; %replace with zeros
y = cumsum(x); %calculate cumulative sum