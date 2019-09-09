function y = interp(y1,y2,x1,x2,x)
%INTERP Summary of this function goes here
%   Detailed explanation goes here
y = y1 + (y2-y1)./(x2-x1).*(x-x1);
end

