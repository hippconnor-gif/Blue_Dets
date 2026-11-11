function y = hat(x)
% HAT triangular window function used for kernels
y = max(0, 1 - abs(x));
end
