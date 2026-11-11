function [s1, prev, sd] = expdecay(s, fRate, decay1, prev)

if isempty(prev)
    prev = zeros(size(s));
end

alpha = exp(-decay1 / fRate);

s1 = alpha .* prev + (1 - alpha) .* s;
sd = s1;
prev = s1;

end
