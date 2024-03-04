function r = bafun(v, f)
%BAFUN  Binary Associative Function
%   r = BAFUN(v,f)
%   v is a vector with at least two elements
%   f is a handle to a function that operates on two numbers

if numel(v) > 2
    r = bafun({f(v{1}, v{2}) v{3:end}}, f);
else
    r = f(v{1}, v{2});
end
