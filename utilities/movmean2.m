function out = movmean2(in)

if ~isvector(in)
    error('input must be a vector!')
    
end
temp = movmean(in,2);
out = temp(2:end);


end