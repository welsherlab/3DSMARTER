function out = ezgrid(in,step)

if ~isvector(in) || ~isnumeric(in)
   error('input 1 must be a vector!') 
end
if ~isscalar(step) || ~isnumeric(step)
    error('input 2 must be a scarlar!')
end

out = floor(min(in)./step).*step:step:ceil(max(in)./step).*step;


end