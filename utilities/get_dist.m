function d = get_dist(P1,P2);

if ~isvector(P1)||~isvector(P2)||length(P1)~=length(P2)||~ismember(length(P1),[2,3])
   error('format of input is not correct! must be vectors of the same length (length = 2 or 3)') 
end
P1 = P1(:);
P2 = P2(:);

d = sqrt(sum((P1-P2).^2));