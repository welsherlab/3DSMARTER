% find the projection of a point on a line

% P,Q: two points that define a line
% R: a point outside the line determined by P and Q

function out = pt_proj_on_line(P,Q,R)

P = P(:); Q = Q(:); R = R(:);

t =  sum((P - Q).*(R - P))./sum((P-Q).^2);

out = t.*(P-Q) + P;