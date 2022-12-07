function [rows,row_ct] = ct_idx3(idx3,dims) % dims: size of each dimension like [3,4,5]
    if ~isvector(dims)
        error('invalid input!')
    end
    if length(dims)~=3
        error('invalid input!')
    end
    
    v3 = reshape(1:prod(dims) ,dims);
    idx1 = zeros(length(idx3(:,1)),1);
    for i = 1:length(idx3(:,1))
    %idx1(i) = v3(idx3( i,:) );
    idx1(i) = v3(idx3( i,1) ,idx3( i,2),idx3( i,3));
    end
    
    N = histcounts(idx1,.5:1:prod(dims)+.5);
    N_idx = find(N~=0);
    rows = zeros(length(N_idx),3);
    row_ct = zeros(length(N_idx),1);
    
    for i = 1:length(N_idx)
        row_ct(i) = N(N_idx(i));
        [a,b,c] = find(v3 == N_idx(i));
        rows(i,:) = [a,b,c];
    end
    
end