% Fundamental matrix fitting code adapted from Peter Kovesi's
% implementation of 8-point fundamental matrix estimation. Need 8 or more
% keypoint correspondences.

function F = fundamentalA_fit(X)
npts = size(X,2);
if npts < 4
    error("not enough points");
end
x1 = X(1:3,:);
x2 = X(4:6,:);

% Build the constraint matrix
A = [x2(1,:)'   x2(2,:)'  x1(1,:)'  x1(2,:)' ones(npts,1) ];
[U,D,V] = svd(A,0); % Under MATLAB use the economy decomposition     
  
% Extract fundamental matrix from the column of V corresponding to
% smallest singular value.
params = V(:,5);
F = [0 0 params(1)
     0 0 params(2)
     params(3) params(4) params(5)];
        
% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,D,V] = svd(F,0);
F = U*diag([D(1,1) D(2,2) 0])*V';
    
end

