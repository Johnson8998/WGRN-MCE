function ge = MCE(data, net)
%Input: 
% data: gene expression dataset, p*n where p is the number of genes and n 
%   is the number of samples.
% net: p*p adjacency matrix of gene network (usually PPI). The diagnal is 1 NOT 0!!!
%Output:
% ge: a length n vector, ge(i) is the MCE of sample i.

[p,n] = size(data);
net(logical(speye(size(net)))) = 1;
[r, c] = find(net);
nx = nnz(net);

ge = zeros(1, n);

for i = 1:n
    p0 = data(:, i); p0 = p0/sum(p0);
    B = sparse(r, c, p0(r), p, p);
    [out, err, nstep] = fun_iteration(net, B, p0);
    
    ge(i) = - dot(p0, log(out(1:p).*out((p+1):end)));
    
    ge(i) = ge(i) - sum(p0(p0>0).*log(p0(p0>0)));
    fprintf('Sample %d/%d: %d iterations with entropy %f, error %f.\n',...
        i, n, nstep, ge(i), err);
end
maxMCE = log(nx); %normalization
ge = ge/maxMCE;
end


%Iteration to solve nonlinear Lagrangian multipliers equations.

function [out, err, n_it] = fun_iteration(A, B, p0)
%Input:
% r,c the nonzeros elements of network, [r,c]=find(net)
% p: the number of genes
% p0: the steady state distribution
%
%Output:
% when |x-x0| < epsilon, return x
% n_it is the count of iteration steps.
%
%Algorithm:
% x = (lambda; theta): length 2p vector.
% lambda0: length p vector, Lagrangian of \sum_j p_ij = 1, i=1:p
% theta0: length p vector, theta0+1 is Lagrangian of \sum_i pi_i*p_ij = pi_j, j=1:
% The transition probability p_ij = exp(theta0_j+lambda0_i/pi_i)
% theta_i = exp(theta0_i), lambda_i = exp(lambda0_i/p_i);
% p_ij = theta_j*lambda_i.


epsilon = 10^-2; Max_n = 10^6;
stop = 0; n_it = 0;
lambda0 = p0;
theta0 = ones(size(p0));
while(~stop)
    n_it = n_it+1;
    
    lambda = 1./(A*theta0);
%     lambda = 1./(A*theta);
%     theta = p0./(B'*lambda0);
    theta = p0./(B'*lambda);
    
    err = norm([lambda-lambda0; theta-theta0], Inf);
    stop = ( (err < epsilon) | (n_it>Max_n) );
    lambda0 = lambda; theta0 = theta;
end
out = [lambda; theta];
end
