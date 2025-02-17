
function [ge, p_max] = MCE_gpu(data_matrix, net_matrix)
%Input: 
% data: gene expression dataset, p*n where p is the number of genes and n 
%   is the number of samples.
% net: p*p adjacency matrix of gene network (usually PPI). The diagnal is 1 NOT 0!!!
%Output:
% ge: a length n vector, ge(i) is the MCE of sample i.
% % % data = load('rnatest2.csv');
% % % net = load('adj_matrix.csv');

data_gpu = gpuArray(data_matrix);
[p,n] = size(data_gpu);

% net_gpu = net_gpu + speye(size(net_gpu));
net_matrix(logical(speye(size(net_matrix)))) = 1;
% net_gpu = sparse(net_gpu);%19892*19892
net_gpu = gpuArray(net_matrix);
[r, c] = find(net_gpu);
nx = nnz(net_gpu);

ge_gpu = gpuArray(zeros(1, n));
% p_max_gpu = cell(2,n);
p_max = cell(2,n);

for i = 1:n
    p0 = data_gpu(:, i); p0 = p0/sum(p0);
    % p0 = full(p0);
    B = sparse(r, c, p0(r), p, p);
    % tic
    [out_gpu, err_gpu, nstep] = Code_5_1_fun_iteration_gpu(net_gpu, B, p0);
    % toc
    
    ge_gpu(i) = - dot(p0, log(out_gpu(1:p).*out_gpu((p+1):end)));
    %ge_gpu(i) = - dot(p0,log(complex(out_gpu(1:p)).*complex(out_gpu((p+1):end))));%有时候会提醒gpu计算需要变成complex
    ge_gpu(i) = ge_gpu(i) - sum(p0(p0>0).*log(p0(p0>0)));

    % p_max_gpu = out_gpu(1:p) * out_gpu((p+1):end)' .* net_gpu; %Out of memory on device.

    % f = fopen('p_max.xls','w');
    % fprintf(f,'%s','source');fprintf(f,'\t');
    % fprintf(f,'%s','target');fprintf(f,'\t');
    % fprintf(f,'%s','weight');fprintf(f,'\r\n');
    % for j = 1:nx
    %     fprintf(f,'%d',r(j));fprintf(f,'\t');
    %     fprintf(f,'%d',c(j));fprintf(f,'\t');
    %     fprintf(f,'%f',out_gpu(r(j))*out_gpu(p+c(j)));fprintf(f,'\r\n');
    % end
    % fclose(f)

    % p_max_gpu{1,i} = out_gpu(1:p);
    % p_max_gpu{2,i} = out_gpu((p+1):end);

    p_max(:,i) = {gather(out_gpu(1:p));gather(out_gpu((p+1):end))};
    % p_max{1,i} = gather(out_gpu(1:p));
    % p_max{2,i} = gather(out_gpu((p+1):end));
    
    fprintf('Sample %d/%d: %d iterations with entropy %f, error %f.\n',...
        i, n, nstep, ge_gpu(i), err_gpu);
end
maxMCE = log(nx); %normalization
ge_gpu = ge_gpu/maxMCE;
ge = gather(ge_gpu);
end


