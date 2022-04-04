function [weight_v] = cos_kernels_weight(Kernels_list,adjmat,dim)
%È¨ÖØ¼ÆËã
%tju cs, bioinformatics. 
% adjmat : binary adjacency matrix
% dim    : dimension (1 - rows, 2 - cols)




num_kernels = size(Kernels_list,3);

v = randn(num_kernels,1);

y = adjmat;
    % Graph based kernel
if dim == 1
        ga = y*y';
else
        ga = y'*y;
end
K_values=[];
ga = Knormalized(ga);
value_ga = get_values(ga);
corr_list=[];
for i=1:num_kernels
	%S=Knormalized(Kernels_list(:,:,i));
	%Kernels_list(:,:,i) = S;
	mean_v = get_values(Kernels_list(:,:,i));
	K_values=[K_values;mean_v];
end

TEMP_V1 = norm(value_ga,2);
TEMP_V2 = K_values*value_ga';

cvx_begin
    variable v(num_kernels,1);
    minimize( norm(v'*K_values,2)*TEMP_V1 - v'*TEMP_V2);
	%subject to
	v >= 0;
	sum(v)==1;
cvx_end

weight_v = v;


end


function cell_v = get_values(K1)
	   
	size_k = size(K1,1);
	cell_v = 0;
	for i=1:size_k
		cell_v = [cell_v,K1(i,i+1:end)];
	end


end