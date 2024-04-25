function [K] = Kernel5D(n1, n2, n3, dt,eps,n,IC,u,eCoarse)  
    es = eCoarse;
    beta=0;
    es2 = es - (beta * es);
    e = es2(IC);
    e = reshape(e,[n1,n2,n3]);
	
    K = exp(-u.^2 .* e.^2);
    K = fftshift(K);
end
