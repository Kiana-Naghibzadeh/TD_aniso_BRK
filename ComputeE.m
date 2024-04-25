function e=ComputeE(ori1,ori2)

load('nf_Coarse.mat')
n=nf;
e = transpose(EulToEP(ori1,ori2, n(:,1),n(:,2),n(:,3)));

end

