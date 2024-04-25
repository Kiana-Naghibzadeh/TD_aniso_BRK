function cval = convolvegrain(ind,val,dims,common,commongroup,current_grain,e,n,u,IC)
maxm=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cval = convolvefamily(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolves the characteristic function of union of grains whose
% level set representation is stored in "ind" and "val" input vars
% with the kernel "KERNEL" (global variable, assumed present).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global Z
% global KERNEL Z; % These global variables *must* have been already
                 % allocated before calling this function! Both are
                 % of dimension "dims".
ZP = -ones(dims);
n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.

%% Convert level set data to volume of fluid representation:
[x,y,z] = ind2sub(dims,ind);
vf = ls2vf3D(int32(x),int32(y),int32(z),val,ZP,dims(1),dims(2),dims(3));
%% Carry out the convolution:
Ku = zeros(dims);
Ku(ind) = vf; % Characteristic function of the union of grains.
M=fftn(Ku);
cval=cell([length(common),1]);
MyDir = '/jet/home/snaghibz/ondemand/DMREF/Ni/Numerical/KernelFFT/min8e-1/';

for ij=1:length(common) 
	if common(ij)==current_grain
		cval{ij}=zeros(n1,n2,n3);
    else
	KERNEL = Kernel5D(n1, n2, n3, 5e-5,1e-1,n,IC,u,e(commongroup(ij),:));
         cval{ij} = real(ifftn(M .* KERNEL));        
	end
end
