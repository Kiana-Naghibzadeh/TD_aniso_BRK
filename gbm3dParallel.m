function gbm3dParallelSerial( grains,shared,sharedgroup, LN, indd,dims,id,e,IC)

    dt=5e-4;
    N = size(grains,1); % Number of grains.
    MyDir3 = '/jet/home/snaghibz/ondemand/DMREF/Ni/Numerical/Data/min1e-1/dt1e-3/';
    maxm =  125.5;


    %% Auxiliary vars.
    % global WORKSPACE Z LONGARRAY;

    n1 = dims(1); % Size of computational grid.
    n2 = dims(2); % Size of computational grid.
    n3 = dims(3); % Size of computational grid.

    WORKSPACE = int32(-ones(dims)); % Needed for dilation.c.
    LONGARRAY = int32(zeros(120000000,1)); % Needed for dilation.c.
    Z = -ones(dims); % Another Workspace var.
                    % Needed in ls2vf3D.c.
    chi1 = (-n1/2:1:n1/2-1)*n1/n1*sqrt(dt);
    chi2 = (-n2/2:1:n2/2-1)*n1/n2*sqrt(dt);
    chi3 = (-n3/2:1:n3/2-1)*n1/n3/4*2.81*sqrt(dt);
    [chi1,chi2,chi3]=meshgrid(chi1,chi2,chi3);
    chi1=permute(chi1,[2 1 3]);
    chi2=permute(chi2,[2 1 3]);
    chi3=permute(chi3,[2 1 3]);
    uK5D = sqrt(chi1.^2 + chi2.^2 + chi3.^2);
    n11=chi1(:)./uK5D(:);
    n12=chi2(:)./uK5D(:);
    n13=chi3(:)./uK5D(:);
    nK5D = [n11, n12, n13];
    %% CONVOLUTION STEP
    for k=1:N
	    fname = sprintf('gcval%d.mat', id(k));
        MyFile1=strcat(MyDir3,fname);

      	 a1 = grains(k,:);
       	 ii=indd{k};
       	 sharedgrouppar = shared{k};
         h = convolvegrainPre(a1{1},a1{2},dims, sharedgrouppar,sharedgroup{k}, k,e,nK5D,uK5D,IC);
       	 b = single(zeros(length(ii),LN(k)));
       	 for i=1:LN(k)
                b(:,i) = h{i}(ii);
         end
    save(MyFile1, 'b','-v6')

    end

end
