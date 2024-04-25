function MyAfterGcvalPSCSerial(t,grains, shared, dataori, id, rmlist, dims)
% Average grain size goal: 1.8155e+03
%%
MyDir2 = '/jet/home/snaghibz/ondemand/DMREF/Ni/Numerical/Data/min1e-1/dt1e-3/';
N=length(grains);

for jj=1:N
    j=id(jj);
    fname = sprintf('gcval%d.mat', j);
    MyFile=strcat(MyDir2,fname);
    load(MyFile);
    grains{jj,3} = double(gcval);
end
%err

%%
ori = 2*pi*rand(N,1);
grain0=cell(1,3);
grain0{1,1}=rmlist;
grain0{1,2}=ones(size(rmlist,1),1);
grain0{1,3}=zeros(size(rmlist,1),1);
tempgrains=[grain0;grains];
tempid=[1;id+1];
tempori=[0;ori];
presence = get_nhd_grains(tempgrains,dims(1)*dims(2)*dims(3));
%%
psi = cell(N,1);
for i=1:N
    psi{i} = zeros(size(grains{i,1}));
    for j=shared{i}
        if j~=i
            m = find(shared{j}==i);
            b = zeros(dims);
            b(grains{j,1}) = grains{j,3}(:,m);
            psi{i} = psi{i} + b(grains{i,1});
        end
    end
end
%%
psi2 = cell(N+1,1);
psi2{1} = zeros(length(rmlist),1);
psi2(2:end)=psi;
%%
for i=1:length(presence)
    a = zeros(length(presence{i,1}),1);
    for j=1:length(presence{i,1})
        p1 = presence{i,1};
        p3 = presence{i,3};
        a(j)=psi2{p1(j)}(p3(j));
    end
    lm = find(a==min(a));
    for j=1:length(presence{i,1})
        p1 = presence{i,1};
        p3 = presence{i,3};
        
        if j==lm
            tempgrains{p1(j),2}(p3(j)) = 1;
        else
            tempgrains{p1(j),2}(p3(j)) = -1;
        end
    end
end
%%
grains=tempgrains(2:end,:);
[grains,id] = removeemptygrains(grains,dims,id);

%%
N= length(grains);
for i=1:N
    grains{i,3} = zeros( length(grains{i,1}) ,1 );
end
%%
WORKSPACE = int32(-ones(dims)); % Needed for dilation.c.
LONGARRAY = int32(zeros(120000000,1)); % Needed for dilation.c.
Z = -ones(dims); % Another Workspace var.
                % Needed in ls2vf3D.c.

  for k=1:N % Loop over grains.
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    cval = grains{k,3}; % Convolution vals. at those pixels.
    Z(ind) = val;      % Lev. set. representation on grid.
    posind = ind(val>0); % Pixels in the interior of grain.
    [x,y,z] = ind2sub(dims,posind);
    [x3,y3,z3] = dilation_fixedbd(int32(x),int32(y),int32(z),5,WORKSPACE,LONGARRAY); % Dilation.
    ind3 = sub2ind(dims,x3,y3,z3);   
    ind2 = rmdilation1(ind3,rmlist );
    val2 = Z(ind2);    % Level set vals.
    Z(ind2) = -1;      % Reset Z back to all -1's.
    Z(ind) = cval - 1; % Convolution values - 1.
    cval2 = Z(ind2);   % Convolution vals - 1.
    Z(ind2) = -1;      % Reset Z back to all -1's.
    grains{k,1} = ind2;   % Refresh grain's data structure.
    grains{k,2} = val2;   % Ditto.
    grains{k,3} = cval2 + 1; % Ditto.
  end % (for k). Loop over grains ends.
%%
clearvars -except dims dataori grains id rmlist t
fname = sprintf('Graint%d.mat', t+1);
save(fname)
%%
end
