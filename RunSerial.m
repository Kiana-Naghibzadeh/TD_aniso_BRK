proc=10;
load('An0exp.mat')
load('eFull.mat')
parpool('local',4)
load('IC0.mat')

for t=0:34
	t
    N = size(grains,1);
    presence = get_nhd_grains(grains,dims(1)*dims(2)*dims(3));
    dataori2 = dataori(2:end,:);
    [shared,sharedgroup, LN, indd] = SharedGroupComp(presence,N,grains, dataori2,id,eFull);    
    load('eFull.mat')
    gbm3dParallel(grains,shared,sharedgroup, LN, indd, dims,id,eFull,IC);
    MyAfterGcvalPSCSerial(t,grains, shared, dataori, id, rmlist, dims); 
    
    clearvars -except t eFull IC
    fname = sprintf('Graint%d.mat', t+1);
    load(fname)
end

