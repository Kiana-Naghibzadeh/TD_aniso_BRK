function EdistHPC(varrr)
    MyDir = '/jet/home/snaghibz/ondemand/DMREF/Ni/Numerical/Code/gbm3dg35/5D/An0-1/energy/t1_Smooth.dream3d'
    Euler= h5read(MyDir,'/DataContainers/TriangleDataContainer/FaceData/EulerAngles');
    nn= h5read(MyDir,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals');
    areaaa= h5read(MyDir,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas');
    FeatID= h5read(MyDir,'/DataContainers/TriangleDataContainer/FaceData/FeatureIds');
    
    FeatID = double(FeatID');
    x = [double(Euler') nn' areaaa' FeatID]; 
    id=1:length(x);
    x=[x, id'];
    [r,~]= find(isnan(x));
    x(unique(r),:)=[];
    
    FeatID(unique(r),:) = [];
    %% Removing free boundary
    aa1 = find(FeatID(:,1) == 0);
    aa2 = find(FeatID(:,2) == 0);
    a = union(aa1,aa2);
    x(a,:) = [];
    %%
    n = x(:,7:9);
    n = num2cell(n,2);
    ori1 = x(:,1:3);
    ori1 = num2cell(ori1,2);
    ori2 = x(:,4:6);
    ori2 = num2cell(ori2,2);
    area = x(:,10);
    %%
    load('final0.mat')
    load('eFull.mat')
    load('nf_Coarse.mat')
    indexx = final(:,1:2);
    m=length(final);
    e=zeros(size(area));
    Loc=zeros(size(area));
    %%
    proc=10;
    alpha = (length(area)/proc);
    startouter = floor((varrr-1)*alpha) + 1;
    endouter = floor(varrr*alpha);
    
    for i=startouter:endouter
          b = find(ismember(indexx,sort([x(i,11), x(i,12)]),'rows'));
          if isempty(b)
                i
                eNew = ComputeE(x(i,1:3), x(i,4:6));
			    m=m+1;
			    eFull(m,:)= eNew;
			    final(m,1:2) = sort([x(i,11), x(i,12)]);
			    indexx = final(:,1:2);
			    b=m;
	      end
          a = repmat(x(i,7:9),length(nf),1);
          distt = a - nf;
          distE = sum(distt.^2,2); 
          h = find(distE == min(distE));
          e(i) = eFull(b(1),h(1));
          Loc(i) = b(1);
    end
    
    
    fname = sprintf('e%d.mat', varrr);
    MyFile1=strcat('/jet/home/snaghibz/ondemand/DMREF/Ni/Numerical/Code/gbm3dg35/5D/An0-1/energy/',fname);
    save(MyFile1,'e','-v7.3')
    
    fname = sprintf('Loc%d.mat', varrr);
    MyFile1=strcat('/jet/home/snaghibz/ondemand/DMREF/Ni/Numerical/Code/gbm3dg35/5D/An0-1/energy/',fname);
    save(MyFile1,'Loc','-v7.3')
    
    save('eFull.mat','eFull','-v7.3')
    save('final0.mat','final','-v7.3')

end
