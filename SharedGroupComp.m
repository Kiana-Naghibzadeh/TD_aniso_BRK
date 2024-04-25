function [shared,sharedgroup, LN, indd] = SharedGroupComp(presence,N, grains, dataori,id,eFull)
    shared = cell([N, 1]);
    indd = cell([N, 1]);
    LN = zeros(N,1);
    load('final0.mat');
    
    for k=1:N % Loop over grains.
        indd{k} = grains{k,1};
        shared{k}=presence(indd{k},1); %local id
        shared{k}=[shared{k}{:}];
        shared{k} = unique(shared{k});
        LN(k)=length(shared{k}); 
    end
    sharedgroup = shared;
    indexx = final(:,1:2);
    m=length(final);
    for k=1:N % Loop over grains.
            LN(k)=length(shared{k}); 
            for i=1:LN(k)
                if k~=shared{k}(i)
                    b = find(ismember(indexx,sort([id(k), id(shared{k}(i))]),'rows'));
                    if isempty(b)
                        eNew = ComputeE(dataori(id(k),:), dataori(id(shared{k}(i)),:));
			m=m+1;
			eFull(m,:)= eNew;
			final(m,1:2) = sort([id(k), id(shared{k}(i))]);
			sharedgroup{k}(i) = m;
			fprintf('Interface %d, %d, is new and it belongs to %d cluster. \n',k,shared{k}(i),m)
			indexx = final(:,1:2);
                    else
                        sharedgroup{k}(i) = b;
                    end
                else 
                    sharedgroup{k}(i) = 0;
                end
            end
    end
save('eFull.mat','eFull','-v7.3')
save('SharedGroup.mat','sharedgroup','-v7.3')
save('final0.mat','final','-v7.3')
end
