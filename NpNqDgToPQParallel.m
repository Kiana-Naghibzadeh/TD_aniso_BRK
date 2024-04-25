function [e] = NpNqDgToPQParallel(dg,np,nq)
%np: normal in the frame of grain p, input is a row 
%nq: normal in frame of grain q, input is a row
%dg: map from grain p to q
    e = zeros(size(np,1),1);

    Pnew=zeros(size(np,1),9);
    Pnew(:,1:3) = np;
    % Constructing a noramalized vector in boundary in the frame of grain P
    a=[(-np(:,2)-np(:,3))./np(:,1), ones(size(np,1),1), ones(size(np,1),1)];
    len = sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2);
    a=a./len;

    Pnew(:,4:6) = a;

    b=cross(Pnew(:,1:3),Pnew(:,4:6));
    len = sqrt(b(:,1).^2 + b(:,2).^2 + b(:,3).^2);
    b=b./len;

    Pnew(:,7:9) = b;

    Qnew = [Pnew(:,1) * dg(1,1) + Pnew(:,2)*dg(1,2) + Pnew(:,3)*dg(1,3), ...
            Pnew(:,1) * dg(2,1) + Pnew(:,2)*dg(2,2) + Pnew(:,3)*dg(2,3), ...
            Pnew(:,1) * dg(3,1) + Pnew(:,2)*dg(3,2) + Pnew(:,3)*dg(3,3), ...
            Pnew(:,4) * dg(1,1) + Pnew(:,5)*dg(1,2) + Pnew(:,6)*dg(1,3),...
            Pnew(:,4) * dg(2,1) + Pnew(:,5)*dg(2,2) + Pnew(:,6)*dg(2,3),...
            Pnew(:,4) * dg(3,1) + Pnew(:,5)*dg(3,2) + Pnew(:,6)*dg(3,3),...
            Pnew(:,7) * dg(1,1) + Pnew(:,8)*dg(1,2) + Pnew(:,9)*dg(1,3),...
            Pnew(:,7) * dg(2,1) + Pnew(:,8)*dg(2,2) + Pnew(:,9)*dg(2,3),...
            Pnew(:,7) * dg(3,1) + Pnew(:,8)*dg(3,2) + Pnew(:,9)*dg(3,3)];
    parfor i=1:size(Qnew,1)
        P = reshape(Pnew(i,:),[3,3]);
        Q = reshape(Qnew(i,:),[3,3]);
        e(i) = GB5DOF(P',Q','Ni');
    end
end
