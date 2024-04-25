function [e] = EulToEP(ori1,ori2, n11,n12,n13)
    e=zeros(length(n11),1);
    g1=[cos(ori1(1))*cos(ori1(3))-sin(ori1(1))*sin(ori1(3))*cos(ori1(2))   sin(ori1(1))*cos(ori1(3))+cos(ori1(1))*sin(ori1(3))*cos(ori1(2))  sin(ori1(3))*sin(ori1(2))
              -cos(ori1(1))*sin(ori1(3))-sin(ori1(1))*cos(ori1(3))*cos(ori1(2))  -sin(ori1(1))*sin(ori1(3))+cos(ori1(1))*cos(ori1(3))*cos(ori1(2))  cos(ori1(3))*sin(ori1(2))
                        sin(ori1(1))*sin(ori1(2))                                       -cos(ori1(1))*sin(ori1(2))                         cos(ori1(2))      ];
    g1=double(g1); 
        
    g2=[cos(ori2(1))*cos(ori2(3))-sin(ori2(1))*sin(ori2(3))*cos(ori2(2))   sin(ori2(1))*cos(ori2(3))+cos(ori2(1))*sin(ori2(3))*cos(ori2(2))  sin(ori2(3))*sin(ori2(2))
          -cos(ori2(1))*sin(ori2(3))-sin(ori2(1))*cos(ori2(3))*cos(ori2(2))  -sin(ori2(1))*sin(ori2(3))+cos(ori2(1))*cos(ori2(3))*cos(ori2(2))  cos(ori2(3))*sin(ori2(2))
                    sin(ori2(1))*sin(ori2(2))                                       -cos(ori2(1))*sin(ori2(2))                         cos(ori2(2))      ];
    g2=double(g2);
        
	O24 = [ 0, -1, 0; 1, 0, 0; 0, 0, 1];
	mm1 = norm(g1*g2'-eye(3));
	mm2 = norm(O24*g1*g2'-eye(3));
        if (mm1<5e-2 || mm2<5e-2)
            e(:,:) = zeros(length(n11),1) + 0.0624 + norm(g1*g2'-eye(3));
        else    
            n=[n11(:), n12(:), n13(:)];
            n1 = g1*n'; %normal in the frame of grain 1
            n2 = g2*n';
            e = NpNqDgToPQParallel(g2*g1',n1',n2');
        end
end

           
