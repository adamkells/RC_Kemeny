for i=1:size(R,1)
    for j=1:size(R,1)
        [i,j]
        cc(i,j)=0;
        for t=0:0.1:100
            tmp=expm(R*t);
            cc(i,j)=cc(i,j)+(tmp(i,j)*P_EQ(j)-P_EQ(i)*P_EQ(j));
        end
    end
end
