for i=1:size(K,1)
    for j=1:size(K,1)
        [i,j]
        cd(i,j)=0;
        for t=0:0.1:100
            tmp=expm(K*t);
            cd(i,j)=cd(i,j)+(tmp(i,j)*eq(j)-eq(i)*eq(j));
        end
    end
end
