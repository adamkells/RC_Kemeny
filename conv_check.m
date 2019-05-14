for ii=0:1
    for jj=[0,2]     
        [ii,jj]
        for choice=1:4
            %choice
            %keyboard
            [kemeny_l(choice),Rtmp]=kemeny_boundary(K,INV_K,eq,cc(:,:,choice),ii,jj);
        end
        if jj==2
            jj=1;
        end
        opt = (ii*3)+(jj+1)
        [a,b]=max(kemeny_l);
        b
        kemeny_l

    end
end