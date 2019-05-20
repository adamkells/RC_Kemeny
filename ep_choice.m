function [end_points]=ep_choice(K,K_eig_R)

secvec=0;
if secvec==0 % MFPT endpoints
    MM=jjhunter(expm(K'));
    MM=MM+MM';
    maximum = max(max(MM));
    [x,y]=find(MM==maximum);
    end_points=[x,y];
elseif secvec==1 % 2nd eigenvector endpoints
    [~,b1] = min(K_eig_R(:,2));
    [~,b2] = max(K_eig_R(:,2));
    end_points=[b1,b2];
end

end