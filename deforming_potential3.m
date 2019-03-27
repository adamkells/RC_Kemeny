% Code by Adam Kells to use Kemeny constant for reaction coordinate
% identification
clear all
close all
% Create graph
% choose a number of edges
% Draw a random pair of states and connect them (if they are not already
% connected) until theres sufficient edges
c=0;
% two options either kemeny (1) or slowest eig (2)
for alpha=0:0.25:1
    for clus_opt=0:1
        c=c+1;
        N=200; %number of states
        x=linspace(0,5*pi,N);
        
        y1=2*sin((x-pi)); %double well
        y2=2*sin((x-pi)/2); %triple well
        y=alpha*y2+(1-alpha)*y1;
        
        y=y-min(y);
        A=10;
        KbT=0.596;
        K=zeros(N);
        for i=1:N-1
            K(i,i+1)=A*exp((y(i+1)-y(i))/2/KbT);
            K(i+1,i)=A*exp((y(i)-y(i+1))/2/KbT);
        end
        for i=1:N
            K(i,i)=-sum(K(:,i));
        end
        K=K';
        % do spectral decomposition
        [Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K');
        kemeny = sum(-1./Keigs(2:end));
        
        % what are the endpoints suggested by the second eigenvector
        [a,b1] = min(K_eig_R(:,2));
        [a,b2] = max(K_eig_R(:,2));
        
        % I want to choose two states at random and define a reaction coordinate
        % based on the commitor probability to reach one state over the other and
        % then search for the two states which maximise the kemeny constant
        % could maybe generalise to three stable states? i would need to consider
        % committor between each pair of states
        
        %%
        % lets choose two states
        tic
        end_points = randi(N,[1,2]);
        kem_max=0;
        num_clusters = 3;
        for i=1%:N-1
            [i]
            for j=N
                end_points=[i,j];
                
                % now for each other state I want to compute the commitor probability to
                % reach one state or the other first, details on what i'm doing are here:
                % www.emma-project.org/v2.2.1/api/generated/msmtools.analysis.committor.html
                [committor]=compute_commit(K,end_points);
                
                % the next bit is to now coarse grain the states based on
                % similarity of committor probability
                % do hummer-szabo clustering of rate matrix K
                kem_local = 0;
                [~,tmp2]=sort(committor);
                for i1=1:N-2
                    for j1=(i1+1):N-1
                        A=zeros(N,num_clusters);
                        A(tmp2(1:i1),1)=1;
                        A(tmp2(i1+1:j1),2)=1;
                        A(tmp2(j1+1:end),3)=1;
                        [R,P_EQ,Aclus]=hummer_szabo_clustering_A(K', eq, A);
                        
                        % analyse eigenvalues vectors of clustered matrix R
                        [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
                        
                        % compute reduced kemeny (for 2 state clustering this is just the
                        % same as the slowest timescale
                        if clus_opt==0
                            kemenyR=sum(-1./Reigs(2));
                        elseif clus_opt==1
                            kemenyR=sum(-1./Reigs(2:end));
                        end
                        
                        %  save the info for the iteration which maximises kemeny (this may
                        %  not be unique), will only save the first choice found
                        if kemenyR>kem_max
                            kem_max = kemenyR;
                            best_wells = end_points;
                            best_split = Aclus;
                        end
                        if kemenyR>kem_local
                            kem_local = kemenyR;
                            local_split(i,j,:,:)=Aclus;
                        end
                        
                    end
                end
                %save for eigenvector end points
                if i==min(b1,b2)
                    if j==max(b1,b2)
                        sec_eig_split = best_split;
                    end
                end
            end
        end
        kemeny_time=toc;
        
        %figure(c)
        subplot(5,2,c)
            
        plot(x(find(best_split(:,1))),y(find(best_split(:,1))),'r','linewidth',2)
        hold on
        plot(x(find(best_split(:,2))),y(find(best_split(:,2))),'b','linewidth',2)
        plot(x(find(best_split(:,3))),y(find(best_split(:,3))),'g','linewidth',2)
        if clus_opt==0
            title(['\alpha = ' num2str(alpha) ', slowest rel'])
        elseif clus_opt==1
            title(['\alpha = ' num2str(alpha) ', Kemeny'])
        end
    end
end