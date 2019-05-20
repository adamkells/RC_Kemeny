function [K,Adj,v] = linear_pot(nstates)

x=linspace(-4*pi,4*pi,nstates/4); % x-axis from -4pi to 4pi

% choice of four difference potentials
% v1 =sin(0.3*(x - 1*pi)); %single well
% %v1 =-sin((x-pi)/2); %double well
% v2=2*sin(1.5*(x)/2-pi/2); %triple well
% v2=-2*sin(1.5*(x)-pi/2);
v=[];
for i=1:4
    v3=0.75*i*sin(0.25*(x)-1.5*pi);
    v3=v3-min(v3);
    v=[v v3];
end
v=-v;
v=v-min(v);
%y=-2*sin((x-pi)/2);
%y=linspace(1,1,length(x)); %flat potential
% v=v1+v2+x/(4*pi);
% % Arr%keyboard
A=1;
KbT=0.596;
for i=1:nstates-1
    K(i,i+1)=A*exp((v(i+1)-v(i))/2/KbT);
    K(i+1,i)=A*exp((v(i)-v(i+1))/2/KbT);
end
for i=1:nstates
    K(i,i)=-sum(K(:,i));
end

Adj=K./K;
Adj(isnan(Adj))=0;

end