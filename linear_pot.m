function [K,Adj,v] = linear_pot(nstates)

x=linspace(-4*pi,4*pi,nstates); % x-axis from -4pi to 4pi

% choice of four difference potentials

v =sin(0.7*(x - 1*pi)); %single well
v =-sin((x-pi)/2); %double well
v=0.5*sin(1.5*(x)/2-pi/2); %triple well
%y=-2*sin((x-pi)/2);
%y=linspace(1,1,length(x)); %flat potential
% y=y1+y2;
% Arrhenius rates: Create the rate matrix for the system
v=v-min(v);
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