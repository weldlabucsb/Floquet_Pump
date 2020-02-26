function RK  = makeRK4(Hmat,dtau,X)
% Makes Runge Kutta 4 time evolution operator. Hmat is Hamiltonian. Tau is
% time step. X is mesh grid.
Imat = speye(length(X));
K=-1i*Hmat;
K1=Imat;
K2=(Imat+0.5*dtau*K);
K3=(Imat+0.5*dtau*K*K2);
K4=(Imat+dtau*K*K3);
k1=K*K1;
k2=K*K2;
k3=K*K3;
k4=K*K4;
RK =Imat+dtau/6*(k1+2*k2+2*k3+k4);

end