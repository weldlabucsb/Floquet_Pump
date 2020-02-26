function Hmatrix= makeHmatrix(k,depth)
%HMATRIX Summary of this function goes here
%   Detailed explanation goes here
% paramsBANDS;
numStates=101;%Number of states
B=depth;


%% First Derivative
C=0;
for ii=2:2:numStates
   C=[C; ii; 0]; 
end
C(end)=[];
D=zeros(numStates,1);
M1=gallery('tridiag',C,D,-C);


%% Second derivative matrix
M2=0;
for ii=2:2:(numStates-1)
    M2=[M2 ii^2 ii^2];
end
M2=sparse(-diag(M2));

%% Cosine Matrix
rowsL=[3:1:numStates];
colsL=[1:1:numStates-2];
valsL=[sqrt(2) ones(1,numStates-3)];
M3=sparse(rowsL,colsL,valsL,numStates,numStates)+...
   sparse(colsL,rowsL,valsL,numStates,numStates);

%% Put it together


Hmatrix=(k^2*speye(numStates)-2*1i*k*M1-M2)-B/4*M3;



end

