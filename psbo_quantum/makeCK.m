function UCK = makeCK(Hmat,dtau,X)
Imat = speye(length(X));
UCK=sparse(...                    % Midpoint Crank Nicolson
            (Imat-1i*0.5*dtau*Hmat)/...
            (Imat+1i*0.5*dtau*Hmat)); 
end