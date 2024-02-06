function [E, i] = Kepler(ECC, MA, tol)
% Function to solve Kepler's equation M = E - e * sinE
%
% Inputs:
%   ECC     : Eccentricity
%   MA      : Mean anomaly [rad] must be 0 <= M <= pi
%   tol     : Threshold tolerance
%
% Outputs:
%   E       : eccentric anomaly
%   i       : Iteration count

if ECC >= 0.99 | ECC < 0                    %Error exit if e is negative or greater than 0.99
    
    error('Input Error: the eccentricity must be contained within 0 < e < 1');
    return

else

    En = MA;                               %initialise algorithm with 'educated guess' 
    i = 0;                                 %initialise iteration count
    fEn = En - ECC * sin(En) - MA;         %initialise F(En)  

    while abs(fEn) > tol                  %iterate while F(En) is higher than threshold value
        
        i= i+1;                             %add to iteration counter 
        fEn = En - ECC * sin(En) - MA;     
        fpEn = 1 - ECC * cos(En);
        Ens = En - fEn/fpEn;                %Newton's method 
        En = Ens;                           %update En with iterated value

    end
end

E = mod(En, 2*pi);                                     %update E with final computed value
