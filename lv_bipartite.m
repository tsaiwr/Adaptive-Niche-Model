function lv =lv_bipartite(t,N,SA,SP,rho,betaA,betaP,gammaAP,gammaPA,h)  

% dimensions: rho: (SA+SP)X1, betaA: SAXSA, betaP: SPXSP, gammaAP: SAXSP, gammaPA: SPXSA

lv = zeros(SA+SP,1);
num = zeros(SA+SP,1);
den = zeros(SA+SP,1);
thetaAP = logical(gammaAP);
thetaPA = thetaAP';
num(1:SA) = gammaAP*N(SA+1:SA+SP);
num(SA+1:SA+SP) = gammaPA*N(1:SA);
den(1:SA) = 1+h*thetaAP*N(SA+1:SA+SP);
den(SA+1:SA+SP) = 1+h*thetaPA*N(1:SA);
numden = num./den;
lv(1:SA) = N(1:SA).*(rho(1:SA) - betaA*N(1:SA) + numden(1:SA));
lv(SA+1:SA+SP) = N(SA+1:SA+SP).*(rho(SA+1:SA+SP) - betaP*N(SA+1:SA+SP) + numden(SA+1:SA+SP));
