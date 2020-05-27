function [gammaAP,gammaPA] = wire_back(gammaAP,gammaPA,gamma_old,idstart,idend_new,idend_old,Type)

if Type==0 
    gammaAP(idstart,idend_new)=0;
    gammaPA(idend_new,idstart)=0;
    gammaAP(idstart,idend_old)=gamma_old;  
    gammaPA(idend_old,idstart)=gamma_old;  
else
    gammaPA(idstart,idend_new)=0;
    gammaAP(idend_new,idstart)=0;
    gammaPA(idstart,idend_old)=gamma_old;  
    gammaAP(idend_old,idstart)=gamma_old; 
end 