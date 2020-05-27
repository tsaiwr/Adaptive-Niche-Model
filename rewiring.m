function [gammaXY,gammaYX,idstart,idend_old,idend_new,gain_old,gamma_old]=rewiring(Type,N,SX,SY,gammaXY,gammaYX,trait,mag_g,sigma_g,trSpan,eta)
% rewiring one link
if (Type==0)  % rewiring X->Y
    idX=randi(SX);
    while isempty(find(gammaXY(idX,:),1))  % if idX has no link to Y at all
        idX=randi(SX);
    end
    
    gain_old=N(idX);  
    
    idY_nz=find(gammaXY(idX,:));    
    idY_old=idY_nz(randi(length(idY_nz)));
    
    gamma_old=gammaXY(idX,idY_old);
    
    idY_new=idY_old;
    deg_old=nnz(gammaXY(:,idY_old));
    pr=1/deg_old^eta;  % rewiring probability 
    if rand>pr   
        idY_z=find(gammaXY(idX,:)==0);
        if ~isempty(idY_z)
        idY_new=idY_z(randi(length(idY_z)));
        % rewire to a new species
        gammaXY(idX,idY_old)=0;
        gammaYX(idY_old,idX)=0;
        gammaXY(idX,idY_new)=mag_g*overlap(trait(idX),trait(SX+idY_new),sigma_g,trSpan);
        gammaYX(idY_new,idX)=gammaXY(idX,idY_new);
        end
    end    
    
    idstart=idX;
    idend_old=idY_old;
    idend_new=idY_new;
       
else     % rewiring Y->X
    idY=randi(SY);
    while isempty(find(gammaYX(idY,:),1))  % if idX has no link to Y at all
        idY=randi(SY);
    end
    
    gain_old=N(idY+SX); 
    
    idX_nz=find(gammaYX(idY,:));
    idX_old=idX_nz(randi(length(idX_nz)));
    
    gamma_old=gammaYX(idY,idX_old);
    
    idX_new=idX_old;
    deg_old=nnz(gammaYX(:,idX_old));
    pr=1/deg_old^eta;
    if rand>pr   
        idX_z=find(gammaYX(idY,:)==0);
        if ~isempty(idX_z)        
        idX_new=idX_z(randi(length(idX_z)));
        % rewire to a new species
        gammaYX(idY,idX_old)=0;
        gammaXY(idX_old,idY)=0;
        gammaYX(idY,idX_new)=mag_g*overlap(trait(SX+idY),trait(idX_new),sigma_g,trSpan);
        gammaXY(idX_new,idY)=gammaYX(idY,idX_new);
        end
    end
    
    idstart=idY;
    idend_old=idX_old;
    idend_new=idX_new;
    
end