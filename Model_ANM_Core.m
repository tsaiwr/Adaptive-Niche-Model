% core of adaptive niche model
sd=rng('shuffle'); 
counter_max = 100000;  % total number of timesteps
to = 0;   
tf = 10;  % integration time within each time step
SA = 56; SP = 56;  % number of species in guilds A and P
S = SA + SP;
trSpan = 1;  % length of niche axis [0,1]
cM = 4/S^0.8;   % connectance

nw = 0.1;    % niche width
int_comp = 0.1;   % intensity of competition
int_mut = 0.1;   % intensity of mutualism 
h = 0.1;    % handling time 
eta = 1;

del_ts=1000;   % snapshots every del_ts steps
temp_snap=cell(1,round(counter_max/del_ts));

%%%%%%%%%%%%%%%%%%%%%%%%%
% intialization

counter=0;
cc=1;
idstart = 1;
Type = 0;
trait = trSpan*rand(S,1);
mu_N0 = 0.1;
mu_rho = 1;
N0 = mu_N0*ones(S,1);  % initial abundance
rho = mu_rho*ones(S,1);  % intrinsic growth rate
N = ones(S,1);

betaA = zeros(SA);
betaP = zeros(SP);

dA = 1;  % intraspecific competition intensity
dP = 1;

for i=1:SA
    for j=i+1:SA        
        betaA(i,j)=int_comp*overlap(trait(i),trait(j),nw,trSpan);   
        betaA(j,i)=betaA(i,j);        
    end
end
betaA = betaA + eye(SA)*dA;
for i=1:SP
    for j=i+1:SP        
        betaP(i,j)=int_comp*overlap(trait(SA+i),trait(SA+j),nw,trSpan);   
        betaP(j,i)=betaP(i,j);        
    end
end
betaP = betaP + eye(SP)*dP;

gammaAP=zeros(SA,SP);
gammaPA=zeros(SP,SA);

for i=1:SA
    for j=1:SP  
        if rand<cM;
            gammaAP(i,j)=int_mut*overlap(trait(i),trait(SA+j),nw,trSpan);   
            gammaPA(j,i)=gammaAP(i,j); 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive rewiring, updating niche proximities, integration

while counter<counter_max

    counter = counter+1;   

    % rewiring, updating niche proximities
    Type=randi(2)-1;  % rewiring in A or P guild                
    [gammaAP,gammaPA,idstart,idend_old,idend_new,gain_old_rw,gamma_old]=rewiring(Type,N0,SA,SP,gammaAP,gammaPA,trait,int_mut,nw,trSpan,eta);

    % integration
    [t N] = ode45(@(t,N)lv_bipartite(t,N,SA,SP,rho,betaA,betaP,gammaAP,gammaPA,h),[to to+tf],N0);
    to = t(end);
    N0=N(end,:)';  % set initial condition for next interval

    % link recovery if abundance does not increase
    if counter>1
        if Type==0
            gain_new_rw=N0(idstart); 
        else
            gain_new_rw=N0(idstart+SA); 
        end
        if gain_old_rw > gain_new_rw     % connect back to previous partner
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
        end
    end

    % output
    if mod(counter,del_ts)==0 && counter>0

        fprintf('\nModel_Core: timestep = %d\n\n', counter);

        if mod(counter,del_ts)==0  
            [nodf,qb,Nm] = cal_structure(gammaAP);
            fprintf('\nnodf = %f, modularity = %f\n\n', nodf, qb);
            fprintf('\n*************************************************\n');                
            temp_snap{cc}=gammaAP;  % network snapshots 
            cc=cc+1;
        end
    end
end

