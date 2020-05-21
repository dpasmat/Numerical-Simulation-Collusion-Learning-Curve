

clear

digitsold=digits(10);  

global eta gamma a b c_o c_n c_upn c_o2 f c_o

syms q_intau q_jntau q_iotau q_jotau q_intau1 q_jntau1 q_iotau1 q_jotau1 eta gamma a b c_o2 c_upn c_o c_n alpha lambda delta k f df

comp=zeros(99,2);

%for k=1:4
for dt=1:99

% Pre-allocate the size of each CELL;
bstar_=cell(1,2);
q_bstars_=cell(1,2);
prof_bstar_i_=cell(1,2);
eqn_t1=cell(1,2);
eqn_t2=cell(1,2);
eqn_t3=cell(1,2);
eqn_t4=cell(1,2);

% Assign parameter values for all coefficients; 
eta=2;
gamma=1;
b=480;
c_upn=40;
a=476;
c_o2=18;
alpha=0.99;
lambda=0.99;

%delta=0.90;
delta=sym(dt)/100

c_o=16;
c_n=20;

k=1;
f = @(x) (4*x)^-k + c_n/c_upn;
df = @(x) -k*(4*x)^(-k-1);

%g = piecewise(x>1,x^(-1) + c_n/c_upn,

% First order conditions; 

for ll=1:4

if ll==1
% o lambda, n alpha 
    eqn_t1{ll} = 2*eta*q_intau + eta*(1 + alpha)*q_jntau + 2*gamma*q_iotau + gamma*(1+lambda)*q_jotau + delta*c_upn*df(q_intau)*q_intau1 + 0 + 0 + 0 == b - c_upn;
    eqn_t2{ll} = eta*(1 + alpha)*q_intau + 2*eta*q_jntau + gamma*(1+lambda)*q_iotau + 2*gamma*q_jotau + 0 + delta*c_upn*df(q_jntau)*q_jntau1 + 0 + 0 == b - c_upn;   
    eqn_t3{ll} = 2*gamma*q_intau + gamma*(1 + alpha)*q_jntau + 2*eta*q_iotau + eta*(1+lambda)*q_jotau + 0 + 0 + 0 + 0 == a - c_o2;   
    eqn_t4{ll} = gamma*(1 + alpha)*q_intau + 2*gamma*q_jntau + eta*(1+lambda)*q_iotau + 2*eta*q_jotau + 0 + 0 + 0 + 0 == a - c_o2;
    eqn_t5{ll} = delta*(0 + 0 + 0 + 0 + 2*eta*q_intau1 + eta*(1 + alpha)*q_jntau1 + 2*gamma*q_iotau1 + gamma*(1+lambda)*q_jotau1 ) == delta*(b - c_upn*f(q_intau));
    eqn_t6{ll} = delta*(0 + 0 + 0 + 0 + eta*(1 + alpha)*q_intau1 + 2*eta*q_jntau1 + gamma*(1+lambda)*q_iotau1 + 2*gamma*q_jotau1 ) == delta*(b - c_upn*f(q_jntau));   
    eqn_t7{ll} = delta*(0 + 0 + 0 + 0 + 2*gamma*q_intau1 + gamma*(1 + alpha)*q_jntau1 + 2*eta*q_iotau1 + eta*(1+lambda)*q_jotau1 ) == delta*(a - c_o);   
    eqn_t8{ll} = delta*(0 + 0 + 0 + 0 + gamma*(1 + alpha)*q_intau1 + 2*gamma*q_jntau1 + eta*(1+lambda)*q_iotau1 + 2*eta*q_jotau1 ) == delta*(a - c_o);
end;


if ll==2    
% o comp, n alpha 
    eqn_t1{ll} = 2*eta*q_intau + eta*(1 + alpha)*q_jntau + 2*gamma*q_iotau + gamma*q_jotau + delta*c_upn*df(q_intau)*q_intau1 + 0 + 0 + 0 == b - c_upn;
        eqn_t2{ll} = q_jntau == double(q_fixjn1);          
    eqn_t3{ll} = 2*gamma*q_intau + gamma*(1 + alpha)*q_jntau + 2*eta*q_iotau + eta*q_jotau + 0 + 0 + 0 + 0 == a - c_o2;   
        eqn_t4{ll} = q_jotau == double(q_fixjo1);
    eqn_t5{ll} = delta*(0 + 0 + 0 + 0 + 2*eta*q_intau1 + eta*(1 + alpha)*q_jntau1 + 2*gamma*q_iotau1 + gamma*q_jotau1 ) == delta*(b - c_upn*f(q_intau));
    eqn_t6{ll} = delta*(0 + 0 + 0 + 0 + eta*(1 + alpha)*q_intau1 + 2*eta*q_jntau1 + gamma*q_iotau1 + 2*gamma*q_jotau1 ) == delta*(b - c_upn*f(q_jntau));   
    eqn_t7{ll} = delta*(0 + 0 + 0 + 0 + 2*gamma*q_intau1 + gamma*(1 + alpha)*q_jntau1 + 2*eta*q_iotau1 + eta*q_jotau1 ) == delta*(a - c_o);   
    eqn_t8{ll} = delta*(0 + 0 + 0 + 0 + gamma*(1 + alpha)*q_intau1 + 2*gamma*q_jntau1 + eta*q_iotau1 + 2*eta*q_jotau1 ) == delta*(a - c_o);
end;       
   

if ll==3
% o alpha, n lambda 
    eqn_t1{ll} = 2*eta*q_intau + eta*(1 + lambda)*q_jntau + 2*gamma*q_iotau + gamma*(1+alpha)*q_jotau + delta*c_upn*df(q_intau)*q_intau1 + 0 + 0 + 0 == b - c_upn;
    eqn_t2{ll} = eta*(1 + lambda)*q_intau + 2*eta*q_jntau + gamma*(1+alpha)*q_iotau + 2*gamma*q_jotau + 0 + delta*c_upn*df(q_jntau)*q_jntau1 + 0 + 0 == b - c_upn;   
    eqn_t3{ll} = 2*gamma*q_intau + gamma*(1 + lambda)*q_jntau + 2*eta*q_iotau + eta*(1+alpha)*q_jotau + 0 + 0 + 0 + 0 == a - c_o2;   
    eqn_t4{ll} = gamma*(1 + lambda)*q_intau + 2*gamma*q_jntau + eta*(1+alpha)*q_iotau + 2*eta*q_jotau + 0 + 0 + 0 + 0 == a - c_o2;
    eqn_t5{ll} = delta*(0 + 0 + 0 + 0 + 2*eta*q_intau1 + eta*(1 + lambda)*q_jntau1 + 2*gamma*q_iotau1 + gamma*(1+alpha)*q_jotau1 ) == delta*(b - c_upn*f(q_intau));
    eqn_t6{ll} = delta*(0 + 0 + 0 + 0 + eta*(1 + lambda)*q_intau1 + 2*eta*q_jntau1 + gamma*(1+alpha)*q_iotau1 + 2*gamma*q_jotau1 ) == delta*(b - c_upn*f(q_jntau));   
    eqn_t7{ll} = delta*(0 + 0 + 0 + 0 + 2*gamma*q_intau1 + gamma*(1 + lambda)*q_jntau1 + 2*eta*q_iotau1 + eta*(1+alpha)*q_jotau1 ) == delta*(a - c_o);   
    eqn_t8{ll} = delta*(0 + 0 + 0 + 0 + gamma*(1 + lambda)*q_intau1 + 2*gamma*q_jntau1 + eta*(1+alpha)*q_iotau1 + 2*eta*q_jotau1 ) == delta*(a - c_o);
end;


if ll==4
% o alpha, n comp 
    eqn_t1{ll} = 2*eta*q_intau + eta*q_jntau + 2*gamma*q_iotau + gamma*(1+alpha)*q_jotau + delta*c_upn*df(q_intau)*q_intau1 + 0 + 0 + 0 == b - c_upn;
        eqn_t2{ll} = q_jntau == double(q_fixjn3); 
    eqn_t3{ll} = 2*gamma*q_intau + gamma*q_jntau + 2*eta*q_iotau + eta*(1+alpha)*q_jotau + 0 + 0 + 0 + 0 == a - c_o2;   
        eqn_t4{ll} = q_jotau == double(q_fixjo3);     
    %eqn_t4{ll} = gamma*q_intau + 2*gamma*q_jntau + eta*(1+alpha)*q_iotau + 2*eta*q_jotau + 0 + 0 + 0 + 0 == a - c_o2;
    eqn_t5{ll} = delta*(0 + 0 + 0 + 0 + 2*eta*q_intau1 + eta*q_jntau1 + 2*gamma*q_iotau1 + gamma*(1+alpha)*q_jotau1 ) == delta*(b - c_upn*f(q_intau));
    eqn_t6{ll} = delta*(0 + 0 + 0 + 0 + eta*q_intau1 + 2*eta*q_jntau1 + gamma*(1+alpha)*q_iotau1 + 2*gamma*q_jotau1 ) == delta*(b - c_upn*f(q_jntau));   
    eqn_t7{ll} = delta*(0 + 0 + 0 + 0 + 2*gamma*q_intau1 + gamma*q_jntau1 + 2*eta*q_iotau1 + eta*(1+alpha)*q_jotau1 ) == delta*(a - c_o);   
    eqn_t8{ll} = delta*(0 + 0 + 0 + 0 + gamma*q_intau1 + 2*gamma*q_jntau1 + eta*(1+alpha)*q_iotau1 + 2*eta*q_jotau1 ) == delta*(a - c_o);
end;


% Solve for vector of eqlm outputs; 
    assume(q_iotau1>0 & q_jotau1>0 & q_intau1>0 & q_jntau1>0 & q_iotau>0 & q_jotau>0 & q_intau>1 & q_jntau>1   )
    bstar_{ll} = solve([eqn_t1{ll}, eqn_t2{ll}, eqn_t3{ll}, eqn_t4{ll}, eqn_t5{ll}, eqn_t6{ll}, eqn_t7{ll}, eqn_t8{ll}], [q_intau, q_jntau, q_iotau, q_jotau, q_intau1, q_jntau1, q_iotau1, q_jotau1], 'Real', true); %, 'Real', true 
    q_bstars_{ll} = [bstar_{ll}.q_intau; bstar_{ll}.q_jntau; bstar_{ll}.q_iotau; bstar_{ll}.q_jotau; bstar_{ll}.q_intau1; bstar_{ll}.q_jntau1; bstar_{ll}.q_iotau1; bstar_{ll}.q_jotau1] 

    
% Save the cooperative eqlm value for j, at tau, in each of the 2 scenarios; 
    
    % o scenario; 
    if ll==1
    fix_scenario_o = q_bstars_{ll}*1;
    q_fixjo1 = fix_scenario_o(4,1);
    q_fixjn1 = fix_scenario_o(2,1);
    end;

    % n scenario;     
    if ll==3
    fix_scenario_n = q_bstars_{ll}*1;
    q_fixjn3 = fix_scenario_n(2,1);
    q_fixjo3 = fix_scenario_n(4,1);    
    end;  
    
% Obtain profits from eqlm outputs;

    prof_bstar_i_{ll} = prof_taus(q_bstars_{ll})

end;

% Display eqlm profits;

celldisp(prof_bstar_i_);

% Compute Delta-Pi terms;
comp(dt,:) = double(prof_comp_t(prof_bstar_i_));

end; 

%x=0.90:0.91:0.1;
%plot(x, comp);