%* *****************************************************************
%* - Function of STAPMAT in Dynamic Solver phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     To solve finite element static equilibrium equations        *
%*                                                                 *
%* - Call procedures:                                              *
%*     ./LDLTFactor.m            - LDLTFactor()                    *
%*     Solve.m                   - Stiff2Sparse()                  *
%*     ./ColSol.m                - ColSol()                        *  
%*     Solve.m                   - WriteDis()                      *
%*     SRC/Mechanics/GetStress.m - GetStress()                     *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%*                                                                 *
%* *****************************************************************

function DynSolve()

global cdata;
global sdata;

NEQ = sdata.NEQ;
NLCASE = cdata.NLCASE;
LINEAR = cdata.LINEAR;

% Applying Boundary Condition
    APPBC();

M = sdata.MASS;
K = sdata.STIFF;
Q = sdata.R;
c1 = 0; c2 = 0;
C = c1.*M + c2.*K;

a0 = sdata.INITX;
Da0 = sdata.INITV;
a0 = zeros(length(K),1);
Da0 = zeros(length(K),1);

rho = sdata.SPRHO;
intimestep = sdata.TSTEP;
time_limit = sdata.CTIME;
err = sdata.ERR;
Gamma1 = sdata.GAMMA1;
Gamma2 = sdata.GAMMA2;
Gamma3 = sdata.GAMMA3;

sdata.STRAIN = zeros(NEQ, NLCASE, 'double');
sdata.STRESS = zeros(NEQ, NLCASE, 'double');

cdata.TIM(4,:) = clock;

LINEAR = 1;
% Solve Dynamic Problem
if(LINEAR)     % Linear Problem
    fprintf("LINEAR Dynamic Solve is Run!\n");

    % 广义alpha法时间积分
    [Time,u,v,A]=Gene_Alpha_Integ_Linear(M,K,C,Q,a0,Da0,rho,intimestep,time_limit,err,Gamma1,Gamma2,Gamma3);
    sdata.DIS = u;
    sdata.OUTT = Time';
    sdata.OUTV = v;
    sdata.OUTA = A;

    % 后处理
    for t = 1:length(sdata.OUTT)
        % Compute Stress
        GetStress(1, t);
        
        % Output Paraview
        OutputVtu(t);

    end

else           % non-Linear Problem

end

cdata.TIM(5, :) = clock;

end % Function



