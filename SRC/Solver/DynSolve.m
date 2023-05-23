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
ENDT = cdata.ENDT;
TSTEP = cdata.TSTEP;

M = sdata.MASS;
K = sdata.STIFF;
Q = sdata.R;
C = zeros('like',K);

a0 = sdata.INITX;
Da0 = sdata.INITV;

rho = sdata.SPRHO;
intimestep = sdata.TSTEP;
time_limit = sdata.CTIME;
err = sdata.ERR;
Gamma1 = sdata.GAMMA1;
Gamma2 = sdata.GAMMA2;
Gamma3 = sdata.GAMMA3;


sdata.DIS = zeros(NEQ, NLCASE, 'double');
sdata.STRAIN = zeros(NEQ, NLCASE, 'double');
sdata.STRESS = zeros(NEQ, NLCASE, 'double');

cdata.TIM(4,:) = clock;


% Solve Dynamic Problem
if(LINEAR)     % Linear Problem
    fprintf("LINEAR Dynamic Solve is Run!");

        % 广义alpha法时间积分
        [Time,u,v,A]=Gene_Alpha_Integ_Linear(M,K,C,Q,a0,Da0,rho,intimestep,time_limit,err,Gamma1,Gamma2,Gamma3);


        % 可视化输出
        View_Output();

else           % non-Linear Problem

end

cdata.TIM(5, :) = clock;

end % Function



