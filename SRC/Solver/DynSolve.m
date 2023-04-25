%* *****************************************************************
%* - Function of STAPMAT in Dynamic Solver phase                           *
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
sdata.DIS = zeros(NEQ, NLCASE, 'double');
sdata.STRAIN = zeros(NEQ, NLCASE, 'double');
sdata.STRESS = zeros(NEQ, NLCASE, 'double');

cdata.TIM(4,:) = clock;

TotalStep = ceil(ENDT/TSTEP);

% Solve Dynamic Problem
if(LINEAR)     % Linear Problem
    fprintf("LINEAR Dynamic Solve is Run!");
    CURT = 0;
    for i = 1:TotalStep
        CURT = CURT + i*TSTEP;
        % 广义alpha法时间积分
        Gene_Alpha_Integ();

        % 可视化输出
        View_Output();
        
        % 输出时间步，后续可设置命令行中输出部分变量以便观看运行情况
        if(i==1 || rem(i,3)==0 || i==TotalStep)   
        fprintf("step = %d;  T = %12.3f", i, CURT);
        end
    end

else           % non-Linear Problem

end

cdata.TIM(5, :) = clock;

end % Function



