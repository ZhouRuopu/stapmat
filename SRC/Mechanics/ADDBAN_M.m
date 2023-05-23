%* *****************************************************************
%* - Function of STAPMAT in Mass phase                             *
%*                                                                 *
%* - Purpose:                                                      *
%*     To assemble element mass into global mass                   *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/ReadTruss.m - Assemble()                            *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     LeiYang Zhao, Yan Liu, Computational Dynamics Group,        *
%*     School of Aerospace Engineering, Tsinghua University,       *
%*     2019.02.22                                                  *
%*                                                                 *
%* *****************************************************************

function ADDBAN_M(M, LM)

% Get global data
global sdata;
MAXA = sdata.MAXA; 
MASS = sdata.MASS;
ND = sdata.NDOF * sdata.NNODE;
for J = 1:ND
    JJ = LM(J);
    if (JJ > 0)
        for I = 1:J
            II = LM(I);
            if (II > 0)
                if (JJ > II) KK = MAXA(JJ) + JJ - II;
                else KK = MAXA(II) + II - JJ; end
                MASS(KK) = MASS(KK) + M(I, J);
            end
        end
    end
end

sdata.MASS = MASS;

end