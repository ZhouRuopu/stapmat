%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function BeamStress(NUM, NG)

% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; AREA = sdata.AREA; LM = sdata.LM;
U = sdata.DIS(:, NUM);

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT        STRESS11_LN       STRESS12_LN       STRESS13_LN       STRESS12_RN       STRESS12_RN       STRESS12_RN\n' ...
    '       NUMBER\n'], NG);

for N = 1:NUME
    MTYPE = MATP(N);
   
%   compute the length of truss element
    DX = XYZ(1, N) - XYZ(4, N);
    DY = XYZ(2, N) - XYZ(5, N);
    DZ = XYZ(3, N) - XYZ(6, N);
    XL2 = DX*DX + DY*DY + DZ*DZ;
    XL = sqrt(XL2);
    
    epsx1 = 0.0;
    epsx2 = 0.0;
    gammaxy1 = 0.0;
    gammaxy2 = 0.0;
    gammaxz1 = 0.0;
    gammaxz2 = 0.0;

    for L = 1:12
        if LM(L,N) > 0
           du0dx = [1-1/XL, 1/XL]* [U(1,LM(L,N)), U(7,LM(L,N))]';
           d2w0dx1 = [6/XL2 + 12/XL2, 4/XL-6/XL, 6/XL2 - 12/XL2, -6/XL+2/XL]* [U(3,LM(L,N)), U(5,LM(L,N)), U(9,LM(L,N)), U(11,LM(L,N))]';
           d2v0dx1 = [6/XL2 + 12/XL2, 4/XL-6/XL, 6/XL2 - 12/XL2, -6/XL+2/XL]* [U(2,LM(L,N)), U(6,LM(L,N)), U(8,LM(L,N)), U(12,LM(L,N))]';
           epsx1 = epsx1 + du0dx - d2w0dx1 * XYZ(3, N) - d2v0dx1 * XYZ(2, N);

           d2w0dx2 = [6/XL2, 4/XL, 6/XL2, 2/XL]* [U(3,LM(L,N)), U(5,LM(L,N)), U(9,LM(L,N)), U(11,LM(L,N))]';
           d2v0dx2 = [6/XL2, 4/XL, 6/XL2, 2/XL]* [U(2,LM(L,N)), U(6,LM(L,N)), U(8,LM(L,N)), U(12,LM(L,N))]';
           epsx2 = epsx2 + du0dx - d2w0dx2 * XYZ(6, N) - d2v0dx2 * XYZ(5, N);

           dthetadx = [1-1/XL, 1/XL]* [U(4,LM(L,N)), U(10,LM(L,N))]';
           gammaxy1 = gammaxy1 - dthetadx*XYZ(3, N);
           gammaxy2 = gammaxy2 - dthetadx*XYZ(6, N);
           gammaxz1 = gammaxz1 + dthetadx*XYZ(2, N);
           gammaxz2 = gammaxz2 + dthetadx*XYZ(5, N);
        end
    end
    
    fprintf(IOUT, ' %10d :  %13.5e    %13.5e    %13.5e          %13.5e    %13.5e    %13.5e \n', N, epsx1, gammaxy1, gammaxz1, epsx2, gammaxy2, gammaxz2);
end

end