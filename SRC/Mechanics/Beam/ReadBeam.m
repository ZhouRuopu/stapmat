%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of truss                       *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadBeam.m - ReadMaterial()                                *
%*     ReadBeam.m - ReadElements()                                *
%*                                                                 *
%* - Called by :                                                   *
%*     ./BeamStiff.m                                              *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function ReadBeam()

% Read Material information
ReadMaterial()

% Read Element information
ReadElements()

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end

% ----------------------- Functions -----------------------------------
% Read Material information
function ReadMaterial()

global cdata;
global sdata;
% Get file pointers
IIN = cdata.IIN;
IOUT = cdata.IOUT;

if (cdata.NPAR(3) == 0) cdata.NPAR(3) = 1; end
fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
fprintf(IOUT, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', ...
    cdata.NPAR(3));
fprintf(IOUT, '  SET       YOUNG''S     CROSS-SECTIONAL\n');
fprintf(IOUT, ' NUMBER     MODULUS          AREA\n');
fprintf(IOUT, '               E              A\n');


% Read material datas
sdata.NUME = cdata.NPAR(2);
sdata.NUMMAT = cdata.NPAR(3);
NUMMAT = cdata.NPAR(3);
sdata.E = zeros(NUMMAT, 1, 'double');
sdata.nu = zeros(NUMMAT, 1, 'double');
sdata.AREA = zeros(NUMMAT, 1, 'double');
sdata.RHO = zeros(NUMMAT, 1, 'double');
sdata.IY = zeros(NUMMAT, 1, 'double');
sdata.IZ = zeros(NUMMAT, 1, 'double');
sdata.IJ = zeros(NUMMAT, 1, 'double');
for I = 1:NUMMAT
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.nu(N) = tmp(3);
    sdata.RHO(N) = tmp(4);
    sdata.AREA(N) = tmp(5);
    sdata.IY(N) = tmp(6);
    sdata.IZ(N) = tmp(7);
    sdata.TJ(N) = tmp(8);
    fprintf(IOUT, '%5d    %12.5e    %12.5e  %14.6e  %14.6e    %12.5e    %12.5e    %12.5e\n', ...
                    N, tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8));
end

sdata.G = sdata.E/(2*(1+sdata.nu));

end

% Read elements information
function ReadElements()

global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n');
fprintf(IOUT, '\n      ELEMENT          NODE          NODE       MATERIAL\n');
fprintf(IOUT, '      NUMBER-N           I             J       SET NUMBER\n');

% Get Position data
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(6, NUME, 'double');
sdata.MATP = zeros(NUME, 1, 'int64');                 % the type of material
sdata.LM = zeros(12, NUME, 'double');                  % connectivity matrix
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
sdata.N1 = zeros(1, NUME, 'double');
sdata.N2 = zeros(1, NUME, 'double');
sdata.N3 = zeros(1, NUME, 'double');
X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
XYZ = sdata.XYZ; MATP = sdata.MATP; LM = sdata.LM;
N1 = sdata.N1; N2 = sdata.N2; N3 = sdata.N3;

for N = 1:NUME
    tmp = str2num(fgetl(IIN));
    I = round(tmp(2));
    J = round(tmp(3));
    N1(1, N) = tmp(4);
    N2(1, N) = tmp(5);
    N3(1, N) = tmp(6);
    MTYPE = round(tmp(7));
    
%   Save element information
    XYZ(1, N) = X(I);
    XYZ(2, N) = Y(I);
    XYZ(3, N) = Z(I);
    XYZ(4, N) = X(J);
    XYZ(5, N) = Y(J);
    XYZ(6, N) = Z(J);
    MATP(N) = MTYPE;

    IJ(N,1) = I;
    IJ(N,2) = J;
    
    fprintf(IOUT, '%10d      %10d    %10d       %5d\n', N, I, J, MTYPE);

%   Compute connectivity matrix
    for L=1:6
        LM(L,N) = ID(L,I);
        LM(L+6,N) = ID(L,J);
    end

%   Updata column heights and bandwidth
    ColHt(LM(:, N))
end
sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;
sdata.N1 = N1; sdata.N2 = N2; sdata.N3 = N3;
sdata.IJ = IJ;

% Clear the memory of X, Y, Z
sdata.X = double(0);
sdata.Y = double(0);
sdata.Z = double(0);

end