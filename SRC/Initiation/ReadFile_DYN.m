%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read input file of STAPMAT                                  *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Initiation/ReadFile.m - InitBasicData()                 *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function ReadFile_DYN(fname)
global fpath;
fname = strcat(fpath, fname);           % Deal the filename

% Get global class
global cdata;
global sdata;

% Open files
cdata.IIN = fopen(fname, 'r');

% Begin Read input file
fprintf('Input phase ...\n\n');

% the first time stamp
cdata.TIM = zeros(5, 6, 'double');
cdata.TIM(1,:) = clock;

IIN = cdata.IIN;
%% Read Control data
cdata.HED = fgetl(IIN);

tmp = str2num(fgetl(IIN));
cdata.NUMNP = int64(tmp(1));
cdata.NUMEG = int64(tmp(2));
cdata.NLCASE = int64(tmp(3));
cdata.MODEX = int64(tmp(4));

if (cdata.NUMNP == 0) return; end

%% Read nodal point data
InitBasicData();
% Define local variables to speed
ID = sdata.ID; X = sdata.X; Y = sdata.Y; Z = sdata.Z;
for i = 1:cdata.NUMNP
    tmp = str2num(fgetl(IIN));
    ID(1, i) = int64(tmp(2));
    ID(2, i) = int64(tmp(3));
    ID(3, i) = int64(tmp(4));
    ID(4, i) = int64(tmp(5));
    ID(5, i) = int64(tmp(6));
    ID(6, i) = int64(tmp(7));
    X(i) = double(tmp(8));
    Y(i) = double(tmp(9));
    Z(i) = double(tmp(10));
end
sdata.ID = ID; sdata.X = X; sdata.Y = Y; sdata.Z = Z;
sdata.BC = ID; sdata.XX = X; sdata.YY = Y; sdata.ZZ = Z;
%% Compute the number of equations
sdata.IDOrigin = ID;
NEQ = 0;
for N=1:cdata.NUMNP
    for I=1:3
        if (ID(I,N) == 0)
            NEQ = NEQ + 1;
            ID(I,N) = NEQ;
        else
            ID(I,N) = 0;
        end
    end
end
sdata.ID = ID;
sdata.NEQ = NEQ;
%% Read load data
% Init control data
NLCASE = cdata.NLCASE;
NUMNP = cdata.NUMNP; NDOF = 6;
sdata.R = zeros(NUMNP*NDOF, NLCASE, 'double');
R = sdata.R;
% Read data
for N = 1:cdata.NLCASE
    tmp = str2num(fgetl(IIN));
    cdata.LL = int64(tmp(1)); cdata.NLOAD = int64(tmp(2));
    NLOAD = cdata.NLOAD;
%   Init load data
    sdata.NOD = zeros(NLOAD, 1, 'int64');
    sdata.IDIRN = zeros(NLOAD, 1, 'int64');
    sdata.FLOAD = zeros(NLOAD, 1, 'double');
    NOD = sdata.NOD; IDIRN = sdata.IDIRN; FLOAD = sdata.FLOAD;
    
%   Read load data
    for I = 1:NLOAD
        tmp = str2num(fgetl(IIN));
        NOD(I) = int64(tmp(1));
        IDIRN(I) = int64(tmp(2));
        FLOAD(I) = double(tmp(3));
    end
    if (cdata.MODEX == 0) return; end
    
%   Compute load vector
    for L = 1:NLOAD
        % II = ID(IDIRN(L), NOD(L));
        % if (II > 0) R(II, N) = R(II, N) + FLOAD(L); end
        R(6*NOD(L)-6+IDIRN(L)) = R(6*NOD(L)-6+IDIRN(L))+FLOAD(L);
    end
    sdata.NOD = NOD; sdata.IDIRN = IDIRN; sdata.FLOAD = FLOAD; sdata.R = R;
end

%% Read Displacement
% % Init control data
% NLCASE = cdata.NLCASE;
% sdata.INITX = zeros(NEQ, NLCASE, 'double');
% INITX = sdata.INITX;
% % Read data
% for N = 1:cdata.NLCASE
%     tmp = str2num(fgetl(IIN));
%     cdata.LX = int64(tmp(1)); cdata.NDIS = int64(tmp(2));
%     NDIS = cdata.NDIS;
% %   Init load data
%     sdata.XNOD = zeros(NDIS, 1, 'int64');
%     sdata.XDIRN = zeros(NDIS, 1, 'int64');
%     sdata.XDIS = zeros(NDIS, 1, 'double');
%     XNOD = sdata.XNOD; XDIRN = sdata.XDIRN; XDIS = sdata.XDIS;
%     
% %   Read load data
%     for I = 1:NDIS
%         tmp = str2num(fgetl(IIN));
%         XNOD(I) = int64(tmp(1));
%         XDIRN(I) = int64(tmp(2));
%         XDIS(I) = double(tmp(3));
%     end
%     if (cdata.MODEX == 0) return; end
%     
% %   Compute load vector
%     for L = 1:NDIS
%         II = ID(XDIRN(L), XNOD(L));
%         if (II > 0) INITX(II, N) = INITX(II, N) + XDIS(L); end
%     end
%     sdata.XNOD = XNOD; sdata.XDIRN = XDIRN; sdata.XDIS = XDIS; sdata.INITX = INITX;
% end
% 
%% Read Velocity
% % Init control data
% NLCASE = cdata.NLCASE;
% sdata.INITV = zeros(NEQ, NLCASE, 'double');
% INITV = sdata.INITV;
% % Read data
% for N = 1:cdata.NLCASE
%     tmp = str2num(fgetl(IIN));
%     cdata.LV = int64(tmp(1)); cdata.NVEL = int64(tmp(2));
%     NVEL = cdata.NVEL;
% %   Init load data
%     sdata.VNOD = zeros(NVEL, 1, 'int64');
%     sdata.VDIRN = zeros(NVEL, 1, 'int64');
%     sdata.VVEL = zeros(NVEL, 1, 'double');
%     VNOD = sdata.VNOD; VDIRN = sdata.VDIRN; VVEL = sdata.VVEL;
%     
% %   Read load data
%     for I = 1:NVEL
%         tmp = str2num(fgetl(IIN));
%         VNOD(I) = int64(tmp(1));
%         VDIRN(I) = int64(tmp(2));
%         VVEL(I) = double(tmp(3));
%     end
%     if (cdata.MODEX == 0) return; end
%     
% %   Compute load vector
%     for L = 1:NVEL
%         II = ID(VDIRN(L), VNOD(L));
%         if (II > 0) INITV(II, N) = INITV(II, N) + VVEL(L); end
%     end
%     sdata.VNOD = VNOD; sdata.VDIRN = VDIRN; sdata.VVEL = VVEL; sdata.INITV = INITV;
% end

%% Read Time step
sdata.SPRHO = 0;
sdata.TSTEP = 0;
sdata.CTIME = 0;
tmp = str2num(fgetl(IIN));
sdata.SPRHO =double(tmp(1)); sdata.TSTEP = double(tmp(2)); sdata.CTIME = double(tmp(3));

%% Read ERROR estimation
sdata.ERR = 0;
sdata.GAMMA1 = 0;
sdata.GAMMA2 = 0;
sdata.GAMMA3 = 0;
tmp = str2num(fgetl(IIN));
sdata.ERR =double(tmp(1)); sdata.GAMMA1 = double(tmp(2)); sdata.GAMMA2 = double(tmp(3)); sdata.GAMMA3 = double(tmp(4));

end

%% Functions


% InitBasicData
function InitBasicData()
global cdata;
global sdata;

cdata.NPAR = zeros(10, 1, 'int64');

sdata.ID = zeros(6,cdata.NUMNP, 'int64');
sdata.X = zeros(cdata.NUMNP, 1, 'double');
sdata.Y = zeros(cdata.NUMNP, 1, 'double');
sdata.Z = zeros(cdata.NUMNP, 1, 'double');
end