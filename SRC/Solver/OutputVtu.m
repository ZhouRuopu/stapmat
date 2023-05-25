%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/DynSolve.m                                       *
%*                                                                 *
%* - Programmed by:                                                *
%*     Rope Zhou                                                   *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* - 备注：各变量储存格式                                          *
%*   U(6,NUMNP);                                                   *
%*   STRESS(NUMNP, 6);                                             *
%*   STRESS(sig11, sig22, sig33, sig12, sig23, sig31);             *
%*                                                                 *
%*                                                                 *
%* *****************************************************************


function OutputVtu(t)
global fname;
global sdata;
global cdata;
global fpath;
fname_out = strcat(fpath, fname);  

NUMNP = cdata.NUMNP; NUME = sdata.NUME;
UOrign = sdata.DIS(:, t); STRAIN = sdata.STRAIN; STRESS = sdata.STRESS;
X = sdata.XX; Y = sdata.YY; Z = sdata.ZZ; IJ = sdata.IJ;
U = zeros(6, NUMNP); MS = zeros(NUMNP,1);   % von Mises Stress

for i = 1:NUMNP
    U(:,i) = UOrign(6*i-5:6*i);
    MS(i) = STRESS(i,1);      % 梁单元中Mises应力等于sig11，其余单元需相应修改；
end


% 创建vtu文件
t = num2str(t);
vtu_name = strcat(fname_out,'_',t,'.vtu');
fp = fopen(vtu_name,"w");

% 开始写入vtu文件
fprintf(fp,'<?xml version="1.0"?>\n');
fprintf(fp,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">\n');
fprintf(fp,'<UnstructuredGrid>\n');

   fprintf(fp,'   <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',NUMNP,NUME);
      % 写入节点信息
      fprintf(fp,'      <Points>\n');
         fprintf(fp,'         <DataArray type="Float32" NumberOfComponents="3" Format="ascii">\n');%节点坐标
         for v=1:NUMNP
             fprintf(fp,'         %3.5f   %3.5f   %3.5f\n',X(v)+U(1,v),Y(v)+U(2,v),Z(v)+U(3,v)); %位移矩阵需要转换后使用
         end
         fprintf(fp,'         </DataArray>\n');
      fprintf(fp,'      </Points>\n');
      % 写入节点值
      fprintf(fp,'      <PointData Scalars="Scalars">\n');
      % Displacement
         fprintf(fp,'          <DataArray type="Float32" Name="Disp" Format="ascii">\n');
         for d=1:NUMNP
             fprintf(fp,'         %3.5f\n',sqrt( U(1,d)*U(1,d)+U(2,d)*U(2,d)+U(3,d)*U(3,d) ) );
         end
         fprintf(fp,'          </DataArray>\n');
      % Mises Stress
         fprintf(fp,'          <DataArray type="Float32" Name="MStress" Format="ascii">\n');
         for s=1:NUMNP
             fprintf(fp,'         %3.5f\n', MS(s));
         end
         fprintf(fp,'          </DataArray>\n');
      fprintf(fp,'      </PointData>\n');

      % 写入单元值
      fprintf(fp,'      <Cells>\n');
         %Cells connectivity
         fprintf(fp,'         <DataArray type="Float32" Name="connectivity" Format="ascii">\n');
         for c=1:NUME
             fprintf(fp,'         %3.5f   %3.5f\n',IJ(c,1)-1,IJ(c,2)-1);
         end
         fprintf(fp,'         </DataArray>\n');
         %Cells offsets
         fprintf(fp,'         <DataArray type="Float32" Name="offsets" Format="ascii">\n');
         for c=1:NUME
             fprintf(fp,'         %d\n',2*c);
         end
         fprintf(fp,'         </DataArray>\n');
         %Cells types
         fprintf(fp,'         <DataArray type="Float32" Name="types" Format="ascii">\n');
         for c=1:NUME
             fprintf(fp,'         3\n');
         end
         fprintf(fp,'         </DataArray>\n');
      fprintf(fp,'      </Cells>\n');

   fprintf(fp,'   </Piece>\n');

fprintf(fp,'</UnstructuredGrid>\n');
fprintf(fp,'</VTKFile>\n');

fclose(fp); 
