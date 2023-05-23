function OutputVtu(i, U, MS)
global fname;
global sdata;
global cdata;

num_n = length(X);
% 创建vtu文件
i = num2string(i);
vtu_name = strcat(fname,'_',i,'.vtu');
fp = fopen(vtu_name,"w");

% 开始写入vtu文件
fprintf(fp,'<?xml version="1.0"?>\n');
fprintf(fp,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">\n');
fprintf(fp,'<UnstructuredGrid>\n');

   fprintf(fp,'   <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',num_n,NUME);
      % 写入节点信息
      fprintf(fp,'      <Points>\n');
         fprintf(fp,'         <DataArray type="Float32" NumberOfComponents="3" Format="ascii">\n');%节点坐标
         for v=1:num_n
             fprintf(fp,'         %3.5f   %3.5f   %3.5f\n',X(v)+U(v,1),Y(v)+U(v,2),Z(v)+U(v,3)); %位移矩阵需要转换后使用
         end
         fprintf(fp,'         </DataArray>\n');
      fprintf(fp,'      </Points>\n');
      % 写入节点值
      fprintf(fp,'      <PointData Scalars="Scalars">\n');
      % Displacement
         fprintf(fp,'          <DataArray type="Float32" Name="Disp" Format="ascii">\n');
         for d=1:num_n
             fprintf(fp,'         %3.5f\n',sqrt( U(d,1)*U(d,1)+U(d,2)*U(d,2)+U(d,3)*U(d,3) ) );
         end
         fprintf(fp,'          </DataArray>\n');
      % Mises Stress
         fprintf(fp,'          <DataArray type="Float32" Name="MStress" Format="ascii">\n');
         for s=1:num_n
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
