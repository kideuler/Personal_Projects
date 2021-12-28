function create_mesh_cpp(n,m,quads,curve)
if nargin < 2
   quads = false; 
end
if nargin < 1
    n = 20;
end
ls = 1/(n-1);
fid = fopen(['random_',num2str(n),'_',num2str(m),'.cpp'],'w');

fprintf(fid,'#include<gmsh.h>\n');
fprintf(fid,'int main(int argc, char **argv)\n{\n');
fprintf(fid,'\tdouble ls = %12.8f;\n\tgmsh::initialize(argc, argv);\n', ls);
t = [0:0.005:2*pi]';
F = curve.param(t);

sz = size(F,1);
for i = 1:sz
   fprintf(fid, '\tgmsh::model::geo::addPoint(%12.8f, %12.8f, 0, ls, %d);\n', F(i,1),F(i,2),i); 
end

%{
fprintf(fid, '\tgmsh::model::geo::addPoint(0, 0, 0, ls, %d);\n',sz+1);
fprintf(fid, '\tgmsh::model::geo::addPoint(1, 0, 0, ls, %d);\n',sz+2);
fprintf(fid, '\tgmsh::model::geo::addPoint(1, 1, 0, ls, %d);\n',sz+3);
fprintf(fid, '\tgmsh::model::geo::addPoint(0, 1, 0, ls, %d);\n\n',sz+4);
%}

%{
fprintf(fid, '\tgmsh::model::geo::addLine(%d,%d,1);\n',sz+1,sz+2);
fprintf(fid, '\tgmsh::model::geo::addLine(%d,%d,2);\n',sz+2,sz+3);
fprintf(fid, '\tgmsh::model::geo::addLine(%d,%d,3);\n',sz+3,sz+4);
fprintf(fid, '\tgmsh::model::geo::addLine(%d,%d,4);\n',sz+4,sz+1);
%}
fprintf(fid, '\tgmsh::model::geo::addSpline({');
for i = 1:sz
   fprintf(fid, '%d,',i);
end
fprintf(fid,'1},1);\n\n');

fprintf(fid,'\tgmsh::model::geo::addCurveLoop({5},1);\n');
%fprintf(fid,'\tgmsh::model::geo::addCurveLoop({1,2,3,4}, 2);\n\n');
fprintf(fid,'\tgmsh::model::geo::addPlaneSurface({1},1);\n\n');
fprintf(fid,'\tgmsh::model::geo::synchronize();\n');
if quads
    fprintf(fid,'\tgmsh::model::mesh::setAlgorithm(2,1,8);\n');
    fprintf(fid,'\tgmsh::model::mesh::setRecombine(2,1);\n');
end
fprintf(fid,'\tgmsh::model::mesh::generate(2);\n');
if quads
    fprintf(fid,['\tgmsh::write("TrainingData/Random_',num2str(n),'_',num2str(m),'_quads.msh");\n']);
else
    fprintf(fid,['\tgmsh::write("TrainingData/Random_',num2str(n),'_',num2str(m),'.msh");\n']);
end


fprintf(fid,'\tgmsh::finalize();\n');
fprintf(fid,'\treturn 0;\n};');

fclose(fid);
end