function build_gmsh(sizes,quads,save)
if nargin < 1
    sizes = 20:5:100;
    quads = 0;
    save = true;
end

if ~contains(cd,'gmsh_api')
   cd gmsh_api 
end

for n = sizes
for m = 1:10
    create_mesh_cpp(n,m,quads,Ellipse)
end
end

% create makefile
fid = fopen('makefile','w');
fprintf(fid,'FC = g++\n');
fprintf(fid,'main:\n');
for n = sizes
    fprintf(fid,'\t$(FC) Ellipse_Hole_%d.cpp -Llib -lgmsh\n',n);
    fprintf(fid,'\t./a.out\n');
end
fprintf(fid,'\t-rm -f *.cpp *.o *.out\n');

fprintf(fid,'.PHONY: clean\nclean:\n');
fprintf(fid,'\t-rm -f *.cpp *.o *.out meshes/*.msh meshes/*.mat meshes/*.cpp\n');

% this actually runs makefile
system('sudo make')

graph_meshes(sizes,save,quads,Ellipse);
system('rm -f meshes/*.msh');
end