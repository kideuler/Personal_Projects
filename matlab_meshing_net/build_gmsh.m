function build_gmsh(sizes,quads,save)
if nargin < 1
    sizes = 20:5:100;
    quads = 0;
    save = true;
end

if ~contains(cd,'gmsh_api')
   cd gmsh_api 
end

H = randi([n-10,n+10]);
rho = rand(1,H) .* logspace(-.5,-5,H);
phi = rand(1,H) .* 2*pi;


r = ones(size(t));
for h=1:H
  r = r + rho(h)*sin(h*t+phi(h));  
end

param = [r.*cos(t), r.*sin(t)]/4 + [0.5,0.5];
curve.param = matlabFunction(param);

for n = sizes
    for m = 1:10
        create_mesh_cpp(n,m,quads,curve)
    end
end

% create makefile
fid = fopen('makefile','w');
fprintf(fid,'FC = g++\n');
fprintf(fid,'main:\n');
for n = sizes
    for m = 1:10
        fprintf(fid,'\t$(FC) random_%d_%d.cpp -Llib -lgmsh\n',n,m);
        fprintf(fid,'\t./a.out\n');
    end
end
fprintf(fid,'\t-rm -f *.cpp *.o *.out\n');

fprintf(fid,'.PHONY: clean\nclean:\n');
fprintf(fid,'\t-rm -f *.cpp *.o *.out meshes/*.msh meshes/*.mat meshes/*.cpp\n');

% this actually runs makefile
system('sudo make')


%gather necessary data
for n = sizes
    for m = 1:10
        [elems, xs] = read_gmsh_file(['TrainingData/random_',num2str(n),'_',num2str(m),'.msh']);
        sibhes = determine_sibling_halfedges(size(xs,1),elems);
        bnd = get_boundary(elems,xs,sibhes);
    end
end

system('rm -f meshes/*.msh');
end