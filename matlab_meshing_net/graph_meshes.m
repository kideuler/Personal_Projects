function graph_meshes(sizes,save_file,quads,curve)
for n = sizes
   if quads
      [elems, xs] = read_gmsh_file(['meshes/Ellipse_hole_',num2str(n),'_quads.msh']);
      [elems, xs, t] = gmsh_fix(elems,xs,quads,curve);
   else
      [elems, xs] = read_gmsh_file(['meshes/Ellipse_hole_',num2str(n),'.msh']); 
      [elems, xs, t] = gmsh_fix(elems,xs,quads,curve);
   end
   
   if n == 20
   if ~ishandle(1)
        h=figure(1);
   else
        h=clf('reset');
   end
   hold on
   graph_hybrid(elems,xs)
   fimplicit(curve.implicit)
   axis equal
   pause(0.01);
   end
   if save_file
       if quads
           save(['meshes/gmsh_ellipse_quad_n',num2str(n),'.mat'],'elems','xs','t');
           system(['git add meshes/gmsh_ellipse_quad_n',num2str(n),'.mat']);
       else
           save(['meshes/gmsh_ellipse_tri_n',num2str(n),'.mat'],'elems','xs','t');
           system(['git add meshes/gmsh_ellipse_tri_n',num2str(n),'.mat']);
       end
   end
   disp(['n = ',num2str(n)]);
end
end