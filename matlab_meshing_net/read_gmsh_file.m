function [elems, xs] = read_gmsh_file(filename)
if nargin < 1
   filename = 'gmsh_api/meshes/Flower_hole_20.msh'; 
end
fid = fopen(filename);

level = 0;
while true
    text = fgetl(fid);
    if (text==-1)
        break;
    end
    
    if level==0
        if strcmp(text,'$Nodes')
            level = 1;
            n = 0;  
            text = fgetl(fid);
            [val,~] = sscanf(text, '%d %d %d %d');
            nv = val(2);
            xs = zeros(nv,2);
        end
    elseif level == 1
        if strcmp(text, '$EndNodes')
            break;
        else
            [val,~] = sscanf(text, '%f');
            if length(val)==3
                n=n+1;
                xs(n,:) = [val(1), val(2)];
            end
        end
    end
end

while true
   text = fgetl(fid);
   if text==-1
      break; 
   end
   
   if level == 1
      if strcmp(text,'$Elements')
         level = 2;
         n = 0;
         text = fgetl(fid);
         [val,~] = sscanf(text, '%d');
         nelems = val(2);
         elems = zeros(nelems,4);
         nz = length(val);
      end
   elseif level == 2
       if strcmp(text, '$EndElements')
           break;
       else
          [val,~] = sscanf(text, '%d');
          if length(val) > 3 && nnz(val) > 3 && nz > 3
             n=n+1;
             if nnz(val) == 4
                elems(n,:) = [val(2), val(3), val(4), 0];
             elseif nnz(val) == 5
                 elems(n,:) = [val(2), val(3), val(4), val(5)];
             end
          end
          nz = length(val);
       end
   end
end
elems = elems(1:n,:);
end