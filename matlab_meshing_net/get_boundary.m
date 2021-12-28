function [segs] = get_boundary(elems,xs,sibhes)

nv = size(xs,1);
nelems = size(elems,1);
segs = zeros(nelems,2);
E = [1,2;2,3;3,1];

% create nodal elem,edge adjacency lists
elem_nodes = zeros(nv,8);
forward_edge = zeros(nv,8);
for i = 1:nelems
   if nnz(elems(i,:)) == 3
      for j = 1:3
         k = elems(i,j);
         elem_nodes(k,nnz(elem_nodes(k,:))+1) = i;
         forward_edge(k,nnz(forward_edge(k,:))+1) = j;
      end
   end
end

% find first segment
for i = 1:nelems
    if nnz(elems(i,:)) == 3
       for j = 1:3
          if sibhes(i,j) == 0
              segs(1,:) = [elems(i,E(j,1)), elems(i,E(j,2))];
          end
       end
    end
end

k = 1;
curr = segs(1,2);
while curr ~= segs(1,1)
    curr = segs(k,2);
    nne = nnz(elem_nodes(curr,:));
    for z = 1:nne
       i = elem_nodes(curr,z);
       j = forward_edge(curr,z);
       if sibhes(i,j) == 0
          k=k+1;
          segs(k,:) = [elems(i,E(j,1)), elems(i,E(j,2))];
       end
    end
end

segs = segs(1:k-1,:);
end