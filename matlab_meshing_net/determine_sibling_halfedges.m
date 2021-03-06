function [sibhes, manifold, oriented] = determine_sibling_halfedges(nv, elems, varargin)
%DETERMINE_SIBLING_HALFEDGES constructs an extended half-edge data structure
%  for a non-oriented or non-manifold surface mesh. It works for both
%  triangle and quadrilateral meshes that are either linear and quadratic.
%
%    [SIBHES, MANIFOLD, ORIENTED] = DETERMINE_SIBLING_HALFEDGES(NV,ELEMS)
%
% Computes mapping from each half-edge to a sibling halfedge for a
% non-oriented or non-manifold surface mesh. The halfedges around an edge
% are in a cyclic order of the elements.
% At output, it returns SIBHES, and also logical variable indicating
% whether the mesh is a manifold, and if so whether it is oriented.
%
% Usage:
%  The following functionscalls set up the extended half-edge data structure:
%  >> sibhes = determine_sibling_halfedges(nv, elems);
%  >> v2he = determine_incident_halfedges(nv, elems, sibhes);
%
%  Then
%  >> heid_1 = v2he(v)
%  gives a half-edge incident on vertex v.
%  >> heid_2 = sibhes(heid2fid(heid_1), heid2leid(heid_1))
%  gives the next half-edge incident around the edge,
%  >> heid_3 = sibhes(heid2fid(heid_2), heid2leid(heid_2))
%  gives the third half-edge incident around the edge.
%  Repeat until heid_n == heid_1.
%
% See also DETERMINE_INCIDENT_HALFEDGES.

%#codegen -args {int32(0), coder.typeof(int32(0), [inf, inf])}
%#codegen determine_sibling_halfedges_usestruct -args
%#codegen {int32(0), coder.typeof(int32(0), [inf,inf]), false}

% Convention: Each half-edge is identified by <face_id,local_edge_id>.
%             We assign 2 bits to local_edge_id.

nvpE = int32(size(elems,2));  % Number of vertices per element
%{
if nvpE==4 || nvpE==8 || nvpE==9
    nepE = int32(4); % Number of edges per element
else
    assert(nvpE==3 || nvpE==6);
    nepE = int32(3); % Number of edges per element
end
%}

%%%%%% TESTING %%%%%
if nvpE==4 || nvpE==8 || nvpE==9 || nvpE==16 || nvpE==25 || nvpE == 36
    nepE = int32(4); % Number of edges per element
else
    assert(nvpE==3 || nvpE==6 || nvpE==10 || nvpE == 15);
    nepE = int32(3); % Number of edges per element
end


next = int32([2,3,4,1]);

%% First, build is_index to store starting position for each vertex.
is_index = zeros(nv+1,1,'int32');
nelems = int32(size(elems,1));
for ii=1:nelems
    if elems(ii,1)==0; nelems = ii-1; break; end
    
    hasthree = nepE~=4 || ~elems(ii,4);
    for j=1:4-int32(hasthree)
        k = elems(ii,j)+1; is_index(k) = is_index(k) + 1;
    end
end

is_index(1) = 1;
for ii=1:nv; is_index(ii+1) = is_index(ii) + is_index(ii+1); end

ne = nelems*nepE;
v2nv = coder.nullcopy(zeros(ne,1,'int32'));  % Vertex to next vertex in each halfedge.
v2he_fid = coder.nullcopy(zeros(ne,1,'int32'));  % Vertex to halfedge.
v2he_leid = coder.nullcopy(zeros(ne,1,'int8'));  % Vertex to halfedge.

for ii=1:nelems
    hasthree = nepE~=4 || ~elems(ii,4);
    for j=1:4-int32(hasthree)
        k = elems(ii,j);
        
        v2nv(is_index(k)) = elems(ii,next(j+int32(hasthree & j==3)));
        v2he_fid(is_index(k)) = ii;  % Vertex to halfedge.
        v2he_leid(is_index(k)) = j;  % Vertex to halfedge.
        is_index(k) = is_index(k) + 1;
    end
end
for ii=nv-1:-1:1; is_index(ii+1) = is_index(ii); end
is_index(1)=1;

%% Set sibhes
if nargin<3 || isempty(varargin{1}) || ~islogical(varargin{1})
    sibhes = zeros(nelems, nepE, 'int32');
elseif islogical(varargin{1})
    sibhes = struct('fid', zeros(size(elems), 'int32'), ...
        'leid', zeros(size(elems), 'int8'));
else
    sibhes = varargin{1};
    assert(size(sibhfs,1)>=nelems && size(sibhfs,2)>=nepE);
    sibhes(:) = int32(0);
end
    
manifold = true; oriented = true;

for ii=int32(1):nelems
    hasthree = nepE~=4 || ~elems(ii,4);
    for jj=1:4-int32(hasthree)
        % Process each edge only once
        if ~isstruct(sibhes) && sibhes(ii,jj) || ...
                isstruct(sibhes) && sibhes.fid(ii,jj)
            continue;
        end
        v = elems(ii,jj); vn = elems(ii,next(jj+int32(hasthree & jj==3)));
                
        first_heid_fid=ii;
        first_heid_leid=jj;
                
        prev_heid_fid=ii;
        prev_heid_leid=jj;
        
        nhes = int32(0);
        
        % LOCATE: Locate index in v2nv(first:last)
        for index = is_index(vn):is_index(vn+1)-1
            if v2nv(index)==v
                if ~isstruct(sibhes)
                    sibhes(prev_heid_fid,prev_heid_leid) = fleids2heid(v2he_fid(index), int32(v2he_leid(index)));
                else
                    sibhes.fid(prev_heid_fid,prev_heid_leid) = v2he_fid(index);
                    sibhes.leid(prev_heid_fid,prev_heid_leid) = v2he_leid(index);
                end
                                
                prev_heid_fid=v2he_fid(index);
                prev_heid_leid=int32(v2he_leid(index));
                
                nhes = nhes+1;
            end
        end
        
        % Check for halfedges in the same orientation
        for index = is_index(v):is_index(v+1)-1
            if v2nv(index)==vn && v2he_fid(index)~=ii
                if ~isstruct(sibhes)
                    sibhes(prev_heid_fid,prev_heid_leid) = fleids2heid(v2he_fid(index), int32(v2he_leid(index)));
                else
                    sibhes.fid(prev_heid_fid,prev_heid_leid) = v2he_fid(index);
                    sibhes.leid(prev_heid_fid,prev_heid_leid) = v2he_leid(index);
                end
                
                prev_heid_fid=v2he_fid(index);
                prev_heid_leid=int32(v2he_leid(index));
                
                nhes = nhes+1;
                oriented = false;
            end
        end
        
        if prev_heid_fid ~= first_heid_fid
            % Close up the cycle
            if ~isstruct(sibhes)
                sibhes(prev_heid_fid,prev_heid_leid) = fleids2heid(first_heid_fid, first_heid_leid);
            else
                sibhes.fid(prev_heid_fid,prev_heid_leid) = first_heid_fid;
                sibhes.leid(prev_heid_fid,prev_heid_leid) = first_heid_leid;
            end
            nhes = nhes+1;
        end
        
        if nargout>1 && manifold && nhes>2
            manifold = false; oriented = false;
        end
    end
end