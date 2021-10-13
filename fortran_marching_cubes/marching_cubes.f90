module marching_cubes
    use init_mesh
    use surface
    use utilities
    use Mesh_smoothing

    implicit NONE
    type mctable
        character(len=12) :: Edge(256)
        integer :: Tris(256,16)
    end type mctable

    type Elemedge
        integer :: elem(3)
        integer :: local(3)
    end type Elemedge
    contains

    function March_cubes(Cubes) result(Surf) ! main marching cubes algorithm
        implicit NONE
        type(mesh),intent(in) :: Cubes
        type(mesh) :: Surf
        type(mctable) :: table
        integer,allocatable :: edge_nodes(:,:),cubeindex(:),Loc_tri(:)
        integer,allocatable :: elem_edge(:,:)
        integer :: nelems, nv,i,E,kk,j
        character(12) :: str

        table = get_table() ! retrieving table for intersection cases
        cubeindex = find_cubeindex(cubes) ! finding index for the table
        allocate(edge_nodes(12,2))
        edge_nodes = transpose(reshape((/1,2, 2,3, 3,4, 4,1, 5,6, 6,7, 7,8, &
        8,5, 1,5, 2,6, 3,7, 4,8/),(/2,12/)))

        nelems = size(Cubes%elems,1) ! number of elements in the cubic mesh
        nv = size(Cubes%xs,1) ! number of nodes in the cubic mesh

        ! allocating arrays
        allocate(Surf%elems(nelems,3))
        allocate(Surf%xs(nv,3))
        allocate(elem_edge(nelems,12))
        allocate(Loc_tri(16))

        E=1
        kk=1
        do i=1,nelems
            str = table%Edge(cubeindex(i))
            ! Looping through edges
            do j = 1,12
                if (Ichar(str(13-j:j)) == 49) then
                    ! for each intersection on the cube the point is linearly interpolated onto the implicit function
                    Surf%xs(kk,:) = interp(0.0,cubes%xs(cubes%elems(i,edge_nodes(j,1)),:), &
                    cubes%xs(cubes%elems(i,edge_nodes(j,2)),:), &
                    isosurface(cubes%xs(cubes%elems(i,edge_nodes(j,1)),1),cubes%xs(cubes%elems(i,edge_nodes(j,1)),2), &
                    cubes%xs(cubes%elems(i,edge_nodes(j,1)),3)), isosurface(cubes%xs(cubes%elems(i,edge_nodes(j,2)),1), &
                    cubes%xs(cubes%elems(i,edge_nodes(j,2)),2),cubes%xs(cubes%elems(i,edge_nodes(j,2)),3)))
                    elem_edge(i,j) = kk
                    kk=kk+1
                end if
            end do
        end do
        Surf%xs = Surf%xs(1:(kk-1),:)

        do i = 1,nelems
            Loc_tri = table%Tris(cubeindex(i),:)
            do j = 1,count(Loc_tri/=0)/3   
                Surf%elems(E,1) = elem_edge(i,Loc_tri(3*(j-1)+1))
                Surf%elems(E,2) = elem_edge(i,Loc_tri(3*(j-1)+2))
                Surf%elems(E,3) = elem_edge(i,Loc_tri(3*(j-1)+3))
                E=E+1
            end do
        end do
        Surf%elems = Surf%elems(1:(E-1),:)

        deallocate(Loc_tri)
        deallocate(elem_edge)
        deallocate(edge_nodes)
        ! subroutine to remaove duplicate nodes, this way the mesh can be used in surface FEM.
        call remove_duplicate_nodes(Surf)
    end function March_cubes

    subroutine Areas(surf)
        implicit NONE
        type(mesh), intent(in) :: surf
        integer :: k,nelems
        real :: ps(3,3)
        real,allocatable :: A(:)

        nelems = size(surf%elems,1)
        allocate(A(nelems))
        do k = 1,nelems
            ps = surf%xs(surf%elems(k,:),:)
            A(k) = Area(ps)
        end do
        print *, 'min = ',minval(A)
        print *, 'max = ',maxval(A)
    end subroutine Areas

    subroutine remove_small_angle_tri(surf)
        implicit NONE
        type(mesh),intent(inout) :: surf
        integer, allocatable :: keep(:)
        integer :: nelems,k,kk=0,nv
        real:: angles(3),l1,l2
        real, PARAMETER :: pi = 3.1415927

        nv = size(surf%xs,1)
        nelems = size(surf%elems,1)
        allocate(keep(nelems))
        do k = 1,nelems
            angles(1) = 180*(1/pi)*Acos(dot_product(surf%xs(surf%elems(k,2),:)-surf%xs(surf%elems(k,1),:), &
            surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,1),:))/(norm2(surf%xs(surf%elems(k,2),:)-&
            surf%xs(surf%elems(k,1),:))*norm2(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,1),:))))

            angles(2) = 180*(1/pi)*Acos(dot_product(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,2),:), &
            surf%xs(surf%elems(k,1),:)-surf%xs(surf%elems(k,2),:))/(norm2(surf%xs(surf%elems(k,2),:)-&
            surf%xs(surf%elems(k,3),:))*norm2(surf%xs(surf%elems(k,2),:)-surf%xs(surf%elems(k,1),:))))

            angles(3) = 180*(1/pi)*Acos(dot_product(surf%xs(surf%elems(k,2),:)-surf%xs(surf%elems(k,3),:), &
            surf%xs(surf%elems(k,1),:)-surf%xs(surf%elems(k,3),:))/(norm2(surf%xs(surf%elems(k,3),:)-&
            surf%xs(surf%elems(k,1),:))*norm2(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,2),:))))

            if (angles(1)<3) then
                l1 = norm2(surf%xs(surf%elems(k,2),:)-surf%xs(surf%elems(k,1),:))
                l2 = norm2(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,1),:))
                if (l1<l2) then
                    surf%xs(surf%elems(k,2),:) = project(surf%xs(surf%elems(k,:),:),2,(/3,1/)) !project 2 onto 1,3
                else
                    surf%xs(surf%elems(k,3),:) = project(surf%xs(surf%elems(k,:),:),3,(/2,1/))!project 3 onto 1,2
                end if
            else if (angles(2)<3) then
                l1 = norm2(surf%xs(surf%elems(k,2),:)-surf%xs(surf%elems(k,1),:))
                l2 = norm2(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,2),:))
                if (l1<l2) then
                    surf%xs(surf%elems(k,1),:) = project(surf%xs(surf%elems(k,:),:),1,(/3,2/))!project 1 onto 2,3
                else
                    surf%xs(surf%elems(k,3),:) = project(surf%xs(surf%elems(k,:),:),3,(/1,2/))!project 3 onto 1,2
                end if
            else if (angles(3)<3) then
                l1 = norm2(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,1),:))
                l2 = norm2(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,2),:))
                if (l1<l2) then
                    surf%xs(surf%elems(k,1),:) = project(surf%xs(surf%elems(k,:),:),1,(/2,3/))!project 1 onto 2,3
                else
                    surf%xs(surf%elems(k,2),:) = project(surf%xs(surf%elems(k,:),:),2,(/1,3/))!project 2 onto 1,3
                end if
            else
                keep(kk+1) = k
                kk=kk+1
            end if 
        end do
        surf%elems = surf%elems(keep(1:kk),:)
        deallocate(keep)
    end subroutine remove_small_angle_tri

    function project(xs,n,line) result(P)
        implicit NONE
        real,intent(in) :: xs(3,3)
        integer, intent(in) :: n,line(2) !node being projected and line projcted on
        real :: P(3),v(3),M(3,3),P0(3)
        integer :: i,j

        v = xs(line(1),:)-xs(line(2),:)
        P0 = xs(line(2),:)
        do i = 1,3
            do j = 1,3
                M(i,j) = v(i)*v(j)
            end do
        end do
        M = M/dot_product(v,v)

        P = matmul(M,xs(n,:)-P0) + P0
    end function project

    subroutine remove_tris_edge(Surf)
        implicit NONE
        type(mesh),intent(inout) :: Surf
        integer :: nelems, nv,k
        real :: ave

        nelems = size(surf%elems,1)
        do k = 1,nelems
            ave = ave + norm2(surf%xs(surf%elems(k,1),:)-surf%xs(surf%elems(k,2),:))
            ave = ave + norm2(surf%xs(surf%elems(k,3),:)-surf%xs(surf%elems(k,2),:))
            ave = ave + norm2(surf%xs(surf%elems(k,1),:)-surf%xs(surf%elems(k,3),:))
        end do
        ave = ave/real(3*nelems)

        do k = 1,nelems
            if (norm2(surf%xs(surf%elems(k,1),:)-surf%xs(surf%elems(k,2),:)) < 0.01*ave) then
                surf%xs(surf%elems(k,2),:) = (surf%xs(surf%elems(k,2),:)+surf%xs(surf%elems(k,1),:))/2.0
                surf%xs(surf%elems(k,1),:) = (surf%xs(surf%elems(k,2),:)+surf%xs(surf%elems(k,1),:))/2.0
            end if
            if (norm2(surf%xs(surf%elems(k,2),:)-surf%xs(surf%elems(k,3),:)) < 0.01*ave) then
                surf%xs(surf%elems(k,2),:) = (surf%xs(surf%elems(k,2),:)+surf%xs(surf%elems(k,3),:))/2.0
                surf%xs(surf%elems(k,3),:) = (surf%xs(surf%elems(k,2),:)+surf%xs(surf%elems(k,3),:))/2.0
            end if
            if (norm2(surf%xs(surf%elems(k,1),:)-surf%xs(surf%elems(k,3),:)) < 0.01*ave) then
                surf%xs(surf%elems(k,3),:) = (surf%xs(surf%elems(k,3),:)+surf%xs(surf%elems(k,1),:))/2.0
                surf%xs(surf%elems(k,1),:) = (surf%xs(surf%elems(k,3),:)+surf%xs(surf%elems(k,1),:))/2.0
            end if
        end do

        call remove_duplicate_nodes(surf)
        call remove_degenerate_elems(surf)
    end subroutine remove_tris_edge

    subroutine remove_degenerate_elems(Surf)
        implicit NONE
        type(mesh),intent(inout) :: Surf
        integer :: nelems,n,kk=0
        integer,allocatable :: keep(:)

        nelems = size(Surf%elems,1)
        allocate(keep(nelems))
        do n = 1,nelems
            if (count(surf%elems(n,:)==surf%elems(n,1))==1 .and. count(surf%elems(n,:)==surf%elems(n,2))==1) then
                keep(kk+1) = n
                kk=kk+1
            end if
        end do
        surf%elems = surf%elems(keep(1:kk),:)
    end subroutine remove_degenerate_elems

    subroutine remove_duplicate_nodes(Surf)
        implicit none
        type(mesh),intent(inout) :: Surf
        integer, allocatable :: map(:),unique(:),o2n(:),nodes_id(:)
        integer :: nelems, nv, kk,n,j,i
        logical :: flag = .false.
        logical, allocatable :: nodes(:)

        nelems = size(Surf%elems,1)
        nv = size(Surf%xs,1)

        allocate(map(nv))
        allocate(unique(nv))

        ! find which nodes are duplicates
        unique(1) = 1
        map(1) =1
        kk=1
        do n = 2,nv
            flag = .false.
            do j = 1,kk
                if (norm2(surf%xs(n,:)-surf%xs(unique(j),:))<0.0001) then
                    map(n) = unique(j)
                    flag = .true.
                    exit
                end if
            end do
            if (flag .eqv. .false.) then
                unique(kk+1) = n
                map(n) = n
                kk=kk+1
            end if
        end do

        ! map duplicate nodes to old ones
        do i = 1,nelems
            do j = 1,3
                surf%elems(i,j) = map(surf%elems(i,j))
            end do
        end do

        ! remove unused nodes
        allocate(o2n(nv))
        allocate(nodes(nv))
        allocate(nodes_id(nv))
        do i = 1,nelems
            do j = 1,3
                nodes(surf%elems(i,j)) = .true.
            end do
        end do
        kk = 1
        do n = 1,nv
            if (nodes(n) .eqv. .true.) then
                nodes_id(kk) = n
                o2n(n) = kk
                kk=kk+1
            end if
        end do
        nodes_id = nodes_id(1:(kk-1))

        surf%xs = surf%xs(nodes_id,:)
        do i = 1,nelems
            do j = 1,3
                surf%elems(i,j) = o2n(surf%elems(i,j))
            end do
        end do



        deallocate(map)
        deallocate(unique)

        deallocate(o2n)
        deallocate(nodes)
        deallocate(nodes_id)
    end subroutine remove_duplicate_nodes

    function find_edges(E,edge,nn) result(conn_edge)
        implicit NONE
        integer, intent(in) :: E,edge,nn
        type(Elemedge) :: conn_edge

        if (edge==1) then
            conn_edge%elem = (/E-1,E-nn-1,E-nn/)
            conn_edge%local = (/5,7,3/)
            return
        else if (edge==2) then
            conn_edge%elem = (/E-1,E-nn**2-1,E-nn**2/)
            conn_edge%local = (/6,8,4/)
            return
        else if (edge==3) then
            conn_edge%elem = (/E-1,E+nn-1,E+nn/)
            conn_edge%local = (/7,5,1/)
            return
        else if (edge==4) then
            conn_edge%elem = (/E-1,E+nn**2-1,E+nn**2/)
            conn_edge%local = (/8,6,2/)
            return
        else if (edge==5) then
            conn_edge%elem = (/E+1,E-nn+1,E-nn/)
            conn_edge%local = (/1,3,7/)
            return
        else if (edge==6) then
            conn_edge%elem = (/E+1,E-nn**2+1,E-nn**2/)
            conn_edge%local = (/2,4,8/)
            return
        else if (edge==7) then
            conn_edge%elem = (/E+1,E+nn+1,E+nn/)
            conn_edge%local = (/3,1,5/)
            return
        else if (edge==8) then
            conn_edge%elem = (/E+1,E+nn+1,E+nn/)
            conn_edge%local = (/4,2,6/)
            return
        else if (edge==9) then
            conn_edge%elem = (/E-nn,E+nn**2,E+nn**2-nn/)
            conn_edge%local = (/12,10,11/)
            return
        else if (edge==10) then
            conn_edge%elem = (/E-nn,E-nn**2,E-nn**2-nn/)
            conn_edge%local = (/11,9,12/)
            return
        else if (edge==11) then
            conn_edge%elem = (/E+nn,E-nn**2,E-nn**2+nn/)
            conn_edge%local = (/10,12,9/)
            return
        else if (edge==12) then
            conn_edge%elem = (/E+nn,E+nn**2,E+nn**2+nn/)
            conn_edge%local = (/9,11,10/)
            return
        end if
    end function find_edges
    

    function interp(isolevel,xs1,xs2,v1,v2) result(xsn)
        real,intent(in)::isolevel,v1,v2,xs1(3),xs2(3)
        real :: xsn(3),mu
        if (abs(isolevel-v1)<0.00001) then
            xsn = xs1
            return
        else if (abs(isolevel-v2)<0.00001) then
            xsn = xs2
            return
        else if (abs(v1-v2)<0.00001) then
            xsn=xs1
            return
        else
            mu = (isolevel-v1)/(v2-v1)
            xsn = xs1 + mu*(xs2-xs1)
            return
        end if
    end function interp

    function find_cubeindex(cubes) result(cubeindex)
        implicit NONE
        type(mesh), intent(in) :: cubes
        integer, allocatable :: cubeindex(:)
        integer :: nelems,k

        nelems = size(cubes%elems,1)

        allocate(cubeindex(nelems))
        do k = 1,nelems
            cubeindex(k) = 1
            if (isosurface(cubes%xs(cubes%elems(k,1),1),cubes%xs(cubes%elems(k,1),2),cubes%xs(cubes%elems(k,1),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+1
            end if
            if (isosurface(cubes%xs(cubes%elems(k,2),1),cubes%xs(cubes%elems(k,2),2),cubes%xs(cubes%elems(k,2),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+2
            end if
            if (isosurface(cubes%xs(cubes%elems(k,3),1),cubes%xs(cubes%elems(k,3),2),cubes%xs(cubes%elems(k,3),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+4
            end if
            if (isosurface(cubes%xs(cubes%elems(k,4),1),cubes%xs(cubes%elems(k,4),2),cubes%xs(cubes%elems(k,4),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+8
            end if
            if (isosurface(cubes%xs(cubes%elems(k,5),1),cubes%xs(cubes%elems(k,5),2),cubes%xs(cubes%elems(k,5),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+16
            end if
            if (isosurface(cubes%xs(cubes%elems(k,6),1),cubes%xs(cubes%elems(k,6),2),cubes%xs(cubes%elems(k,6),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+32
            end if
            if (isosurface(cubes%xs(cubes%elems(k,7),1),cubes%xs(cubes%elems(k,7),2),cubes%xs(cubes%elems(k,7),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+64
            end if
            if (isosurface(cubes%xs(cubes%elems(k,8),1),cubes%xs(cubes%elems(k,8),2),cubes%xs(cubes%elems(k,8),3)) < 0.0) then
                cubeindex(k)=cubeindex(k)+128
            end if
        end do
    end function find_cubeindex



    function get_table() result(tables)
        implicit NONE
        type(mctable) tables

        tables%Edge = (/'000000000000',&
        '000100001001',&
        '001000000011',&
        '001100001010',&
        '010000000110',&
        '010100001111',&
        '011000000101',&
        '011100001100',&
        '100000001100',&
        '100100000101',&
        '101000001111',&
        '101100000110',&
        '110000001010',&
        '110100000011',&
        '111000001001',&
        '111100000000',&
        '000110010000',&
        '000010011001',&
        '001110010011',&
        '001010011010',&
        '010110010110',&
        '010010011111',&
        '011110010101',&
        '011010011100',&
        '100110011100',&
        '100010010101',&
        '101110011111',&
        '101010010110',&
        '110110011010',&
        '110010010011',&
        '111110011001',&
        '111010010000',&
        '001000110000',&
        '001100111001',&
        '000000110011',&
        '000100111010',&
        '011000110110',&
        '011100111111',&
        '010000110101',&
        '010100111100',&
        '101000111100',&
        '101100110101',&
        '100000111111',&
        '100100110110',&
        '111000111010',&
        '111100110011',&
        '110000111001',&
        '110100110000',&
        '001110100000',&
        '001010101001',&
        '000110100011',&
        '000010101010',&
        '011110100110',&
        '011010101111',&
        '010110100101',&
        '010010101100',&
        '101110101100',&
        '101010100101',&
        '100110101111',&
        '100010100110',&
        '111110101010',&
        '111010100011',&
        '110110101001',&
        '110010100000',&
        '010001100000',&
        '010101101001',&
        '011001100011',&
        '011101101010',&
        '000001100110',&
        '000101101111',&
        '001001100101',&
        '001101101100',&
        '110001101100',&
        '110101100101',&
        '111001101111',&
        '111101100110',&
        '100001101010',&
        '100101100011',&
        '101001101001',&
        '101101100000',&
        '010111110000',&
        '010011111001',&
        '011111110011',&
        '011011111010',&
        '000111110110',&
        '000011111111',&
        '001111110101',&
        '001011111100',&
        '110111111100',&
        '110011110101',&
        '111111111111',&
        '111011110110',&
        '100111111010',&
        '100011110011',&
        '101111111001',&
        '101011110000',&
        '011001010000',&
        '011101011001',&
        '010001010011',&
        '010101011010',&
        '001001010110',&
        '001101011111',&
        '000001010101',&
        '000101011100',&
        '111001011100',&
        '111101010101',&
        '110001011111',&
        '110101010110',&
        '101001011010',&
        '101101010011',&
        '100001011001',&
        '100101010000',&
        '011111000000',&
        '011011001001',&
        '010111000011',&
        '010011001010',&
        '001111000110',&
        '001011001111',&
        '000111000101',&
        '000011001100',&
        '111111001100',&
        '111011000101',&
        '110111001111',&
        '110011000110',&
        '101111001010',&
        '101011000011',&
        '100111001001',&
        '100011000000',&
        '100011000000',&
        '100111001001',&
        '101011000011',&
        '101111001010',&
        '110011000110',&
        '110111001111',&
        '111011000101',&
        '111111001100',&
        '000011001100',&
        '000111000101',&
        '001011001111',&
        '001111000110',&
        '010011001010',&
        '010111000011',&
        '011011001001',&
        '011111000000',&
        '100101010000',&
        '100001011001',&
        '101101010011',&
        '101001011010',&
        '110101010110',&
        '110001011111',&
        '111101010101',&
        '111001011100',&
        '000101011100',&
        '000001010101',&
        '001101011111',&
        '001001010110',&
        '010101011010',&
        '010001010011',&
        '011101011001',&
        '011001010000',&
        '101011110000',&
        '101111111001',&
        '100011110011',&
        '100111111010',&
        '111011110110',&
        '111111111111',&
        '110011110101',&
        '110111111100',&
        '001011111100',&
        '001111110101',&
        '000011111111',&
        '000111110110',&
        '011011111010',&
        '011111110011',&
        '010011111001',&
        '010111110000',&
        '101101100000',&
        '101001101001',&
        '100101100011',&
        '100001101010',&
        '111101100110',&
        '111001101111',&
        '110101100101',&
        '110001101100',&
        '001101101100',&
        '001001100101',&
        '000101101111',&
        '000001100110',&
        '011101101010',&
        '011001100011',&
        '010101101001',&
        '010001100000',&
        '110010100000',&
        '110110101001',&
        '111010100011',&
        '111110101010',&
        '100010100110',&
        '100110101111',&
        '101010100101',&
        '101110101100',&
        '010010101100',&
        '010110100101',&
        '011010101111',&
        '011110100110',&
        '000010101010',&
        '000110100011',&
        '001010101001',&
        '001110100000',&
        '110100110000',&
        '110000111001',&
        '111100110011',&
        '111000111010',&
        '100100110110',&
        '100000111111',&
        '101100110101',&
        '101000111100',&
        '010100111100',&
        '010000110101',&
        '011100111111',&
        '011000110110',&
        '000100111010',&
        '000000110011',&
        '001100111001',&
        '001000110000',&
        '111010010000',&
        '111110011001',&
        '110010010011',&
        '110110011010',&
        '101010010110',&
        '101110011111',&
        '100010010101',&
        '100110011100',&
        '011010011100',&
        '011110010101',&
        '010010011111',&
        '010110010110',&
        '001010011010',&
        '001110010011',&
        '000010011001',&
        '000110010000',&
        '111100000000',&
        '111000001001',&
        '110100000011',&
        '110000001010',&
        '101100000110',&
        '101000001111',&
        '100100000101',&
        '100000001100',&
        '011100001100',&
        '011000000101',&
        '010100001111',&
        '010000000110',&
        '001100001010',&
        '001000000011',&
        '000100001001',&
        '000000000000'/)

        tables%Tris = transpose(reshape((/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        1,9,4,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        1,2,10,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        2,9,4,10,9,2,0,0,0,0,0,0,0,0,0,0,&
        2,3,11,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        1,9,4,2,3,11,0,0,0,0,0,0,0,0,0,0,&
        10,3,11,1,3,10,0,0,0,0,0,0,0,0,0,0,&
        3,9,4,3,11,9,11,10,9,0,0,0,0,0,0,0,&
        4,12,3,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        1,12,3,9,12,1,0,0,0,0,0,0,0,0,0,0,&
        2,10,1,3,4,12,0,0,0,0,0,0,0,0,0,0,&
        2,12,3,2,10,12,10,9,12,0,0,0,0,0,0,0,&
        4,11,2,12,11,4,0,0,0,0,0,0,0,0,0,0,&
        1,11,2,1,9,11,9,12,11,0,0,0,0,0,0,0,&
        4,10,1,4,12,10,12,11,10,0,0,0,0,0,0,0,&
        10,9,11,11,9,12,0,0,0,0,0,0,0,0,0,0,&
        5,8,9,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        5,4,1,8,4,5,0,0,0,0,0,0,0,0,0,0,&
        1,2,10,9,5,8,0,0,0,0,0,0,0,0,0,0,&
        5,2,10,5,8,2,8,4,2,0,0,0,0,0,0,0,&
        2,3,11,9,5,8,0,0,0,0,0,0,0,0,0,0,&
        4,5,8,4,1,5,2,3,11,0,0,0,0,0,0,0,&
        10,3,11,10,1,3,9,5,8,0,0,0,0,0,0,0,&
        3,11,10,3,10,8,3,8,4,8,10,5,0,0,0,0,&
        9,5,8,4,12,3,0,0,0,0,0,0,0,0,0,0,&
        12,5,8,12,3,5,3,1,5,0,0,0,0,0,0,0,&
        10,1,2,9,5,8,3,4,12,0,0,0,0,0,0,0,&
        5,8,12,10,5,12,10,12,3,10,3,2,0,0,0,0,&
        4,11,2,4,12,11,8,9,5,0,0,0,0,0,0,0,&
        2,12,11,2,5,12,2,1,5,8,12,5,0,0,0,0,&
        5,8,9,10,1,12,10,12,11,12,1,4,0,0,0,0,&
        5,8,12,5,12,10,10,12,11,0,0,0,0,0,0,0,&
        10,6,5,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        10,6,5,1,9,4,0,0,0,0,0,0,0,0,0,0,&
        1,6,5,2,6,1,0,0,0,0,0,0,0,0,0,0,&
        9,6,5,9,4,6,4,2,6,0,0,0,0,0,0,0,&
        2,3,11,10,6,5,0,0,0,0,0,0,0,0,0,0,&
        4,1,9,2,3,11,5,10,6,0,0,0,0,0,0,0,&
        6,3,11,6,5,3,5,1,3,0,0,0,0,0,0,0,&
        3,11,6,4,3,6,4,6,5,4,5,9,0,0,0,0,&
        10,6,5,3,4,12,0,0,0,0,0,0,0,0,0,0,&
        1,12,3,1,9,12,5,10,6,0,0,0,0,0,0,0,&
        1,6,5,1,2,6,3,4,12,0,0,0,0,0,0,0,&
        3,2,6,3,6,9,3,9,12,5,9,6,0,0,0,0,&
        11,4,12,11,2,4,10,6,5,0,0,0,0,0,0,0,&
        5,10,6,1,9,2,9,11,2,9,12,11,0,0,0,0,&
        6,5,1,6,1,12,6,12,11,12,1,4,0,0,0,0,&
        6,5,9,6,9,11,11,9,12,0,0,0,0,0,0,0,&
        10,8,9,6,8,10,0,0,0,0,0,0,0,0,0,0,&
        10,4,1,10,6,4,6,8,4,0,0,0,0,0,0,0,&
        1,8,9,1,2,8,2,6,8,0,0,0,0,0,0,0,&
        2,6,4,4,6,8,0,0,0,0,0,0,0,0,0,0,&
        10,8,9,10,6,8,11,2,3,0,0,0,0,0,0,0,&
        11,2,3,10,6,1,6,4,1,6,8,4,0,0,0,0,&
        9,1,3,9,3,6,9,6,8,11,6,3,0,0,0,0,&
        3,11,6,3,6,4,4,6,8,0,0,0,0,0,0,0,&
        8,10,6,8,9,10,4,12,3,0,0,0,0,0,0,0,&
        10,6,8,10,8,3,10,3,1,3,8,12,0,0,0,0,&
        3,4,12,1,2,9,2,8,9,2,6,8,0,0,0,0,&
        12,3,2,12,2,8,8,2,6,0,0,0,0,0,0,0,&
        10,6,9,9,6,8,11,2,4,11,4,12,0,0,0,0,&
        6,8,1,6,1,10,8,12,1,2,1,11,12,11,1,0,&
        12,11,1,12,1,4,11,6,1,9,1,8,6,8,1,0,&
        12,11,6,8,12,6,0,0,0,0,0,0,0,0,0,0,&
        11,7,6,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        1,9,4,6,11,7,0,0,0,0,0,0,0,0,0,0,&
        10,1,2,6,11,7,0,0,0,0,0,0,0,0,0,0,&
        2,9,4,2,10,9,6,11,7,0,0,0,0,0,0,0,&
        2,7,6,3,7,2,0,0,0,0,0,0,0,0,0,0,&
        2,7,6,2,3,7,4,1,9,0,0,0,0,0,0,0,&
        10,7,6,10,1,7,1,3,7,0,0,0,0,0,0,0,&
        6,10,9,6,9,3,6,3,7,4,3,9,0,0,0,0,&
        3,4,12,11,7,6,0,0,0,0,0,0,0,0,0,0,&
        12,1,9,12,3,1,11,7,6,0,0,0,0,0,0,0,&
        1,2,10,3,4,12,6,11,7,0,0,0,0,0,0,0,&
        6,11,7,2,10,3,10,12,3,10,9,12,0,0,0,0,&
        7,4,12,7,6,4,6,2,4,0,0,0,0,0,0,0,&
        1,9,12,1,12,6,1,6,2,6,12,7,0,0,0,0,&
        4,12,7,1,4,7,1,7,6,1,6,10,0,0,0,0,&
        7,6,10,7,10,12,12,10,9,0,0,0,0,0,0,0,&
        6,11,7,5,8,9,0,0,0,0,0,0,0,0,0,0,&
        5,4,1,5,8,4,7,6,11,0,0,0,0,0,0,0,&
        2,10,1,6,11,7,9,5,8,0,0,0,0,0,0,0,&
        11,7,6,2,10,8,2,8,4,8,10,5,0,0,0,0,&
        7,2,3,7,6,2,5,8,9,0,0,0,0,0,0,0,&
        2,3,6,6,3,7,4,1,5,4,5,8,0,0,0,0,&
        9,5,8,10,1,6,1,7,6,1,3,7,0,0,0,0,&
        8,4,10,8,10,5,4,3,10,6,10,7,3,7,10,0,&
        4,12,3,8,9,5,11,7,6,0,0,0,0,0,0,0,&
        6,11,7,5,8,3,5,3,1,3,8,12,0,0,0,0,&
        1,2,10,5,8,9,3,4,12,6,11,7,0,0,0,0,&
        10,3,2,10,12,3,10,5,12,8,12,5,6,11,7,0,&
        9,5,8,4,12,6,4,6,2,6,12,7,0,0,0,0,&
        6,2,12,6,12,7,2,1,12,8,12,5,1,5,12,0,&
        1,6,10,1,7,6,1,4,7,12,7,4,9,5,8,0,&
        7,6,10,7,10,12,5,8,10,8,12,10,0,0,0,0,&
        11,5,10,7,5,11,0,0,0,0,0,0,0,0,0,0,&
        5,11,7,5,10,11,1,9,4,0,0,0,0,0,0,0,&
        11,1,2,11,7,1,7,5,1,0,0,0,0,0,0,0,&
        9,4,2,9,2,7,9,7,5,7,2,11,0,0,0,0,&
        2,5,10,2,3,5,3,7,5,0,0,0,0,0,0,0,&
        4,1,9,2,3,10,3,5,10,3,7,5,0,0,0,0,&
        1,3,5,5,3,7,0,0,0,0,0,0,0,0,0,0,&
        9,4,3,9,3,5,5,3,7,0,0,0,0,0,0,0,&
        11,5,10,11,7,5,12,3,4,0,0,0,0,0,0,0,&
        1,9,3,3,9,12,5,10,11,5,11,7,0,0,0,0,&
        4,12,3,1,2,7,1,7,5,7,2,11,0,0,0,0,&
        7,5,2,7,2,11,5,9,2,3,2,12,9,12,2,0,&
        10,7,5,10,4,7,10,2,4,12,7,4,0,0,0,0,&
        9,12,2,9,2,1,12,7,2,10,2,5,7,5,2,0,&
        4,12,7,4,7,1,1,7,5,0,0,0,0,0,0,0,&
        7,5,9,12,7,9,0,0,0,0,0,0,0,0,0,0,&
        8,11,7,8,9,11,9,10,11,0,0,0,0,0,0,0,&
        1,8,4,1,11,8,1,10,11,7,8,11,0,0,0,0,&
        11,7,8,2,11,8,2,8,9,2,9,1,0,0,0,0,&
        11,7,8,11,8,2,2,8,4,0,0,0,0,0,0,0,&
        2,3,7,2,7,9,2,9,10,9,7,8,0,0,0,0,&
        3,7,10,3,10,2,7,8,10,1,10,4,8,4,10,0,&
        8,9,1,8,1,7,7,1,3,0,0,0,0,0,0,0,&
        8,4,3,7,8,3,0,0,0,0,0,0,0,0,0,0,&
        3,4,12,11,7,9,11,9,10,9,7,8,0,0,0,0,&
        3,1,8,3,8,12,1,10,8,7,8,11,10,11,8,0,&
        2,9,1,2,8,9,2,11,8,7,8,11,3,4,12,0,&
        12,3,2,12,2,8,11,7,2,7,8,2,0,0,0,0,&
        9,10,7,9,7,8,10,2,7,12,7,4,2,4,7,0,&
        1,10,2,12,7,8,0,0,0,0,0,0,0,0,0,0,&
        8,9,1,8,1,7,4,12,1,12,7,1,0,0,0,0,&
        8,12,7,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        8,7,12,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        4,1,9,12,8,7,0,0,0,0,0,0,0,0,0,0,&
        1,2,10,12,8,7,0,0,0,0,0,0,0,0,0,0,&
        9,2,10,9,4,2,12,8,7,0,0,0,0,0,0,0,&
        11,2,3,7,12,8,0,0,0,0,0,0,0,0,0,0,&
        2,3,11,4,1,9,7,12,8,0,0,0,0,0,0,0,&
        3,10,1,3,11,10,7,12,8,0,0,0,0,0,0,0,&
        7,12,8,3,11,4,11,9,4,11,10,9,0,0,0,0,&
        8,3,4,7,3,8,0,0,0,0,0,0,0,0,0,0,&
        8,1,9,8,7,1,7,3,1,0,0,0,0,0,0,0,&
        3,8,7,3,4,8,1,2,10,0,0,0,0,0,0,0,&
        2,7,3,2,9,7,2,10,9,9,8,7,0,0,0,0,&
        11,8,7,11,2,8,2,4,8,0,0,0,0,0,0,0,&
        11,8,7,2,8,11,2,9,8,2,1,9,0,0,0,0,&
        1,4,8,1,8,11,1,11,10,7,11,8,0,0,0,0,&
        8,7,11,8,11,9,9,11,10,0,0,0,0,0,0,0,&
        7,9,5,12,9,7,0,0,0,0,0,0,0,0,0,0,&
        4,7,12,4,1,7,1,5,7,0,0,0,0,0,0,0,&
        9,7,12,9,5,7,10,1,2,0,0,0,0,0,0,0,&
        10,5,7,10,7,4,10,4,2,12,4,7,0,0,0,0,&
        7,9,5,7,12,9,3,11,2,0,0,0,0,0,0,0,&
        2,3,11,4,1,12,1,7,12,1,5,7,0,0,0,0,&
        5,12,9,5,7,12,1,3,10,3,11,10,0,0,0,0,&
        11,10,4,11,4,3,10,5,4,12,4,7,5,7,4,0,&
        9,3,4,9,5,3,5,7,3,0,0,0,0,0,0,0,&
        1,5,3,5,7,3,0,0,0,0,0,0,0,0,0,0,&
        2,10,1,3,4,5,3,5,7,5,4,9,0,0,0,0,&
        2,10,5,2,5,3,3,5,7,0,0,0,0,0,0,0,&
        9,2,4,9,7,2,9,5,7,7,11,2,0,0,0,0,&
        11,2,1,11,1,7,7,1,5,0,0,0,0,0,0,0,&
        5,7,4,5,4,9,7,11,4,1,4,10,11,10,4,0,&
        11,10,5,7,11,5,0,0,0,0,0,0,0,0,0,0,&
        5,10,6,8,7,12,0,0,0,0,0,0,0,0,0,0,&
        1,9,4,5,10,6,12,8,7,0,0,0,0,0,0,0,&
        6,1,2,6,5,1,8,7,12,0,0,0,0,0,0,0,&
        12,8,7,9,4,5,4,6,5,4,2,6,0,0,0,0,&
        10,6,5,11,2,3,8,7,12,0,0,0,0,0,0,0,&
        7,12,8,2,3,11,1,9,4,5,10,6,0,0,0,0,&
        8,7,12,6,5,11,5,3,11,5,1,3,0,0,0,0,&
        4,5,9,4,6,5,4,3,6,11,6,3,12,8,7,0,&
        8,3,4,8,7,3,6,5,10,0,0,0,0,0,0,0,&
        10,6,5,1,9,7,1,7,3,7,9,8,0,0,0,0,&
        4,7,3,4,8,7,2,6,1,6,5,1,0,0,0,0,&
        7,3,9,7,9,8,3,2,9,5,9,6,2,6,9,0,&
        10,6,5,11,2,7,2,8,7,2,4,8,0,0,0,0,&
        2,7,11,2,8,7,2,1,8,9,8,1,10,6,5,0,&
        5,1,11,5,11,6,1,4,11,7,11,8,4,8,11,0,&
        8,7,11,8,11,9,6,5,11,5,9,11,0,0,0,0,&
        7,10,6,7,12,10,12,9,10,0,0,0,0,0,0,0,&
        4,7,12,1,7,4,1,6,7,1,10,6,0,0,0,0,&
        1,12,9,1,6,12,1,2,6,6,7,12,0,0,0,0,&
        7,12,4,7,4,6,6,4,2,0,0,0,0,0,0,0,&
        2,3,11,10,6,12,10,12,9,12,6,7,0,0,0,0,&
        1,12,4,1,7,12,1,10,7,6,7,10,2,3,11,0,&
        12,9,6,12,6,7,9,1,6,11,6,3,1,3,6,0,&
        7,12,4,7,4,6,3,11,4,11,6,4,0,0,0,0,&
        6,9,10,6,3,9,6,7,3,4,9,3,0,0,0,0,&
        10,6,7,10,7,1,1,7,3,0,0,0,0,0,0,0,&
        2,6,9,2,9,1,6,7,9,4,9,3,7,3,9,0,&
        2,6,7,3,2,7,0,0,0,0,0,0,0,0,0,0,&
        2,4,7,2,7,11,4,9,7,6,7,10,9,10,7,0,&
        11,2,1,11,1,7,10,6,1,6,7,1,0,0,0,0,&
        1,4,9,6,7,11,0,0,0,0,0,0,0,0,0,0,&
        11,6,7,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        12,6,11,8,6,12,0,0,0,0,0,0,0,0,0,0,&
        12,6,11,12,8,6,9,4,1,0,0,0,0,0,0,0,&
        6,12,8,6,11,12,2,10,1,0,0,0,0,0,0,0,&
        11,8,6,11,12,8,10,9,2,9,4,2,0,0,0,0,&
        12,2,3,12,8,2,8,6,2,0,0,0,0,0,0,0,&
        1,9,4,2,3,8,2,8,6,8,3,12,0,0,0,0,&
        10,8,6,10,3,8,10,1,3,3,12,8,0,0,0,0,&
        8,6,3,8,3,12,6,10,3,4,3,9,10,9,3,0,&
        3,6,11,3,4,6,4,8,6,0,0,0,0,0,0,0,&
        9,3,1,9,6,3,9,8,6,11,3,6,0,0,0,0,&
        10,1,2,6,11,4,6,4,8,4,11,3,0,0,0,0,&
        10,9,3,10,3,2,9,8,3,11,3,6,8,6,3,0,&
        2,4,6,4,8,6,0,0,0,0,0,0,0,0,0,0,&
        1,9,8,1,8,2,2,8,6,0,0,0,0,0,0,0,&
        10,1,4,10,4,6,6,4,8,0,0,0,0,0,0,0,&
        10,9,8,6,10,8,0,0,0,0,0,0,0,0,0,0,&
        6,9,5,6,11,9,11,12,9,0,0,0,0,0,0,0,&
        6,1,5,6,12,1,6,11,12,12,4,1,0,0,0,0,&
        1,2,10,9,5,11,9,11,12,11,5,6,0,0,0,0,&
        11,12,5,11,5,6,12,4,5,10,5,2,4,2,5,0,&
        3,6,2,3,9,6,3,12,9,5,6,9,0,0,0,0,&
        1,5,12,1,12,4,5,6,12,3,12,2,6,2,12,0,&
        1,3,6,1,6,10,3,12,6,5,6,9,12,9,6,0,&
        10,5,6,3,12,4,0,0,0,0,0,0,0,0,0,0,&
        3,6,11,4,6,3,4,5,6,4,9,5,0,0,0,0,&
        6,11,3,6,3,5,5,3,1,0,0,0,0,0,0,0,&
        4,11,3,4,6,11,4,9,6,5,6,9,1,2,10,0,&
        6,11,3,6,3,5,2,10,3,10,5,3,0,0,0,0,&
        9,5,6,9,6,4,4,6,2,0,0,0,0,0,0,0,&
        1,5,6,2,1,6,0,0,0,0,0,0,0,0,0,0,&
        9,5,6,9,6,4,10,1,6,1,4,6,0,0,0,0,&
        10,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        5,12,8,5,10,12,10,11,12,0,0,0,0,0,0,0,&
        1,9,4,5,10,8,10,12,8,10,11,12,0,0,0,0,&
        2,11,12,2,12,5,2,5,1,8,5,12,0,0,0,0,&
        4,2,5,4,5,9,2,11,5,8,5,12,11,12,5,0,&
        5,12,8,10,12,5,10,3,12,10,2,3,0,0,0,0,&
        10,8,5,10,12,8,10,2,12,3,12,2,1,9,4,0,&
        12,8,5,12,5,3,3,5,1,0,0,0,0,0,0,0,&
        12,8,5,12,5,3,9,4,5,4,3,5,0,0,0,0,&
        3,10,11,3,8,10,3,4,8,8,5,10,0,0,0,0,&
        10,11,8,10,8,5,11,3,8,9,8,1,3,1,8,0,&
        4,8,11,4,11,3,8,5,11,2,11,1,5,1,11,0,&
        2,11,3,9,8,5,0,0,0,0,0,0,0,0,0,0,&
        5,10,2,5,2,8,8,2,4,0,0,0,0,0,0,0,&
        5,10,2,5,2,8,1,9,2,9,8,2,0,0,0,0,&
        5,1,4,8,5,4,0,0,0,0,0,0,0,0,0,0,&
        5,9,8,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        10,11,9,11,12,9,0,0,0,0,0,0,0,0,0,0,&
        4,1,10,4,10,12,12,10,11,0,0,0,0,0,0,0,&
        1,2,11,1,11,9,9,11,12,0,0,0,0,0,0,0,&
        4,2,11,12,4,11,0,0,0,0,0,0,0,0,0,0,&
        2,3,12,2,12,10,10,12,9,0,0,0,0,0,0,0,&
        4,1,10,4,10,12,2,3,10,3,12,10,0,0,0,0,&
        1,3,12,9,1,12,0,0,0,0,0,0,0,0,0,0,&
        4,3,12,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        3,4,9,3,9,11,11,9,10,0,0,0,0,0,0,0,&
        10,11,3,1,10,3,0,0,0,0,0,0,0,0,0,0,&
        3,4,9,3,9,11,1,2,9,2,11,9,0,0,0,0,&
        2,11,3,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        2,4,9,10,2,9,0,0,0,0,0,0,0,0,0,0,&
        1,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        1,4,9,0,0,0,0,0,0,0,0,0,0,0,0,0,&
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/), (/16, 256/)))

    end function get_table

end module marching_cubes