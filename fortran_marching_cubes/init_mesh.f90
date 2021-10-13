MODULE init_mesh
    implicit NONE
    ! data structure for the mesh
    type mesh
        integer, allocatable :: elems(:,:)
        real, allocatable :: xs(:,:)
        integer :: npoints
    end type mesh
    contains
    
    function gen_cubes(npoints) result(cubes) ! generate the initial cube mesh
        implicit NONE

        ! declaring variables
        integer, intent(in) :: npoints
        type(mesh) :: cubes
        integer :: i,j,k,nelems,nv,kk,ii,jj,eid
        real :: bar(npoints)

        nelems = (npoints-1)**3 ! number of elements in connectivity table
        nv = npoints**3 ! total number of points in mesh

        allocate(cubes%elems(nelems,8))
        allocate(cubes%xs(nv,3))
        cubes%npoints = npoints

        ! initial set of equally spaced points [0, 1]
        do k = 1,npoints
            bar(k) = real(k-1)/real(npoints-1)
        end do

        ! forming the set of points in the cube mesh by taking the tensor product of the equally spaced points
        ! in bar
        k = 1
        do i = 1,npoints
            do j = 1,npoints
                do kk = 1,npoints
                    cubes%xs(k,:) = (/bar(i), bar(j), bar(kk)/)
                    k=k+1
                end do
            end do
        end do

        ! forming the connectivity table using CGNS format
        ! Geometry formatting information can be found here https://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid
        eid=1
        do k = 1,(npoints-1)
            do j = 1,(npoints-1)
                ii = (k-1)*npoints*npoints+(j-1)*npoints+1
                jj = ii + npoints
                do i = 1,(npoints-1)
                    cubes%elems(eid,:) =(/jj,jj+npoints**2,ii+npoints**2,ii, jj+1,jj+npoints**2+1,ii+npoints**2+1,ii+1/)
                    eid=eid+1
                    ii=ii+1
                    jj=jj+1
                end do
            end do
        end do
    end function gen_cubes
end module init_mesh