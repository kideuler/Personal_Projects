module utilities
    use init_mesh
    implicit none
    contains

    subroutine Convert2obj(Surf)
        implicit NONE
        type(mesh),intent(in) :: Surf
        integer :: nelems,nv,k

        nelems = size(Surf%elems,1)
        nv = size(Surf%xs,1)
        
        open(1,file = 'test.obj', status = 'unknown')
        do k = 1,nv
            write(1,*) 'v ', Surf%xs(k,:)
        end do
        write(1,*) ''
        do k = 1,nelems
            write(1,*) 'f ', Surf%elems(k,:)
        end do

    end subroutine Convert2obj


    subroutine print_array(A)
        implicit NONE
        real, intent(in),dimension(:,:) :: A
        integer :: m,n,i
        m = size(A,1)
        n = size(A,2)
        do i = 1,m
            print*, A(i,:)
        end do
    end subroutine print_array

    subroutine print_array_int(A)
        implicit NONE
        integer, intent(in),dimension(:,:) :: A
        integer :: m,n,i
        m = size(A,1)
        n = size(A,2)
        do i = 1,m
            print*, A(i,:)
        end do
    end subroutine print_array_int

    subroutine angles(mesh1)
        implicit NONE
        type(mesh),intent(in) :: mesh1
        integer :: nelems,k
        real,allocatable :: ang(:,:)
        real :: mx(3),mn(3)
        REAL, PARAMETER :: Pi = 3.1415927


        nelems = size(mesh1%elems,1)
        allocate(ang(nelems,3))
        do k = 1, nelems
            ang(k,1) = 180*(1/pi)*Acos(dot_product(mesh1%xs(mesh1%elems(k,2),:)-mesh1%xs(mesh1%elems(k,1),:), &
            mesh1%xs(mesh1%elems(k,3),:)-mesh1%xs(mesh1%elems(k,1),:))/(norm2(mesh1%xs(mesh1%elems(k,2),:)-&
            mesh1%xs(mesh1%elems(k,1),:))*&
            norm2(mesh1%xs(mesh1%elems(k,3),:)-mesh1%xs(mesh1%elems(k,1),:))))

            ang(k,2) = 180*(1/pi)*Acos(dot_product(mesh1%xs(mesh1%elems(k,1),:)-mesh1%xs(mesh1%elems(k,2),:), &
            mesh1%xs(mesh1%elems(k,3),:)-mesh1%xs(mesh1%elems(k,2),:))/(norm2(mesh1%xs(mesh1%elems(k,1),:)-&
            mesh1%xs(mesh1%elems(k,2),:))*&
            norm2(mesh1%xs(mesh1%elems(k,3),:)-mesh1%xs(mesh1%elems(k,2),:))))

            ang(k,3) = 180*(1/pi)*Acos(dot_product(mesh1%xs(mesh1%elems(k,1),:)-mesh1%xs(mesh1%elems(k,3),:), &
            mesh1%xs(mesh1%elems(k,2),:)-mesh1%xs(mesh1%elems(k,3),:))/(norm2(mesh1%xs(mesh1%elems(k,1),:)-&
            mesh1%xs(mesh1%elems(k,3),:))*&
            norm2(mesh1%xs(mesh1%elems(k,2),:)-mesh1%xs(mesh1%elems(k,3),:))))
        end do

        mn(1) = minval(ang(:,1))
        mn(2) = minval(ang(:,2))
        mn(3) = minval(ang(:,3))

        mx(1) = maxval(ang(:,1))
        mx(2) = maxval(ang(:,2))
        mx(3) = maxval(ang(:,3))

        print *, 'min = ', mn
        print *, 'max = ', mx
        deallocate(ang)
    end subroutine angles

end module utilities