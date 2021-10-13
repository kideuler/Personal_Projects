module surface
    implicit NONE
    contains

    function isosurface(x,y,z) result(S)
        implicit NONE
        real, intent(in) :: x,y,z
        real :: S
        S = (sqrt((x-0.5)**2 + (y-0.5)**2) - 0.3)**2 + (z-0.5)**2 - 0.01;
    end function isosurface

    function isosurface_grad(xs) result(grad)
        implicit NONE
        real :: grad(3),x,y,z
        real,intent(in) :: xs(3)

        x=xs(1)
        y=xs(2)
        z=xs(3)
        
    end function isosurface_grad

end module surface