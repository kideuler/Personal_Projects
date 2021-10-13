program main
    use init_mesh
    use marching_cubes
    use utilities
    implicit none

    type(mesh) :: cubes,surf

    ! generates the initial cube mesh
    cubes = gen_cubes(20)

    ! from the cube mesh generates the surface mesh
    surf = March_cubes(cubes)

    call Convert2obj(surf)
end program main