program main
    use geom
    use delta
    use fftw3
    use omp_lib
    implicit none
    integer :: i,j,k,n
    integer :: ii, jj, kk
    integer :: Ntri
    real(kind=4) :: nhx,nhy,nhz,v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z
    real(kind=4) :: e1, e2, e3
    real(kind=4) :: Volume = 0.0
    real(kind=4) :: maxE = -9999.0
    real(kind=4) :: minE = 9999.0
    real(kind=8) :: avgE = 0.0
    character(len=50) :: filename 

    integer :: cntE = 0
    character(len=80) :: header
    integer(kind=2) :: attribute_byte_count

    ! Triangulated geometry arrays
    real, allocatable, dimension(:) :: Atri, xl, yl, zl
    real, allocatable, dimension(:,:) :: nhat_F

    ! Eulerian field
    integer :: Nx, Ny, Nz
    real :: Lx, Ly, Lz
    real :: dx, dy, dz
    real, allocatable, dimension(:) :: xc, xm, yc, ym, zc, zm

    ! Fourier metrics for fast-Poisson solver
    real, allocatable, dimension(:) :: lmb_x_on_dx2, lmb_y_on_dy2, lmb_z_on_dz2
    real :: PI = 2.d0*dasin(1.d0) 
    type(C_PTR) :: fftw_plan_fwd
    type(C_PTR) :: fftw_plan_bwd
    type(C_PTR) :: ptr1, ptr2, ptr3, ptr4
    complex(C_DOUBLE_COMPLEX), pointer :: rhs_hat(:,:,:), phihat(:,:,:)
    real(C_DOUBLE), pointer :: phi(:,:,:), rhs(:,:,:) ! Poisson equation

    ! For Poisson reconstruction
    real, allocatable, dimension(:,:,:) :: Gx, Gy, Gz
    real :: weight, drx, dry, drz, rx, ry, rz
    integer :: nsup = 2
    real :: DI = -1.0
    real :: phimin, phimax, inv_dphi

    ! For OMP acceleration
    integer :: void, nthreads


    void =  fftw_init_threads()


    !$omp parallel
    !$omp master
    nthreads = OMP_GET_NUM_THREADS() 
    !$omp end master
    !$omp end parallel

    call  fftw_plan_with_nthreads(nthreads)

    ! filename = "stl/unitSphere_N3.stl"
    ! Nx = 64
    ! Ny = 64
    ! Nz = 64
    ! Lx = 10.0
    ! Ly = 10.0
    ! Lz = 10.0

    ! filename = "stl/unitSphere_N4.stl"
    ! Nx = 128
    ! Ny = 128
    ! Nz = 128
    ! Lx = 10.0
    ! Ly = 10.0
    ! Lz = 10.0

    ! filename = "stl/unitSphere_N5.stl"
    ! Nx = 256
    ! Ny = 256
    ! Nz = 256
    ! Lx = 10.0
    ! Ly = 10.0
    ! Lz = 10.0

    ! filename = "stl/unitSphere_N6.stl"
    ! Nx = 512
    ! Ny = 512
    ! Nz = 512
    ! Lx = 10.0
    ! Ly = 10.0
    ! Lz = 10.0

    ! filename = "stl/ellipsoid_N32.stl"
    ! Nx = 128
    ! Ny = 128
    ! Nz = 128
    ! Lx = 8.0
    ! Ly = 8.0
    ! Lz = 8.0

    filename = "stl/256_ellipse.stl"
    Nx = 256
    Ny = 256
    Nz = 256
    Lx = 1.0
    Ly = 1.0
    Lz = 1.0

    ! ---------------------------- Generate Eulerian grid ---------------------------------
    ! Grid spacing
    dx = Lx / Nx
    dy = Ly / Ny
    dz = Lz / Nz


    ! Allocate memory for spatial locations

    allocate( xc( 1:Nx  )  )
    allocate( xm( 1:Nx  )  )

    allocate( yc( 1:Ny  )  )
    allocate( ym( 1:Ny  )  )

    allocate( zc( 1:Nz  )  )
    allocate( zm( 1:Nz  )  )

    do i = 1,Nx
        xc(i) = dble(i - 1) * dx
        xm(i) = ( dble(i-1) + 0.5 ) * dx
    enddo

    do i = 1,Ny
        yc(i) = dble(i - 1) * dy
        ym(i) = ( dble(i-1) + 0.5 ) * dy
    enddo

    do i = 1,Nz
        zc(i) = dble(i - 1) * dz
        zm(i) = ( dble(i-1) + 0.5 ) * dz
    enddo

    allocate( Gx(1:Nx, 1:Ny, 1:Nz) ) ; allocate( Gy(1:Nx, 1:Ny, 1:Nz) ) ; allocate( Gz(1:Nx, 1:Ny, 1:Nz) )

    ! ---------------------------- FFT solver memory ---------------------------------

    ! Modified wavenumbers
    allocate( lmb_x_on_dx2(Nx/2+1) ) ; allocate( lmb_y_on_dy2(Ny) ) ; allocate( lmb_z_on_dz2(Nz) )

    do i = 1,Nx/2+1
        lmb_x_on_dx2(i) = ( 2.0 * cos(2.0*PI*(real(i)-1.0) / Nx ) - 2.0 ) / dx**2
    enddo

    do j = 1,Ny
        lmb_y_on_dy2(j) = ( 2.0 * cos(2.0*PI*(real(j)-1.0) / Ny ) - 2.0 ) / dy**2
    enddo

    do k = 1,Nz
        lmb_z_on_dz2(k) = ( 2.0 * cos(2.0*PI*(real(k)-1.0) / Nz ) - 2.0 ) / dz**2
    enddo

        ! Solution array real(Nx,Nz)
    ptr1 = fftw_alloc_real(int(Nx * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr1, phi, [Nx,Ny,Nz])

    ptr2 = fftw_alloc_complex(int((Nx/2+1) * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr2, phihat, [Nx/2+1,Ny,Nz])

    ptr3 = fftw_alloc_complex(int((Nx/2+1) * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr3, rhs_hat, [Nx/2+1,Ny,Nz])

    ptr4 = fftw_alloc_real(int(Nx * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr4, rhs, [Nx,Ny,Nz])

    write(*,*) "Generating FFT plan..."

    ! 3D real to complex plan
    fftw_plan_fwd = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, rhs(:,:,:), rhs_hat(:,:,:), FFTW_ESTIMATE)

    ! 3D complex to real transform
    fftw_plan_bwd = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, phihat(:,:,:), phi(:,:,:), FFTW_ESTIMATE)


    !----------- First read the STL file ------------------------------------------------------------
    open(unit=10, file=filename, access="stream", status="old", form="unformatted")
  
    read(10) header             ! Read 80-byte header
    read(10) Ntri               ! Read number of triangles

    allocate( Atri(1:Ntri) ) ; allocate( nhat_F(3,1:Ntri) )
    allocate( xl(1:Ntri) ) ;  allocate( yl(1:Ntri) ) ;  allocate( zl(1:Ntri) )

    write(*,*) "Reading in .stl file with num_tri = ", Ntri

    do n = 1, Ntri
        read(10) nhx, nhy, nhz
        read(10) v1x, v1y, v1z           
        read(10) v2x, v2y, v2z            
        read(10) v3x, v3y, v3z             
        read(10) attribute_byte_count ! Attribute byte count, dummy/un-used

        ! Normalize the normal vector
        nhat_F(1,n) = dble (nhx) ; nhat_F(2,n) = dble (nhy) ; nhat_F(3,n) = dble (nhz)

        nhat_F(:,n) = nhat_F(:,n) / norm2( nhat_F(:,n) )

        ! Compute Lagrangian marker as triangle centroid
        xl(n) = dble( ( v1x + v2x + v3x ) / 3.0 )
        yl(n) = dble( ( v1y + v2y + v3y ) / 3.0 )
        zl(n) = dble( ( v1z + v2z + v3z ) / 3.0 )

        ! Compute triangle area
        Atri(n) = dble ( compute_area_tri(v1x,v1y,v1z, v2x,v2y,v2z, v3x,v3y,v3z)  )

        ! Edge lengths for diagnostics
        e1 = dist3d(v1x,v1y,v1z, v2x,v2y,v2z)
        e2 = dist3d(v1x,v1y,v1z, v3x,v3y,v3z)
        e3 = dist3d(v3x,v3y,v3z, v2x,v2y,v2z)

        maxE = max(maxE, e1, e2, e3)
        minE = min(minE, e1, e2, e3)

        avgE = avgE + dble( e1 + e2 + e3 )
        cntE = cntE + 3

        ! Volume of Lagrangian object: accumulate signed Tetrahedral volume
        Volume =  Volume    +  v1x * (v2y*v3z - v2z*v3y) &
                            +  v2x * (v3y*v1z - v3z*v1y) &
                            +  v3x * (v1y*v2z - v1z*v2y)
    end do
    close(10)

    Volume = Volume / 6.0 ! normalization
    avgE = avgE / dble(cntE)

    write(*,*) "Volume of triangulated surface: ", Volume
    write(*,*) "min, avg, max edge/dx: ", minE/dx, avgE/dx, maxE/dx 


    !--------------- Re-centre geometry if needed -------------------------
    do n = 1,Ntri
        xl(n) = xl(n) + 0.5*Lx
        yl(n) = yl(n) + 0.5*Ly
        zl(n) = zl(n) + 0.5*Lz
    enddo


    write(*,*) "Finished initialisation...beginning Lagrangian loop"

    ! Begin Lagrangian ----> Eulerian interpolation
    do n = 1,Ntri

        !-------------------------- X-gradient --------------------------------------------------------
        i = nint (xl(n) / dx ) + 1
        j = floor(yl(n) / dy ) + 1
        k = floor(zl(n) / dz ) + 1

        do kk = k-nsup,k+nsup

            rz = ( zl(n) - zm(kk) ) / dz
            drz = deltaFunc(rz)

            do jj = j-nsup,j+nsup

                ry = ( yl(n) - ym(jj) ) / dy
                dry = deltaFunc(ry)

                do ii = i-nsup,i+nsup

                    rx = ( xl(n) - xc(ii) ) / dx
                    drx = deltaFunc(rx)

                    weight = drx * dry * drz

                    Gx(ii,jj,kk) = Gx(ii,jj,kk) + DI * nhat_F(1,n)  * weight * Atri(n) / (dx*dy*dz)

                enddo
            enddo
        enddo


        !-------------------------- Y-gradient --------------------------------------------------------
        i = floor(xl(n) / dx ) + 1
        j = nint (yl(n) / dy ) + 1
        k = floor(zl(n) / dz ) + 1

        do kk = k-nsup,k+nsup

            rz = ( zl(n) - zm(kk) ) / dz
            drz = deltaFunc(rz)

            do jj = j-nsup,j+nsup

                ry = ( yl(n) - yc(jj) ) / dy
                dry = deltaFunc(ry)

                do ii = i-nsup,i+nsup

                    rx = ( xl(n) - xm(ii) ) / dx
                    drx = deltaFunc(rx)

                    weight = drx * dry * drz

                    Gy(ii,jj,kk) = Gy(ii,jj,kk) + DI * nhat_F(2,n)  * weight * Atri(n) / (dx*dy*dz)

                enddo
            enddo
        enddo

        !-------------------------- Z-gradient --------------------------------------------------------
        i = floor(xl(n) / dx ) + 1
        j = floor(yl(n) / dy ) + 1
        k = nint (zl(n) / dz ) + 1

        do kk = k-nsup,k+nsup

            rz = ( zl(n) - zc(kk) ) / dz
            drz = deltaFunc(rz)

            do jj = j-nsup,j+nsup

                ry = ( yl(n) - ym(jj) ) / dy
                dry = deltaFunc(ry)

                do ii = i-nsup,i+nsup

                    rx = ( xl(n) - xm(ii) ) / dx
                    drx = deltaFunc(rx)

                    weight = drx * dry * drz

                    Gz(ii,jj,kk) = Gz(ii,jj,kk) + DI * nhat_F(3,n)  * weight * Atri(n) / (dx*dy*dz)

                enddo
            enddo
        enddo

    enddo

    write(*,*) "Finished Lagrangian loop...beginning Poisson solve"


    ! Build RHS to the Poisson equation
    do k = 1,Nz-1
        do j = 1,Ny-1
            do i = 1,Nx-1
                rhs(i,j,k) = ( Gx(i+1,j,k) - Gx(i,j,k) ) / dx   + &
                             ( Gy(i,j+1,k) - Gy(i,j,k) ) / dy   + &
                             ( Gz(i,j,k+1) - Gz(i,j,k) ) / dz   
            enddo
        enddo
    enddo


    call dfftw_execute_dft_r2c(fftw_plan_fwd, rhs(:,:,:), rhs_hat(:,:,:))

    do k = 1,Nz
        do j = 1,Ny
            do i = 1,Nx/2+1

                rhs_hat(i,j,k) =  rhs_hat(i,j,k) / dble(Nx * Ny * Nz) ! Normalisation for back-transform
                phihat(i,j,k) = rhs_hat(i,j,k) / ( lmb_x_on_dx2(i) + lmb_y_on_dy2(j) + lmb_z_on_dz2(k) )

            enddo
        enddo
    enddo

    ! Arbitrary
    phihat(1,1,1) = dble(Volume) / (Lx * Ly * Lz)

    call dfftw_execute_dft_c2r(fftw_plan_bwd, phihat(:,:,:), phi(:,:,:))

    write(*,*) "Finished Poisson solve, normalising and writing data..."

    ! ! This is a non-volume-conserving operation
    ! !
    ! ! Re-normalisation to [0,1]
    ! phimin = minval(phi) ; phimax = maxval(phi) 
    ! inv_dphi = 1.0 / (phimax - phimin)
    ! do k = 1,Nz
    !     do j = 1,Ny
    !         do i = 1,Nx
    !             phi(i,j,k) = (phi(i,j,k) - phimin ) * inv_dphi
    !         enddo
    !     enddo
    ! enddo

    !call write3DField(phi,Nx,Ny,Nz,dx,dy,dz,'output')

    !     ! Re-normalisation to [0,1]
    ! phimin = minval(phi) ; phimax = maxval(phi) 
    ! inv_dphi = 0.0
    ! ! Compute and back out the mean
    ! do k = 1,Nz
    !     do j = 1,Ny
    !         do i = 1,Nx
    !             inv_dphi = inv_dphi + phi(i,j,k)*dx*dy*dz
    !         enddo
    !     enddo
    ! enddo


    ! write(*,*) "Volume is: ", inv_dphi
    ! write(*,*) "phimin, phimax is: ", phimin, phimax

    call write3DField(phi,Nx,Ny,Nz,dx,dy,dz,'output')


    deallocate(Atri) ; deallocate(nhat_F)
    deallocate(xl) ; deallocate(yl) ; deallocate(zl)

    deallocate(xc) ; deallocate(xm) ;  deallocate(yc) ; deallocate(ym) ;  deallocate(zc) ; deallocate(zm) 

    deallocate(lmb_x_on_dx2) ; deallocate(lmb_y_on_dy2) ; deallocate(lmb_z_on_dz2)
    call fftw_destroy_plan(fftw_plan_fwd)
    call fftw_destroy_plan(fftw_plan_bwd)
    call fftw_free(ptr1) ! phi
    call fftw_free(ptr2) ! phihat
    call fftw_free(ptr3) ! rhs_hat
    call fftw_free(ptr4) ! rhs
    call fftw_cleanup_threads

    deallocate(Gx) ; deallocate(Gy) ; deallocate(Gz)


end program main