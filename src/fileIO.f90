subroutine write3DField(array,Nx,Ny,Nz,dx,dy,dz,fn)
    ! Write array(Nx,Ny) field as binary data
    implicit none
    integer :: Nx, Ny, Nz
    real(8), intent(in) :: dx, dy, dz   ! Grid spacing: (dx, dy)
    character(len=50) :: filename 
    real, dimension(Nx,Ny,Nz) :: array
    character(len=*) :: fn

    write(filename, '(A,".dat")') trim(fn)

    open(unit=10, file=filename, status='replace', access='stream', action='write')
    write(10) array
    close(10)

    call write_3d_xdmf(fn, Nx, Ny, Nz, dx, dy, dz) ! Writing header for viewing in Paraview
end subroutine write3DField

! Write .xmf files to read binary dumps in Paraview
subroutine write_3d_xdmf(fpref,nx,ny,nz,dx,dy,dz)
    implicit none
    character(len=50) :: filename  ! xdmf file name
    character(len=50) :: bfilename ! binary filename
    character(len=*) :: fpref
    integer, intent(in)   :: nx, ny, nz   
    real(8), intent(in) :: dx, dy, dz

    character(len=32) :: nx_str, ny_str, nz_str
    integer :: unit       

    write(nx_str, '(I0)') Nx
    write(ny_str, '(I0)') Ny
    write(nz_str, '(I0)') Nz


    write(filename, '(A,".xmf")') trim(fpref)
    write(bfilename, '(A,".dat")') trim(fpref)


    ! Open the file
    open(newunit=unit, file=filename, status='replace', action='write')
    ! Write the XDMF header
    write(unit, '(A)') '<?xml version="1.0" ?>'
    write(unit, '(A)') '<Xdmf Version="3.0">'
    write(unit, '(A)') '  <Domain>'
    write(unit, '(A)') '    <Grid Name="StructuredGrid" GridType="Uniform">'
    write(unit, '(A)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
                       trim(adjustl(nz_str)) // ' ' // &
                       trim(adjustl(ny_str)) // ' ' // &
                       trim(adjustl(nx_str)) // '"/>'
    write(unit, '(A)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
    write(unit, '(A)') '      <DataItem Format="XML" Dimensions="3">0.0 0.0 0.0</DataItem>'

    write(unit, '(A,F12.6, A, F12.6, A, F12.6, A)') '        <DataItem Format="XML" Dimensions="3">',&
                 dx, ' ', dy, ' ', dz, '</DataItem>'

    write(unit, '(A)') '      </Geometry>'
    write(unit, '(A)') '      <Attribute Name="' // fpref // '" AttributeType="Scalar" Center="Node">'
    write(unit, '(A)') '        <DataItem Format="Binary" Endian="Little" DataType="Float" Precision="8" Dimensions="' // &
                                trim(adjustl(nz_str)) // ' ' // &
                                trim(adjustl(ny_str)) // ' ' // &
                                trim(adjustl(nx_str)) // '">' // trim(adjustl(bfilename)) // '</DataItem>'

    write(unit, '(A)') '      </Attribute>'
    write(unit, '(A)') '    </Grid>'
    write(unit, '(A)') '  </Domain>'
    write(unit, '(A)') '</Xdmf>'

    ! Close the file
    close(unit)

end subroutine write_3d_xdmf