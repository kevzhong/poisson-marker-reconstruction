module geom
    implicit none
    contains

    function dist3d(x1,y1,z1, x2,y2,z2) result(dist)
        implicit none
        real(kind=4) :: x1, y1, z1, x2, y2, z2, dist

        dist = sqrt ( (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2 )

    end function dist3d

    function compute_area_tri(x1,y1,z1, x2,y2,z2, x3,y3,z3) result(area)
        implicit none
        real(kind=4) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, area
        real(kind=4) :: a, b, c, s

        a = dist3d(x1,y1,z1,x2,y2,z2)
        b = dist3d(x2,y2,z2,x3,y3,z3)
        c = dist3d(x3,y3,z3,x1,y1,z1)


        s = (a + b + c) / 2   
        area = sqrt ( (s*(s-a) * (s-b)*(s-c)) )

    end function compute_area_tri


end module geom