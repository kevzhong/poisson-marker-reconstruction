module delta
    implicit none
    contains

    function deltaFunc(r) result(rdelta)
        implicit none
        real :: r, rdelta

        rdelta = brDelta(r)
        !rdelta = peskinDelta(r)

    end function deltaFunc


    ! Brackbill & Ruppel (1986) delta function
    function brDelta(r) result(rdelta)
        implicit none
        real :: r, rdelta

        if( ( abs(r) .ge. 0.0 ) .and. ( abs(r) .le. 1.0)  ) then 
            rdelta = 2.0/3.0 - r**2 + 0.5 * abs(r)**3
        elseif ( (abs(r) .gt. 1.0) .and. (abs(r) .le. 2.0) ) then
            rdelta = ( 2.0 - abs(r))**3 / 6.0
        else
            rdelta = 0.0
        endif
    end function brDelta

    ! Peksin (1977) delta function
    function peskinDelta(r) result(rdelta)
        implicit none
        real :: r, rdelta
        real :: PI = 2.d0*dasin(1.d0) 

        if ( abs(r) .lt. 2.0 ) then
            rdelta = 0.25 * ( 1.0 + cos(0.5*PI*r) )
        else
            rdelta = 0.0
        endif
    end function peskinDelta

end module delta