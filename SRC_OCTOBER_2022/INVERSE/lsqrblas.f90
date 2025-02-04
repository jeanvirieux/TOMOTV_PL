module s_lsqr_blas
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     File lsqrblas.f90   (double precision)  bascule en simple precision
!
!     This file contains the following BLAS routines
!        dcopy, ddot, dnrm2, dscal
!     required by subroutines LSQR and Acheck.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SCOPY copies a vector X to a vector Y.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!    The routine uses unrolled loops for increments equal to one.
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of elements in DX and DY.
!
!    Input, real ( kind = 4 ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries of DX.
!
!    Output, real ( kind = 4 ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries of DY.


      subroutine  scopy(n,sx,incx,sy,incy)

      implicit none
      real sx(*),sy(*)
      integer i,incx,incy,ix,iy,m,n

      if ( n <= 0 ) then
         return
      end if

      if ( incx == 1 .and. incy == 1 ) then

         m = mod ( n, 7 )

         if ( m /= 0 ) then
            sy(1:m) = sx(1:m)
         end if

         do i = m+1, n, 7
            sy(i) = sx(i)
            sy(i + 1) = sx(i + 1)
            sy(i + 2) = sx(i + 2)
            sy(i + 3) = sx(i + 3)
            sy(i + 4) = sx(i + 4)
            sy(i + 5) = sx(i + 5)
            sy(i + 6) = sx(i + 6)
         end do

        else

           if ( 0 <= incx ) then
              ix = 1
           else
              ix = ( -n + 1 ) * incx + 1
           end if

           if ( 0 <= incy ) then
              iy = 1
           else
              iy = ( -n + 1 ) * incy + 1
           end if

           do i = 1, n
              sy(iy) = sx(ix)
              ix = ix + incx
              iy = iy + incy
           end do
        end if
        return
end subroutine scopy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! SDOT forms the dot product of two vectors.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!    This routine uses unrolled loops for increments equal to one.
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = 4 ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries in DX.
!
!    Input, real ( kind = 4 ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries in DY.
!
!    Output, real ( kind = 4 ) DDOT, the sum of the product of the 
!    corresponding entries of DX and DY.


      real function sdot(n,sx,incx,sy,incy)

      implicit         none
      real sx(*),sy(*),stemp
      integer          i,incx,incy,ix,iy,m,n

      sdot = 0.0
      stemp = 0.0
      if ( n <= 0 ) then
         return
      end if

!  Code for unequal increments or equal increments
!  not equal to 1.

      if ( incx /= 1 .or. incy /= 1 ) then

         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         if ( 0 <= incy ) then
            iy = 1
         else
            iy = ( - n + 1 ) * incy + 1
         end if

         do i = 1, n
            stemp = stemp + sx(ix) * sy(iy)
            ix = ix + incx
            iy = iy + incy
         end do

!  Code for both increments equal to 1.

        else

           m = mod ( n, 5 )

           do i = 1, m
              stemp = stemp + sx(i) * sy(i)
           end do

           do i = m+1, n, 5
              stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) + sx(i+2)*sy(i+2) &
                                          + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
           end do

        end if

        sdot = stemp
        return
end function sdot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************80
!
!! SNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!    This routine uses single precision real arithmetic.
!     SNRM2 ( X ) = sqrt ( X' * X )
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Sven Hammarling
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 4 ) DNRM2, the Euclidean norm of X.
!

      real function snrm2 ( n, x, incx)
      implicit         none
      integer          ix,n,incx
      real x(*), ssq,absxi,norm,scale

      if ( n < 1 .or. incx < 1 ) then
         norm  = 0.
      else if ( n == 1 ) then
         norm  = abs ( x(1) )
      else
         scale = 0.
         ssq = 1.

         do ix = 1, 1 + ( n - 1 )*incx, incx
            if ( x(ix) /= 0. ) then
               absxi = abs ( x(ix) )
               if ( scale < absxi ) then
                  ssq = 1. + ssq * ( scale / absxi )**2
                  scale = absxi
               else
                  ssq = ssq + ( absxi / scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt ( ssq )
      end if

      snrm2 = norm
      return
end function snrm2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SSCAL scales a vector by a constant.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!
!  Modified:
!    08 April 1999
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) SA, the multiplier.
!
!    Input/output, real ( kind = 4 ) X(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive entries of X.
!

      subroutine  sscal(n,sa,x,incx)

      implicit none

      integer i
      integer incx
      integer ix
      integer m
      integer n
      real sa
      real x(*)

      if ( n <= 0 ) then
         return
      else if ( incx == 1 ) then
         m = mod ( n, 5 )
         x(1:m) = sa * x(1:m)

         do i = m+1, n, 5
            x(i)   = sa * x(i)
            x(i+1) = sa * x(i+1)
            x(i+2) = sa * x(i+2)
            x(i+3) = sa * x(i+3)
            x(i+4) = sa * x(i+4)
         end do
      else
         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         do i = 1, n
            x(ix) = sa * x(ix)
            ix = ix + incx
         end do

      end if

      return
end subroutine sscal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module s_lsqr_blas
