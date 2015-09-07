      subroutine lagint(y, dy, x, x2, x3, xo, m, n, l)
      integer m, n, l, i, j, k
      real*8 x2(m,n), x3(m,n), xo(l), y(l), dy(l)
      real*8 xi(m*n), yi(m*n)
      integer x, n0, nl, nr
      real*8 xin(x), yin(x), xi1,xi2,fi1,fi2

Cf2py intent(in) x, x2, x3, xo
Cf2py intent(out) y(l), dy(l)

C x - number of points = order + 1
C order of the interpolant must be <= max(m,n)
      if(x>max(m,n)) then
         x = max(m,n)
      endif
      
      do k=1,l

      do i=1,m*n
        if(m.eq.1) then
           xi(i) = x2(1,i)
           yi(i) = x3(1,i)
         else
           xi(i) = x2(i,1)
           yi(i) = x3(i,1)
         endif
      enddo
C     

C nearest nod in the grid right to xo 
      do i=1,m*n
        if(xo(k).le.xi(i)) then       
           n0 = i;
           exit
         endif
      enddo
C number of elements left to xo
      nl = n0 - 1
C number of elements right to xo
      nr = m*n - n0 + 1
C      nr = m*n - n0
C cut x points around xo
      if(floor(dble(x)/2)>nl) then
        xin = xi(1:x)
        yin = yi(1:x)
      endif
      if(ceiling(dble(x)/2)>nr) then
        xin = xi(m*n-x+1:m*n)
        yin = yi(m*n-x+1:m*n)
      endif
      if((floor(dble(x)/2)<=nl).and.(ceiling(dble(x)/2)<=nr)) then
        xin = xi(n0-floor(dble(x)/2):n0+ceiling(dble(x)/2)-1)
        yin = yi(n0-floor(dble(x)/2):n0+ceiling(dble(x)/2)-1)
      endif
      
C Subroutine performing the Lagrange interpolation with the
C Aitken method. y: interpolated value. dy: error estimated.

      do  i = 1, x
        do  j = 1, x-i
          xi1 = xin(j)
          xi2 = xin(j+i)
          fi1 = yin(j)
          fi2 = yin(j+1)
          yin(j) = (xo(k)-xi1)/(xi2-xi1)*fi2+(xo(k)-xi2)/(xi1-xi2)*fi1
        enddo
      enddo
      

      y(k) = yin(1) 
      dy(k) = (abs(y(k)-fi1)+abs(y(k)-fi2))/2.d0

      enddo

      return
      end


      subroutine lagintt(y, dy, x, x2, x3, xo, m, n, l)
      integer m, n, l, i, j, k
      real*8 x2(m,n), x3(m,n), xo(l), y(l), dy(l)
      real*8 xi(m*n), yi(m*n)
      integer x, n0, nl, nr
      real*8 xin(x), yin(x), xi1,xi2,fi1,fi2

Cf2py intent(in) x, x2, x3, xo
Cf2py intent(out) y(l), dy(l)

C x - number of points = order + 1
C order of the interpolant must be <= max(m,n)
      if(x>max(m,n)) then
         x = max(m,n)
      endif
      
      do k=1,l

      do i=1,m*n
        if(m.eq.1) then
           xi(i) = x2(1,i)
           yi(i) = x3(1,i)
         else
           xi(i) = x2(i,1)
           yi(i) = x3(i,1)
         endif
      enddo
C     

C nearest nod in the grid right to xo 
      do i=1,m*n
        if(xo(k).le.xi(i)) then       
           n0 = i;
           exit
         endif
      enddo

C number of points to cut from the left-hand side
      nl = floor(dble(x))/2
C number of points to cut from the right-hand side
      nr = ceiling(dble(x))/2
C check/correct bounds:
      if (size(xi(1:n0-1)) < nl) then
          nr = nr + nl - size(xi(1:n0-1))
          nl = size(xi(1:n0-1))
      endif
      if (size(xi(n0:m*n)) < nr) then
          nl = nl + nr - size(xi(n0:m*n))
          nr = size(xi(n0:m*n))
      endif
            
C cut the proper piece:
      xin = xi(n0-nl:n0+nr-1)
      yin = yi(n0-nl:n0+nr-1)

      
C Subroutine performing the Lagrange interpolation with the
C Aitken method. y: interpolated value. dy: error estimated.

      do  i = 1, x
        do  j = 1, x-i
          xi1 = xin(j)
          xi2 = xin(j+i)
          fi1 = yin(j)
          fi2 = yin(j+1)
          yin(j) = (xo(k)-xi1)/(xi2-xi1)*fi2+(xo(k)-xi2)/(xi1-xi2)*fi1
        enddo
      enddo
      

      y(k) = yin(1) 
      dy(k) = (abs(y(k)-fi1)+abs(y(k)-fi2))/2.d0

      enddo

      return
      end


      subroutine iau_xys00a_fort ( X, Y, S, DATE1, DATE2 )
*f2py intent(in) DATE1, DATE2
*f2py intent(out) X, Y, S
*+
*  - - - - - - - - - - -
*   i a u _ X Y S 0 0 A
*  - - - - - - - - - - -
*
*  For a given TT date, compute the X,Y coordinates of the Celestial
*  Intermediate Pole and the CIO locator s, using the IAU 2000A
*  precession-nutation model.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1,DATE2   d    TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     X,Y           d    Celestial Intermediate Pole (Note 2)
*     S             d    the CIO locator s (Note 2)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The Celestial Intermediate Pole coordinates are the x,y components
*     of the unit vector in the Geocentric Celestial Reference System.
*
*  3) The CIO locator s (in radians) positions the Celestial
*     Intermediate Origin on the equator of the CIP.
*
*  4) A faster, but slightly less accurate result (about 1 mas for X,Y),
*     can be obtained by using instead the iau_XYS00B routine.
*
*  Called:
*     iau_PNM00A   classical NPB matrix, IAU 2000A
*     iau_BPN2XY   extract CIP X,Y coordinates from NPB matrix
*     iau_S00      the CIO locator s, given X,Y, IAU 2000A
*
*  Reference:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*  This revision:  2006 November 13
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, X, Y, S

      DOUBLE PRECISION RBPN(3,3)

      DOUBLE PRECISION iau_S00

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Form the bias-precession-nutation matrix, IAU 2000A.
      CALL iau_PNM00A ( DATE1, DATE2, RBPN )

*  Extract X,Y.
      CALL iau_BPN2XY ( RBPN, X, Y )

*  Obtain s.
      S = iau_S00 ( DATE1, DATE2, X, Y )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      SUBROUTINE iau_PNM00A ( DATE1, DATE2, RBPN )
*+
*  - - - - - - - - - - -
*   i a u _ P N M 0 0 A
*  - - - - - - - - - - -
*
*  Form the matrix of precession-nutation for a given date (including
*  frame bias), equinox-based, IAU 2000A model.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1,DATE2    d       TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     RBPN         d(3,3)    classical NPB matrix (Note 2)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The matrix operates in the sense V(date) = RBPN * V(GCRS), where
*     the p-vector V(date) is with respect to the true equatorial triad
*     of date DATE1+DATE2 and the p-vector V(GCRS) is with respect to
*     the Geocentric Celestial Reference System (IAU, 2000).
*
*  3) A faster, but slightly less accurate result (about 1 mas), can be
*     obtained by using instead the iau_PNM00B routine.
*
*  Called:
*     iau_PN00A    bias/precession/nutation, IAU 2000A
*
*  Reference:
*
*     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
*     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
*     (2000)
*
*  This revision:  2009 December 21
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, RBPN(3,3)

      DOUBLE PRECISION DPSI, DEPS, EPSA, RB(3,3), RP(3,3), RBP(3,3),
     :                 RN(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Obtain the required matrix (discarding other results).
      CALL iau_PN00A ( DATE1, DATE2,
     :                 DPSI, DEPS, EPSA, RB, RP, RBP, RN, RBPN )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      SUBROUTINE iau_PN00A ( DATE1, DATE2,
     :                       DPSI, DEPS, EPSA, RB, RP, RBP, RN, RBPN )
*+
*  - - - - - - - - - -
*   i a u _ P N 0 0 A
*  - - - - - - - - - -
*
*  Precession-nutation, IAU 2000A model:  a multi-purpose routine,
*  supporting classical (equinox-based) use directly and CIO-based
*  use indirectly.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1,DATE2   d       TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     DPSI,DEPS     d       nutation (Note 2)
*     EPSA          d       mean obliquity (Note 3)
*     RB          d(3,3)    frame bias matrix (Note 4)
*     RP          d(3,3)    precession matrix (Note 5)
*     RBP         d(3,3)    bias-precession matrix (Note 6)
*     RN          d(3,3)    nutation matrix (Note 7)
*     RBPN        d(3,3)    GCRS-to-true matrix (Notes 8,9)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The nutation components (luni-solar + planetary, IAU 2000A) in
*     longitude and obliquity are in radians and with respect to the
*     equinox and ecliptic of date.  Free core nutation is omitted;  for
*     the utmost accuracy, use the iau_PN00 routine, where the nutation
*     components are caller-specified.  For faster but slightly less
*     accurate results, use the iau_PN00B routine.
*
*  3) The mean obliquity is consistent with the IAU 2000 precession.
*
*  4) The matrix RB transforms vectors from GCRS to J2000.0 mean equator
*     and equinox by applying frame bias.
*
*  5) The matrix RP transforms vectors from J2000.0 mean equator and
*     equinox to mean equator and equinox of date by applying
*     precession.
*
*  6) The matrix RBP transforms vectors from GCRS to mean equator and
*     equinox of date by applying frame bias then precession.  It is the
*     product RP x RB.
*
*  7) The matrix RN transforms vectors from mean equator and equinox of
*     date to true equator and equinox of date by applying the nutation
*     (luni-solar + planetary).
*
*  8) The matrix RBPN transforms vectors from GCRS to true equator and
*     equinox of date.  It is the product RN x RBP, applying frame bias,
*     precession and nutation in that order.
*
*  9) The X,Y,Z coordinates of the IAU 2000A Celestial Intermediate Pole
*     are elements (3,1-3) of the matrix RBPN.
*
*  Called:
*     iau_NUT00A   nutation, IAU 2000A
*     iau_PN00     bias/precession/nutation results, IAU 2000
*
*  Reference:
*
*     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
*     "Expressions for the Celestial Intermediate Pole and Celestial
*     Ephemeris Origin consistent with the IAU 2000A precession-nutation
*     model", Astron.Astrophys. 400, 1145-1154 (2003).
*
*     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
*          intermediate origin" (CIO) by IAU 2006 Resolution 2.
*
*  This revision:  2010 January 18
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS, EPSA,
     :                 RB(3,3), RP(3,3), RBP(3,3), RN(3,3), RBPN(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Nutation.
      CALL iau_NUT00A ( DATE1, DATE2, DPSI, DEPS )

*  Remaining results.
      CALL iau_PN00 ( DATE1, DATE2, DPSI, DEPS,
     :                EPSA, RB, RP, RBP, RN, RBPN )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      SUBROUTINE iau_NUT00A ( DATE1, DATE2, DPSI, DEPS )
*+
*  - - - - - - - - - - -
*   i a u _ N U T 0 0 A
*  - - - - - - - - - - -
*
*  Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
*  with free core nutation omitted).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2   d    TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     DPSI,DEPS     d    nutation, luni-solar + planetary (Note 2)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The nutation components in longitude and obliquity are in radians
*     and with respect to the equinox and ecliptic of date.  The
*     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
*     value of 84381.448 arcsec.
*
*     Both the luni-solar and planetary nutations are included.  The
*     latter are due to direct planetary nutations and the perturbations
*     of the lunar and terrestrial orbits.
*
*  3) The routine computes the MHB2000 nutation series with the
*     associated corrections for planetary nutations.  It is an
*     implementation of the nutation part of the IAU 2000A precession-
*     nutation model, formally adopted by the IAU General Assembly in
*     2000, namely MHB2000 (Mathews et al. 2002), but with the free core
*     nutation (FCN - see Note 4) omitted.
*
*  4) The full MHB2000 model also contains contributions to the
*     nutations in longitude and obliquity due to the free-excitation of
*     the free-core-nutation during the period 1979-2000.  These FCN
*     terms, which are time-dependent and unpredictable, are NOT
*     included in the present routine and, if required, must be
*     independently computed.  With the FCN corrections included, the
*     present routine delivers a pole which is at current epochs
*     accurate to a few hundred microarcseconds.  The omission of FCN
*     introduces further errors of about that size.
*
*  5) The present routine provides classical nutation.  The MHB2000
*     algorithm, from which it is adapted, deals also with (i) the
*     offsets between the GCRS and mean poles and (ii) the adjustments
*     in longitude and obliquity due to the changed precession rates.
*     These additional functions, namely frame bias and precession
*     adjustments, are supported by the SOFA routines iau_BI00 and
*     iau_PR00.
*
*  6) The MHB2000 algorithm also provides "total" nutations, comprising
*     the arithmetic sum of the frame bias, precession adjustments,
*     luni-solar nutation and planetary nutation.  These total nutations
*     can be used in combination with an existing IAU 1976 precession
*     implementation, such as iau_PMAT76, to deliver GCRS-to-true
*     predictions of sub-mas accuracy at current epochs.  However, there
*     are three shortcomings in the MHB2000 model that must be taken
*     into account if more accurate or definitive results are required
*     (see Wallace 2002):
*
*       (i) The MHB2000 total nutations are simply arithmetic sums,
*           yet in reality the various components are successive Euler
*           rotations.  This slight lack of rigor leads to cross terms
*           that exceed 1 mas after a century.  The rigorous procedure
*           is to form the GCRS-to-true rotation matrix by applying the
*           bias, precession and nutation in that order.
*
*      (ii) Although the precession adjustments are stated to be with
*           respect to Lieske et al. (1977), the MHB2000 model does
*           not specify which set of Euler angles are to be used and
*           how the adjustments are to be applied.  The most literal and
*           straightforward procedure is to adopt the 4-rotation
*           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR to
*           psi_A and DEPSPR to both omega_A and eps_A.
*
*     (iii) The MHB2000 model predates the determination by Chapront
*           et al. (2002) of a 14.6 mas displacement between the J2000.0
*           mean equinox and the origin of the ICRS frame.  It should,
*           however, be noted that neglecting this displacement when
*           calculating star coordinates does not lead to a 14.6 mas
*           change in right ascension, only a small second-order
*           distortion in the pattern of the precession-nutation effect.
*
*     For these reasons, the SOFA routines do not generate the "total
*     nutations" directly, though they can of course easily be generated
*     by calling iau_BI00, iau_PR00 and the present routine and adding
*     the results.
*
*  7) The MHB2000 model contains 41 instances where the same frequency
*     appears multiple times, of which 38 are duplicates and three are
*     triplicates.  To keep the present code close to the original MHB
*     algorithm, this small inefficiency has not been corrected.
*
*  Called:
*     iau_FAL03    mean anomaly of the Moon
*     iau_FAF03    mean argument of the latitude of the Moon
*     iau_FAOM03   mean longitude of the Moon's ascending node
*     iau_FAME03   mean longitude of Mercury
*     iau_FAVE03   mean longitude of Venus
*     iau_FAE03    mean longitude of Earth
*     iau_FAMA03   mean longitude of Mars
*     iau_FAJU03   mean longitude of Jupiter
*     iau_FASA03   mean longitude of Saturn
*     iau_FAUR03   mean longitude of Uranus
*     iau_FAPA03   general accumulated precession in longitude
*
*  References:
*
*     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
*     Astron.Astrophys. 387, 700
*
*     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
*     Astron.Astrophys. 58, 1-16
*
*     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
*     107, B4.  The MHB_2000 code itself was obtained on 9th September
*     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*     Wallace, P.T., "Software for Implementing the IAU 2000
*     Resolutions", in IERS Workshop 5.1 (2002)
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Units of 0.1 microarcsecond to radians
      DOUBLE PRECISION U2R
      PARAMETER ( U2R = DAS2R/1D7 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG,
     :                 AL, ALSU, AF, AD, AOM, ALME, ALVE, ALEA, ALMA,
     :                 ALJU, ALSA, ALUR, ALNE, APA, DPSILS, DEPSLS,
     :                 DPSIPL, DEPSPL
      DOUBLE PRECISION iau_FAL03, iau_FAF03, iau_FAOM03, iau_FAME03,
     :                 iau_FAVE03, iau_FAE03, iau_FAMA03, iau_FAJU03,
     :                 iau_FASA03, iau_FAUR03, iau_FAPA03

*  -------------------------
*  Luni-Solar nutation model
*  -------------------------

*  Number of terms in the luni-solar nutation model
      INTEGER NLS
      PARAMETER ( NLS = 678 )

*  Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)

*  Longitude and obliquity coefficients
      DOUBLE PRECISION CLS(6,NLS)

*  ---------------
*  Planetary terms
*  ---------------

*  Number of terms in the planetary nutation model
      INTEGER NPL
      PARAMETER ( NPL = 687 )

*  Coefficients for fundamental arguments
      INTEGER NAPL(14,NPL)

*  Longitude and obliquity coefficients
      INTEGER ICPL(4,NPL)

*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------

*
*  Luni-Solar argument multipliers
*               L     L'    F     D     Om

      DATA ( ( NALS(I,J), I=1,5 ), J=  1, 10 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 11, 20 ) /
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 21, 30 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 31, 40 ) /
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 41, 50 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 51, 60 ) /
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 61, 70 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 71, 80 ) /
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2,
     :         -2,    0,    2,    4,    2,
     :         -1,    0,    4,    0,    2,
     :          2,    0,    2,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 81, 90 ) /
     :          2,    0,    2,    2,    2,
     :          1,    0,    0,    2,    1,
     :          3,    0,    0,    0,    0,
     :          3,    0,    2,   -2,    2,
     :          0,    0,    4,   -2,    2,
     :          0,    1,    2,    0,    1,
     :          0,    0,   -2,    2,    1,
     :          0,    0,    2,   -2,    3,
     :         -1,    0,    0,    4,    0,
     :          2,    0,   -2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 91,100 ) /
     :         -2,    0,    0,    4,    0,
     :         -1,   -1,    0,    2,    1,
     :         -1,    0,    0,    1,    1,
     :          0,    1,    0,    0,    2,
     :          0,    0,   -2,    0,    1,
     :          0,   -1,    2,    0,    1,
     :          0,    0,    2,   -1,    2,
     :          0,    0,    2,    4,    2,
     :         -2,   -1,    0,    2,    0,
     :          1,    1,    0,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=101,110 ) /
     :         -1,    1,    0,    2,    0,
     :         -1,    1,    0,    1,    2,
     :          1,   -1,    0,    0,    1,
     :          1,   -1,    2,    2,    2,
     :         -1,    1,    2,    2,    2,
     :          3,    0,    2,    0,    1,
     :          0,    1,   -2,    2,    0,
     :         -1,    0,    0,   -2,    1,
     :          0,    1,    2,    2,    2,
     :         -1,   -1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=111,120 ) /
     :          0,   -1,    0,    0,    2,
     :          1,    0,    2,   -4,    1,
     :         -1,    0,   -2,    2,    0,
     :          0,   -1,    2,    2,    1,
     :          2,   -1,    2,    0,    2,
     :          0,    0,    0,    2,    2,
     :          1,   -1,    2,    0,    1,
     :         -1,    1,    2,    0,    2,
     :          0,    1,    0,    2,    0,
     :          0,   -1,   -2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=121,130 ) /
     :          0,    3,    2,   -2,    2,
     :          0,    0,    0,    1,    1,
     :         -1,    0,    2,    2,    0,
     :          2,    1,    2,    0,    2,
     :          1,    1,    0,    0,    1,
     :          1,    1,    2,    0,    1,
     :          2,    0,    0,    2,    0,
     :          1,    0,   -2,    2,    0,
     :         -1,    0,    0,    2,    2,
     :          0,    1,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=131,140 ) /
     :          0,    1,    0,   -2,    1,
     :         -1,    0,    2,   -2,    2,
     :          0,    0,    0,   -1,    1,
     :         -1,    1,    0,    0,    1,
     :          1,    0,    2,   -1,    2,
     :          1,   -1,    0,    2,    0,
     :          0,    0,    0,    4,    0,
     :          1,    0,    2,    1,    2,
     :          0,    0,    2,    1,    1,
     :          1,    0,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=141,150 ) /
     :         -1,    0,    2,    4,    1,
     :          1,    0,   -2,    0,    1,
     :          1,    1,    2,   -2,    1,
     :          0,    0,    2,    2,    0,
     :         -1,    0,    2,   -1,    1,
     :         -2,    0,    2,    2,    1,
     :          4,    0,    2,    0,    2,
     :          2,   -1,    0,    0,    0,
     :          2,    1,    2,   -2,    2,
     :          0,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=151,160 ) /
     :          1,    0,    4,   -2,    2,
     :         -1,   -1,    0,    0,    1,
     :          0,    1,    0,    2,    1,
     :         -2,    0,    2,    4,    1,
     :          2,    0,    2,    0,    0,
     :          1,    0,    0,    1,    0,
     :         -1,    0,    0,    4,    1,
     :         -1,    0,    4,    0,    1,
     :          2,    0,    2,    2,    1,
     :          0,    0,    2,   -3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=161,170 ) /
     :         -1,   -2,    0,    2,    0,
     :          2,    1,    0,    0,    0,
     :          0,    0,    4,    0,    2,
     :          0,    0,    0,    0,    3,
     :          0,    3,    0,    0,    0,
     :          0,    0,    2,   -4,    1,
     :          0,   -1,    0,    2,    1,
     :          0,    0,    0,    4,    1,
     :         -1,   -1,    2,    4,    2,
     :          1,    0,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=171,180 ) /
     :         -2,    2,    0,    2,    0,
     :         -2,   -1,    2,    0,    1,
     :         -2,    0,    0,    2,    2,
     :         -1,   -1,    2,    0,    2,
     :          0,    0,    4,   -2,    1,
     :          3,    0,    2,   -2,    1,
     :         -2,   -1,    0,    2,    1,
     :          1,    0,    0,   -1,    1,
     :          0,   -2,    0,    2,    0,
     :         -2,    0,    0,    4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=181,190 ) /
     :         -3,    0,    0,    0,    1,
     :          1,    1,    2,    2,    2,
     :          0,    0,    2,    4,    1,
     :          3,    0,    2,    2,    2,
     :         -1,    1,    2,   -2,    1,
     :          2,    0,    0,   -4,    1,
     :          0,    0,    0,   -2,    2,
     :          2,    0,    2,   -4,    1,
     :         -1,    1,    0,    2,    1,
     :          0,    0,    2,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=191,200 ) /
     :          0,   -2,    2,    2,    2,
     :          2,    0,    0,    2,    1,
     :          4,    0,    2,   -2,    2,
     :          2,    0,    0,   -2,    2,
     :          0,    2,    0,    0,    1,
     :          1,    0,    0,   -4,    1,
     :          0,    2,    2,   -2,    1,
     :         -3,    0,    0,    4,    0,
     :         -1,    1,    2,    0,    1,
     :         -1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=201,210 ) /
     :         -1,   -2,    2,    2,    2,
     :         -2,   -1,    2,    4,    2,
     :          1,   -1,    2,    2,    1,
     :         -2,    1,    0,    2,    0,
     :         -2,    1,    2,    0,    1,
     :          2,    1,    0,   -2,    1,
     :         -3,    0,    2,    0,    1,
     :         -2,    0,    2,   -2,    1,
     :         -1,    1,    0,    2,    2,
     :          0,   -1,    2,   -1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=211,220 ) /
     :         -1,    0,    4,   -2,    2,
     :          0,   -2,    2,    0,    2,
     :         -1,    0,    2,    1,    2,
     :          2,    0,    0,    0,    2,
     :          0,    0,    2,    0,    3,
     :         -2,    0,    4,    0,    2,
     :         -1,    0,   -2,    0,    1,
     :         -1,    1,    2,    2,    1,
     :          3,    0,    0,    0,    1,
     :         -1,    0,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=221,230 ) /
     :          2,   -1,    2,    0,    1,
     :          0,    1,    2,    2,    1,
     :          0,   -1,    2,    4,    2,
     :          2,   -1,    2,    2,    2,
     :          0,    2,   -2,    2,    0,
     :         -1,   -1,    2,   -1,    1,
     :          0,   -2,    0,    0,    1,
     :          1,    0,    2,   -4,    2,
     :          1,   -1,    0,   -2,    1,
     :         -1,   -1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=231,240 ) /
     :          1,   -1,    2,   -2,    2,
     :         -2,   -1,    0,    4,    0,
     :         -1,    0,    0,    3,    0,
     :         -2,   -1,    2,    2,    2,
     :          0,    2,    2,    0,    2,
     :          1,    1,    0,    2,    0,
     :          2,    0,    2,   -1,    2,
     :          1,    0,    2,    1,    1,
     :          4,    0,    0,    0,    0,
     :          2,    1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=241,250 ) /
     :          3,   -1,    2,    0,    2,
     :         -2,    2,    0,    2,    1,
     :          1,    0,    2,   -3,    1,
     :          1,    1,    2,   -4,    1,
     :         -1,   -1,    2,   -2,    1,
     :          0,   -1,    0,   -1,    1,
     :          0,   -1,    0,   -2,    1,
     :         -2,    0,    0,    0,    2,
     :         -2,    0,   -2,    2,    0,
     :         -1,    0,   -2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=251,260 ) /
     :          1,   -2,    0,    0,    0,
     :          0,    1,    0,    1,    1,
     :         -1,    2,    0,    2,    0,
     :          1,   -1,    2,   -2,    1,
     :          1,    2,    2,   -2,    2,
     :          2,   -1,    2,   -2,    2,
     :          1,    0,    2,   -1,    1,
     :          2,    1,    2,   -2,    1,
     :         -2,    0,    0,   -2,    1,
     :          1,   -2,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=261,270 ) /
     :          0,    1,    2,    1,    1,
     :          1,    0,    4,   -2,    1,
     :         -2,    0,    4,    2,    2,
     :          1,    1,    2,    1,    2,
     :          1,    0,    0,    4,    0,
     :          1,    0,    2,    2,    0,
     :          2,    0,    2,    1,    2,
     :          3,    1,    2,    0,    2,
     :          4,    0,    2,    0,    1,
     :         -2,   -1,    2,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=271,280 ) /
     :          0,    1,   -2,    2,    1,
     :          1,    0,   -2,    1,    0,
     :          0,   -1,   -2,    2,    1,
     :          2,   -1,    0,   -2,    1,
     :         -1,    0,    2,   -1,    2,
     :          1,    0,    2,   -3,    2,
     :          0,    1,    2,   -2,    3,
     :          0,    0,    2,   -3,    1,
     :         -1,    0,   -2,    2,    1,
     :          0,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=281,290 ) /
     :         -2,    1,    0,    0,    1,
     :         -1,    0,    0,   -1,    1,
     :          2,    0,    2,   -4,    2,
     :          0,    0,    4,   -4,    4,
     :          0,    0,    4,   -4,    2,
     :         -1,   -2,    0,    2,    1,
     :         -2,    0,    0,    3,    0,
     :          1,    0,   -2,    2,    1,
     :         -3,    0,    2,    2,    2,
     :         -3,    0,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=291,300 ) /
     :         -2,    0,    2,    2,    0,
     :          2,   -1,    0,    0,    1,
     :         -2,    1,    2,    2,    2,
     :          1,    1,    0,    1,    0,
     :          0,    1,    4,   -2,    2,
     :         -1,    1,    0,   -2,    1,
     :          0,    0,    0,   -4,    1,
     :          1,   -1,    0,    2,    1,
     :          1,    1,    0,    2,    1,
     :         -1,    2,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=301,310 ) /
     :          3,    1,    2,   -2,    2,
     :          0,   -1,    0,    4,    0,
     :          2,   -1,    0,    2,    0,
     :          0,    0,    4,    0,    1,
     :          2,    0,    4,   -2,    2,
     :         -1,   -1,    2,    4,    1,
     :          1,    0,    0,    4,    1,
     :          1,   -2,    2,    2,    2,
     :          0,    0,    2,    3,    2,
     :         -1,    1,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=311,320 ) /
     :          3,    0,    0,    2,    0,
     :         -1,    0,    4,    2,    2,
     :          1,    1,    2,    2,    1,
     :         -2,    0,    2,    6,    2,
     :          2,    1,    2,    2,    2,
     :         -1,    0,    2,    6,    2,
     :          1,    0,    2,    4,    1,
     :          2,    0,    2,    4,    2,
     :          1,    1,   -2,    1,    0,
     :         -3,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=321,330 ) /
     :          2,    0,   -2,    0,    2,
     :         -1,    0,    0,    1,    2,
     :         -4,    0,    2,    2,    1,
     :         -1,   -1,    0,    1,    0,
     :          0,    0,   -2,    2,    2,
     :          1,    0,    0,   -1,    2,
     :          0,   -1,    2,   -2,    3,
     :         -2,    1,    2,    0,    0,
     :          0,    0,    2,   -2,    4,
     :         -2,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=331,340 ) /
     :         -2,    0,   -2,    4,    0,
     :          0,   -2,   -2,    2,    0,
     :          1,    2,    0,   -2,    1,
     :          3,    0,    0,   -4,    1,
     :         -1,    1,    2,   -2,    2,
     :          1,   -1,    2,   -4,    1,
     :          1,    1,    0,   -2,    2,
     :         -3,    0,    2,    0,    0,
     :         -3,    0,    2,    0,    2,
     :         -2,    0,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=341,350 ) /
     :          0,    0,   -2,    1,    0,
     :         -3,    0,    0,    2,    1,
     :         -1,   -1,   -2,    2,    0,
     :          0,    1,    2,   -4,    1,
     :          2,    1,    0,   -4,    1,
     :          0,    2,    0,   -2,    1,
     :          1,    0,    0,   -3,    1,
     :         -2,    0,    2,   -2,    2,
     :         -2,   -1,    0,    0,    1,
     :         -4,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=351,360 ) /
     :          1,    1,    0,   -4,    1,
     :         -1,    0,    2,   -4,    1,
     :          0,    0,    4,   -4,    1,
     :          0,    3,    2,   -2,    2,
     :         -3,   -1,    0,    4,    0,
     :         -3,    0,    0,    4,    1,
     :          1,   -1,   -2,    2,    0,
     :         -1,   -1,    0,    2,    2,
     :          1,   -2,    0,    0,    1,
     :          1,   -1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=361,370 ) /
     :          0,    0,    0,    1,    2,
     :         -1,   -1,    2,    0,    0,
     :          1,   -2,    2,   -2,    2,
     :          0,   -1,    2,   -1,    1,
     :         -1,    0,    2,    0,    3,
     :          1,    1,    0,    0,    2,
     :         -1,    1,    2,    0,    0,
     :          1,    2,    0,    0,    0,
     :         -1,    2,    2,    0,    2,
     :         -1,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=371,380 ) /
     :          3,    0,    2,   -4,    2,
     :          1,    2,    2,   -2,    1,
     :          1,    0,    4,   -4,    2,
     :         -2,   -1,    0,    4,    1,
     :          0,   -1,    0,    2,    2,
     :         -2,    1,    0,    4,    0,
     :         -2,   -1,    2,    2,    1,
     :          2,    0,   -2,    2,    0,
     :          1,    0,    0,    1,    1,
     :          0,    1,    0,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=381,390 ) /
     :          1,   -1,    2,   -1,    2,
     :         -2,    0,    4,    0,    1,
     :          2,    1,    0,    0,    1,
     :          0,    1,    2,    0,    0,
     :          0,   -1,    4,   -2,    2,
     :          0,    0,    4,   -2,    4,
     :          0,    2,    2,    0,    1,
     :         -3,    0,    0,    6,    0,
     :         -1,   -1,    0,    4,    1,
     :          1,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=391,400 ) /
     :         -1,    0,    0,    4,    2,
     :         -1,   -2,    2,    2,    1,
     :         -1,    0,    0,   -2,    2,
     :          1,    0,   -2,   -2,    1,
     :          0,    0,   -2,   -2,    1,
     :         -2,    0,   -2,    0,    1,
     :          0,    0,    0,    3,    1,
     :          0,    0,    0,    3,    0,
     :         -1,    1,    0,    4,    0,
     :         -1,   -1,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=401,410 ) /
     :         -2,    0,    2,    3,    2,
     :          1,    0,    0,    2,    2,
     :          0,   -1,    2,    1,    2,
     :          3,   -1,    0,    0,    0,
     :          2,    0,    0,    1,    0,
     :          1,   -1,    2,    0,    0,
     :          0,    0,    2,    1,    0,
     :          1,    0,    2,    0,    3,
     :          3,    1,    0,    0,    0,
     :          3,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=411,420 ) /
     :          2,    0,    2,   -1,    1,
     :          1,    1,    2,    0,    0,
     :          0,    0,    4,   -1,    2,
     :          1,    2,    2,    0,    2,
     :         -2,    0,    0,    6,    0,
     :          0,   -1,    0,    4,    1,
     :         -2,   -1,    2,    4,    1,
     :          0,   -2,    2,    2,    1,
     :          0,   -1,    2,    2,    0,
     :         -1,    0,    2,    3,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=421,430 ) /
     :         -2,    1,    2,    4,    2,
     :          2,    0,    0,    2,    2,
     :          2,   -2,    2,    0,    2,
     :         -1,    1,    2,    3,    2,
     :          3,    0,    2,   -1,    2,
     :          4,    0,    2,   -2,    1,
     :         -1,    0,    0,    6,    0,
     :         -1,   -2,    2,    4,    2,
     :         -3,    0,    2,    6,    2,
     :         -1,    0,    2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=431,440 ) /
     :          3,    0,    0,    2,    1,
     :          3,   -1,    2,    0,    1,
     :          3,    0,    2,    0,    0,
     :          1,    0,    4,    0,    2,
     :          5,    0,    2,   -2,    2,
     :          0,   -1,    2,    4,    1,
     :          2,   -1,    2,    2,    1,
     :          0,    1,    2,    4,    2,
     :          1,   -1,    2,    4,    2,
     :          3,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=441,450 ) /
     :          3,    0,    2,    2,    1,
     :          5,    0,    2,    0,    2,
     :          0,    0,    2,    6,    2,
     :          4,    0,    2,    2,    2,
     :          0,   -1,    1,   -1,    1,
     :         -1,    0,    1,    0,    3,
     :          0,   -2,    2,   -2,    3,
     :          1,    0,   -1,    0,    1,
     :          2,   -2,    0,   -2,    1,
     :         -1,    0,    1,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=451,460 ) /
     :         -1,    0,    1,    0,    1,
     :         -1,   -1,    2,   -1,    2,
     :         -2,    2,    0,    2,    2,
     :         -1,    0,    1,    0,    0,
     :         -4,    1,    2,    2,    2,
     :         -3,    0,    2,    1,    1,
     :         -2,   -1,    2,    0,    2,
     :          1,    0,   -2,    1,    1,
     :          2,   -1,   -2,    0,    1,
     :         -4,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=461,470 ) /
     :         -3,    1,    0,    3,    0,
     :         -1,    0,   -1,    2,    0,
     :          0,   -2,    0,    0,    2,
     :          0,   -2,    0,    0,    2,
     :         -3,    0,    0,    3,    0,
     :         -2,   -1,    0,    2,    2,
     :         -1,    0,   -2,    3,    0,
     :         -4,    0,    0,    4,    0,
     :          2,    1,   -2,    0,    1,
     :          2,   -1,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=471,480 ) /
     :          0,    0,    1,   -1,    0,
     :         -1,    2,    0,    1,    0,
     :         -2,    1,    2,    0,    2,
     :          1,    1,    0,   -1,    1,
     :          1,    0,    1,   -2,    1,
     :          0,    2,    0,    0,    2,
     :          1,   -1,    2,   -3,    1,
     :         -1,    1,    2,   -1,    1,
     :         -2,    0,    4,   -2,    2,
     :         -2,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=481,490 ) /
     :         -2,   -2,    0,    2,    1,
     :         -2,    0,   -2,    4,    0,
     :          1,    2,    2,   -4,    1,
     :          1,    1,    2,   -4,    2,
     :         -1,    2,    2,   -2,    1,
     :          2,    0,    0,   -3,    1,
     :         -1,    2,    0,    0,    1,
     :          0,    0,    0,   -2,    0,
     :         -1,   -1,    2,   -2,    2,
     :         -1,    1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=491,500 ) /
     :          0,    0,    0,   -1,    2,
     :         -2,    1,    0,    1,    0,
     :          1,   -2,    0,   -2,    1,
     :          1,    0,   -2,    0,    2,
     :         -3,    1,    0,    2,    0,
     :         -1,    1,   -2,    2,    0,
     :         -1,   -1,    0,    0,    2,
     :         -3,    0,    0,    2,    0,
     :         -3,   -1,    0,    2,    0,
     :          2,    0,    2,   -6,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=501,510 ) /
     :          0,    1,    2,   -4,    2,
     :          2,    0,    0,   -4,    2,
     :         -2,    1,    2,   -2,    1,
     :          0,   -1,    2,   -4,    1,
     :          0,    1,    0,   -2,    2,
     :         -1,    0,    0,   -2,    0,
     :          2,    0,   -2,   -2,    1,
     :         -4,    0,    2,    0,    1,
     :         -1,   -1,    0,   -1,    1,
     :          0,    0,   -2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=511,520 ) /
     :         -3,    0,    0,    1,    0,
     :         -1,    0,   -2,    1,    0,
     :         -2,    0,   -2,    2,    1,
     :          0,    0,   -4,    2,    0,
     :         -2,   -1,   -2,    2,    0,
     :          1,    0,    2,   -6,    1,
     :         -1,    0,    2,   -4,    2,
     :          1,    0,    0,   -4,    2,
     :          2,    1,    2,   -4,    2,
     :          2,    1,    2,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=521,530 ) /
     :          0,    1,    4,   -4,    4,
     :          0,    1,    4,   -4,    2,
     :         -1,   -1,   -2,    4,    0,
     :         -1,   -3,    0,    2,    0,
     :         -1,    0,   -2,    4,    1,
     :         -2,   -1,    0,    3,    0,
     :          0,    0,   -2,    3,    0,
     :         -2,    0,    0,    3,    1,
     :          0,   -1,    0,    1,    0,
     :         -3,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=531,540 ) /
     :          1,    1,   -2,    2,    0,
     :         -1,    1,    0,    2,    2,
     :          1,   -2,    2,   -2,    1,
     :          0,    0,    1,    0,    2,
     :          0,    0,    1,    0,    1,
     :          0,    0,    1,    0,    0,
     :         -1,    2,    0,    2,    1,
     :          0,    0,    2,    0,    2,
     :         -2,    0,    2,    0,    2,
     :          2,    0,    0,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=541,550 ) /
     :          3,    0,    0,   -2,    1,
     :          1,    0,    2,   -2,    3,
     :          1,    2,    0,    0,    1,
     :          2,    0,    2,   -3,    2,
     :         -1,    1,    4,   -2,    2,
     :         -2,   -2,    0,    4,    0,
     :          0,   -3,    0,    2,    0,
     :          0,    0,   -2,    4,    0,
     :         -1,   -1,    0,    3,    0,
     :         -2,    0,    0,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=551,560 ) /
     :         -1,    0,    0,    3,    1,
     :          2,   -2,    0,    0,    0,
     :          1,   -1,    0,    1,    0,
     :         -1,    0,    0,    2,    0,
     :          0,   -2,    2,    0,    1,
     :         -1,    0,    1,    2,    1,
     :         -1,    1,    0,    3,    0,
     :         -1,   -1,    2,    1,    2,
     :          0,   -1,    2,    0,    0,
     :         -2,    1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=561,570 ) /
     :          2,   -2,    2,   -2,    2,
     :          1,    1,    0,    1,    1,
     :          1,    0,    1,    0,    1,
     :          1,    0,    1,    0,    0,
     :          0,    2,    0,    2,    0,
     :          2,   -1,    2,   -2,    1,
     :          0,   -1,    4,   -2,    1,
     :          0,    0,    4,   -2,    3,
     :          0,    1,    4,   -2,    1,
     :          4,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=571,580 ) /
     :          2,    2,    2,   -2,    2,
     :          2,    0,    4,   -4,    2,
     :         -1,   -2,    0,    4,    0,
     :         -1,   -3,    2,    2,    2,
     :         -3,    0,    2,    4,    2,
     :         -3,    0,    2,   -2,    1,
     :         -1,   -1,    0,   -2,    1,
     :         -3,    0,    0,    0,    2,
     :         -3,    0,   -2,    2,    0,
     :          0,    1,    0,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=581,590 ) /
     :         -2,    1,    0,   -2,    1,
     :         -4,    0,    0,    0,    1,
     :         -1,    0,    0,   -4,    1,
     :         -3,    0,    0,   -2,    1,
     :          0,    0,    0,    3,    2,
     :         -1,    1,    0,    4,    1,
     :          1,   -2,    2,    0,    1,
     :          0,    1,    0,    3,    0,
     :         -1,    0,    2,    2,    3,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=591,600 ) /
     :         -2,    0,    2,    2,    2,
     :         -1,    1,    2,    2,    0,
     :          3,    0,    0,    0,    2,
     :          2,    1,    0,    1,    0,
     :          2,   -1,    2,   -1,    2,
     :          0,    0,    2,    0,    1,
     :          0,    0,    3,    0,    3,
     :          0,    0,    3,    0,    2,
     :         -1,    2,    2,    2,    1,
     :         -1,    0,    4,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=601,610 ) /
     :          1,    2,    2,    0,    1,
     :          3,    1,    2,   -2,    1,
     :          1,    1,    4,   -2,    2,
     :         -2,   -1,    0,    6,    0,
     :          0,   -2,    0,    4,    0,
     :         -2,    0,    0,    6,    1,
     :         -2,   -2,    2,    4,    2,
     :          0,   -3,    2,    2,    2,
     :          0,    0,    0,    4,    2,
     :         -1,   -1,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=611,620 ) /
     :         -2,    0,    2,    4,    0,
     :          2,   -1,    0,    2,    1,
     :          1,    0,    0,    3,    0,
     :          0,    1,    0,    4,    1,
     :          0,    1,    0,    4,    0,
     :          1,   -1,    2,    1,    2,
     :          0,    0,    2,    2,    3,
     :          1,    0,    2,    2,    2,
     :         -1,    0,    2,    2,    2,
     :         -2,    0,    4,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=621,630 ) /
     :          2,    1,    0,    2,    1,
     :          2,    1,    0,    2,    0,
     :          2,   -1,    2,    0,    0,
     :          1,    0,    2,    1,    0,
     :          0,    1,    2,    2,    0,
     :          2,    0,    2,    0,    3,
     :          3,    0,    2,    0,    2,
     :          1,    0,    2,    0,    2,
     :          1,    0,    3,    0,    3,
     :          1,    1,    2,    1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=631,640 ) /
     :          0,    2,    2,    2,    2,
     :          2,    1,    2,    0,    0,
     :          2,    0,    4,   -2,    1,
     :          4,    1,    2,   -2,    2,
     :         -1,   -1,    0,    6,    0,
     :         -3,   -1,    2,    6,    2,
     :         -1,    0,    0,    6,    1,
     :         -3,    0,    2,    6,    1,
     :          1,   -1,    0,    4,    1,
     :          1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=641,650 ) /
     :         -2,    0,    2,    5,    2,
     :          1,   -2,    2,    2,    1,
     :          3,   -1,    0,    2,    0,
     :          1,   -1,    2,    2,    0,
     :          0,    0,    2,    3,    1,
     :         -1,    1,    2,    4,    1,
     :          0,    1,    2,    3,    2,
     :         -1,    0,    4,    2,    1,
     :          2,    0,    2,    1,    1,
     :          5,    0,    0,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=651,660 ) /
     :          2,    1,    2,    1,    2,
     :          1,    0,    4,    0,    1,
     :          3,    1,    2,    0,    1,
     :          3,    0,    4,   -2,    2,
     :         -2,   -1,    2,    6,    2,
     :          0,    0,    0,    6,    0,
     :          0,   -2,    2,    4,    2,
     :         -2,    0,    2,    6,    1,
     :          2,    0,    0,    4,    1,
     :          2,    0,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=661,670 ) /
     :          2,   -2,    2,    2,    2,
     :          0,    0,    2,    4,    0,
     :          1,    0,    2,    3,    2,
     :          4,    0,    0,    2,    0,
     :          2,    0,    2,    2,    0,
     :          0,    0,    4,    2,    2,
     :          4,   -1,    2,    0,    2,
     :          3,    0,    2,    1,    2,
     :          2,    1,    2,    2,    1,
     :          4,    1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=671,678 ) /
     :         -1,   -1,    2,    6,    2,
     :         -1,    0,    2,    6,    1,
     :          1,   -1,    2,    4,    1,
     :          1,    1,    2,    4,    2,
     :          3,    1,    2,    2,    2,
     :          5,    0,    2,    0,    1,
     :          2,   -1,    2,    4,    2,
     :          2,    0,    2,    4,    1 /

*
*  Luni-Solar nutation coefficients, unit 1e-7 arcsec
*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
*

      DATA ( ( CLS(I,J), I=1,6 ), J=  1, 10 ) /
     : -172064161D0, -174666D0,  33386D0, 92052331D0,  9086D0, 15377D0,
     :  -13170906D0,   -1675D0, -13696D0,  5730336D0, -3015D0, -4587D0,
     :   -2276413D0,    -234D0,   2796D0,   978459D0,  -485D0,  1374D0,
     :    2074554D0,     207D0,   -698D0,  -897492D0,   470D0,  -291D0,
     :    1475877D0,   -3633D0,  11817D0,    73871D0,  -184D0, -1924D0,
     :    -516821D0,    1226D0,   -524D0,   224386D0,  -677D0,  -174D0,
     :     711159D0,      73D0,   -872D0,    -6750D0,     0D0,   358D0,
     :    -387298D0,    -367D0,    380D0,   200728D0,    18D0,   318D0,
     :    -301461D0,     -36D0,    816D0,   129025D0,   -63D0,   367D0,
     :     215829D0,    -494D0,    111D0,   -95929D0,   299D0,   132D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 11, 20 ) /
     :     128227D0,     137D0,    181D0,   -68982D0,    -9D0,    39D0,
     :     123457D0,      11D0,     19D0,   -53311D0,    32D0,    -4D0,
     :     156994D0,      10D0,   -168D0,    -1235D0,     0D0,    82D0,
     :      63110D0,      63D0,     27D0,   -33228D0,     0D0,    -9D0,
     :     -57976D0,     -63D0,   -189D0,    31429D0,     0D0,   -75D0,
     :     -59641D0,     -11D0,    149D0,    25543D0,   -11D0,    66D0,
     :     -51613D0,     -42D0,    129D0,    26366D0,     0D0,    78D0,
     :      45893D0,      50D0,     31D0,   -24236D0,   -10D0,    20D0,
     :      63384D0,      11D0,   -150D0,    -1220D0,     0D0,    29D0,
     :     -38571D0,      -1D0,    158D0,    16452D0,   -11D0,    68D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 21, 30 ) /
     :      32481D0,       0D0,      0D0,   -13870D0,     0D0,     0D0,
     :     -47722D0,       0D0,    -18D0,      477D0,     0D0,   -25D0,
     :     -31046D0,      -1D0,    131D0,    13238D0,   -11D0,    59D0,
     :      28593D0,       0D0,     -1D0,   -12338D0,    10D0,    -3D0,
     :      20441D0,      21D0,     10D0,   -10758D0,     0D0,    -3D0,
     :      29243D0,       0D0,    -74D0,     -609D0,     0D0,    13D0,
     :      25887D0,       0D0,    -66D0,     -550D0,     0D0,    11D0,
     :     -14053D0,     -25D0,     79D0,     8551D0,    -2D0,   -45D0,
     :      15164D0,      10D0,     11D0,    -8001D0,     0D0,    -1D0,
     :     -15794D0,      72D0,    -16D0,     6850D0,   -42D0,    -5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 31, 40 ) /
     :      21783D0,       0D0,     13D0,     -167D0,     0D0,    13D0,
     :     -12873D0,     -10D0,    -37D0,     6953D0,     0D0,   -14D0,
     :     -12654D0,      11D0,     63D0,     6415D0,     0D0,    26D0,
     :     -10204D0,       0D0,     25D0,     5222D0,     0D0,    15D0,
     :      16707D0,     -85D0,    -10D0,      168D0,    -1D0,    10D0,
     :      -7691D0,       0D0,     44D0,     3268D0,     0D0,    19D0,
     :     -11024D0,       0D0,    -14D0,      104D0,     0D0,     2D0,
     :       7566D0,     -21D0,    -11D0,    -3250D0,     0D0,    -5D0,
     :      -6637D0,     -11D0,     25D0,     3353D0,     0D0,    14D0,
     :      -7141D0,      21D0,      8D0,     3070D0,     0D0,     4D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 41, 50 ) /
     :      -6302D0,     -11D0,      2D0,     3272D0,     0D0,     4D0,
     :       5800D0,      10D0,      2D0,    -3045D0,     0D0,    -1D0,
     :       6443D0,       0D0,     -7D0,    -2768D0,     0D0,    -4D0,
     :      -5774D0,     -11D0,    -15D0,     3041D0,     0D0,    -5D0,
     :      -5350D0,       0D0,     21D0,     2695D0,     0D0,    12D0,
     :      -4752D0,     -11D0,     -3D0,     2719D0,     0D0,    -3D0,
     :      -4940D0,     -11D0,    -21D0,     2720D0,     0D0,    -9D0,
     :       7350D0,       0D0,     -8D0,      -51D0,     0D0,     4D0,
     :       4065D0,       0D0,      6D0,    -2206D0,     0D0,     1D0,
     :       6579D0,       0D0,    -24D0,     -199D0,     0D0,     2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 51, 60 ) /
     :       3579D0,       0D0,      5D0,    -1900D0,     0D0,     1D0,
     :       4725D0,       0D0,     -6D0,      -41D0,     0D0,     3D0,
     :      -3075D0,       0D0,     -2D0,     1313D0,     0D0,    -1D0,
     :      -2904D0,       0D0,     15D0,     1233D0,     0D0,     7D0,
     :       4348D0,       0D0,    -10D0,      -81D0,     0D0,     2D0,
     :      -2878D0,       0D0,      8D0,     1232D0,     0D0,     4D0,
     :      -4230D0,       0D0,      5D0,      -20D0,     0D0,    -2D0,
     :      -2819D0,       0D0,      7D0,     1207D0,     0D0,     3D0,
     :      -4056D0,       0D0,      5D0,       40D0,     0D0,    -2D0,
     :      -2647D0,       0D0,     11D0,     1129D0,     0D0,     5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 61, 70 ) /
     :      -2294D0,       0D0,    -10D0,     1266D0,     0D0,    -4D0,
     :       2481D0,       0D0,     -7D0,    -1062D0,     0D0,    -3D0,
     :       2179D0,       0D0,     -2D0,    -1129D0,     0D0,    -2D0,
     :       3276D0,       0D0,      1D0,       -9D0,     0D0,     0D0,
     :      -3389D0,       0D0,      5D0,       35D0,     0D0,    -2D0,
     :       3339D0,       0D0,    -13D0,     -107D0,     0D0,     1D0,
     :      -1987D0,       0D0,     -6D0,     1073D0,     0D0,    -2D0,
     :      -1981D0,       0D0,      0D0,      854D0,     0D0,     0D0,
     :       4026D0,       0D0,   -353D0,     -553D0,     0D0,  -139D0,
     :       1660D0,       0D0,     -5D0,     -710D0,     0D0,    -2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 71, 80 ) /
     :      -1521D0,       0D0,      9D0,      647D0,     0D0,     4D0,
     :       1314D0,       0D0,      0D0,     -700D0,     0D0,     0D0,
     :      -1283D0,       0D0,      0D0,      672D0,     0D0,     0D0,
     :      -1331D0,       0D0,      8D0,      663D0,     0D0,     4D0,
     :       1383D0,       0D0,     -2D0,     -594D0,     0D0,    -2D0,
     :       1405D0,       0D0,      4D0,     -610D0,     0D0,     2D0,
     :       1290D0,       0D0,      0D0,     -556D0,     0D0,     0D0,
     :      -1214D0,       0D0,      5D0,      518D0,     0D0,     2D0,
     :       1146D0,       0D0,     -3D0,     -490D0,     0D0,    -1D0,
     :       1019D0,       0D0,     -1D0,     -527D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 81, 90 ) /
     :      -1100D0,       0D0,      9D0,      465D0,     0D0,     4D0,
     :       -970D0,       0D0,      2D0,      496D0,     0D0,     1D0,
     :       1575D0,       0D0,     -6D0,      -50D0,     0D0,     0D0,
     :        934D0,       0D0,     -3D0,     -399D0,     0D0,    -1D0,
     :        922D0,       0D0,     -1D0,     -395D0,     0D0,    -1D0,
     :        815D0,       0D0,     -1D0,     -422D0,     0D0,    -1D0,
     :        834D0,       0D0,      2D0,     -440D0,     0D0,     1D0,
     :       1248D0,       0D0,      0D0,     -170D0,     0D0,     1D0,
     :       1338D0,       0D0,     -5D0,      -39D0,     0D0,     0D0,
     :        716D0,       0D0,     -2D0,     -389D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 91,100 ) /
     :       1282D0,       0D0,     -3D0,      -23D0,     0D0,     1D0,
     :        742D0,       0D0,      1D0,     -391D0,     0D0,     0D0,
     :       1020D0,       0D0,    -25D0,     -495D0,     0D0,   -10D0,
     :        715D0,       0D0,     -4D0,     -326D0,     0D0,     2D0,
     :       -666D0,       0D0,     -3D0,      369D0,     0D0,    -1D0,
     :       -667D0,       0D0,      1D0,      346D0,     0D0,     1D0,
     :       -704D0,       0D0,      0D0,      304D0,     0D0,     0D0,
     :       -694D0,       0D0,      5D0,      294D0,     0D0,     2D0,
     :      -1014D0,       0D0,     -1D0,        4D0,     0D0,    -1D0,
     :       -585D0,       0D0,     -2D0,      316D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=101,110 ) /
     :       -949D0,       0D0,      1D0,        8D0,     0D0,    -1D0,
     :       -595D0,       0D0,      0D0,      258D0,     0D0,     0D0,
     :        528D0,       0D0,      0D0,     -279D0,     0D0,     0D0,
     :       -590D0,       0D0,      4D0,      252D0,     0D0,     2D0,
     :        570D0,       0D0,     -2D0,     -244D0,     0D0,    -1D0,
     :       -502D0,       0D0,      3D0,      250D0,     0D0,     2D0,
     :       -875D0,       0D0,      1D0,       29D0,     0D0,     0D0,
     :       -492D0,       0D0,     -3D0,      275D0,     0D0,    -1D0,
     :        535D0,       0D0,     -2D0,     -228D0,     0D0,    -1D0,
     :       -467D0,       0D0,      1D0,      240D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=111,120 ) /
     :        591D0,       0D0,      0D0,     -253D0,     0D0,     0D0,
     :       -453D0,       0D0,     -1D0,      244D0,     0D0,    -1D0,
     :        766D0,       0D0,      1D0,        9D0,     0D0,     0D0,
     :       -446D0,       0D0,      2D0,      225D0,     0D0,     1D0,
     :       -488D0,       0D0,      2D0,      207D0,     0D0,     1D0,
     :       -468D0,       0D0,      0D0,      201D0,     0D0,     0D0,
     :       -421D0,       0D0,      1D0,      216D0,     0D0,     1D0,
     :        463D0,       0D0,      0D0,     -200D0,     0D0,     0D0,
     :       -673D0,       0D0,      2D0,       14D0,     0D0,     0D0,
     :        658D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=121,130 ) /
     :       -438D0,       0D0,      0D0,      188D0,     0D0,     0D0,
     :       -390D0,       0D0,      0D0,      205D0,     0D0,     0D0,
     :        639D0,     -11D0,     -2D0,      -19D0,     0D0,     0D0,
     :        412D0,       0D0,     -2D0,     -176D0,     0D0,    -1D0,
     :       -361D0,       0D0,      0D0,      189D0,     0D0,     0D0,
     :        360D0,       0D0,     -1D0,     -185D0,     0D0,    -1D0,
     :        588D0,       0D0,     -3D0,      -24D0,     0D0,     0D0,
     :       -578D0,       0D0,      1D0,        5D0,     0D0,     0D0,
     :       -396D0,       0D0,      0D0,      171D0,     0D0,     0D0,
     :        565D0,       0D0,     -1D0,       -6D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=131,140 ) /
     :       -335D0,       0D0,     -1D0,      184D0,     0D0,    -1D0,
     :        357D0,       0D0,      1D0,     -154D0,     0D0,     0D0,
     :        321D0,       0D0,      1D0,     -174D0,     0D0,     0D0,
     :       -301D0,       0D0,     -1D0,      162D0,     0D0,     0D0,
     :       -334D0,       0D0,      0D0,      144D0,     0D0,     0D0,
     :        493D0,       0D0,     -2D0,      -15D0,     0D0,     0D0,
     :        494D0,       0D0,     -2D0,      -19D0,     0D0,     0D0,
     :        337D0,       0D0,     -1D0,     -143D0,     0D0,    -1D0,
     :        280D0,       0D0,     -1D0,     -144D0,     0D0,     0D0,
     :        309D0,       0D0,      1D0,     -134D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=141,150 ) /
     :       -263D0,       0D0,      2D0,      131D0,     0D0,     1D0,
     :        253D0,       0D0,      1D0,     -138D0,     0D0,     0D0,
     :        245D0,       0D0,      0D0,     -128D0,     0D0,     0D0,
     :        416D0,       0D0,     -2D0,      -17D0,     0D0,     0D0,
     :       -229D0,       0D0,      0D0,      128D0,     0D0,     0D0,
     :        231D0,       0D0,      0D0,     -120D0,     0D0,     0D0,
     :       -259D0,       0D0,      2D0,      109D0,     0D0,     1D0,
     :        375D0,       0D0,     -1D0,       -8D0,     0D0,     0D0,
     :        252D0,       0D0,      0D0,     -108D0,     0D0,     0D0,
     :       -245D0,       0D0,      1D0,      104D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=151,160 ) /
     :        243D0,       0D0,     -1D0,     -104D0,     0D0,     0D0,
     :        208D0,       0D0,      1D0,     -112D0,     0D0,     0D0,
     :        199D0,       0D0,      0D0,     -102D0,     0D0,     0D0,
     :       -208D0,       0D0,      1D0,      105D0,     0D0,     0D0,
     :        335D0,       0D0,     -2D0,      -14D0,     0D0,     0D0,
     :       -325D0,       0D0,      1D0,        7D0,     0D0,     0D0,
     :       -187D0,       0D0,      0D0,       96D0,     0D0,     0D0,
     :        197D0,       0D0,     -1D0,     -100D0,     0D0,     0D0,
     :       -192D0,       0D0,      2D0,       94D0,     0D0,     1D0,
     :       -188D0,       0D0,      0D0,       83D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=161,170 ) /
     :        276D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -286D0,       0D0,      1D0,        6D0,     0D0,     0D0,
     :        186D0,       0D0,     -1D0,      -79D0,     0D0,     0D0,
     :       -219D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :        276D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :       -153D0,       0D0,     -1D0,       84D0,     0D0,     0D0,
     :       -156D0,       0D0,      0D0,       81D0,     0D0,     0D0,
     :       -154D0,       0D0,      1D0,       78D0,     0D0,     0D0,
     :       -174D0,       0D0,      1D0,       75D0,     0D0,     0D0,
     :       -163D0,       0D0,      2D0,       69D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=171,180 ) /
     :       -228D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         91D0,       0D0,     -4D0,      -54D0,     0D0,    -2D0,
     :        175D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -159D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        141D0,       0D0,      0D0,      -72D0,     0D0,     0D0,
     :        147D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -132D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        159D0,       0D0,    -28D0,      -54D0,     0D0,    11D0,
     :        213D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        123D0,       0D0,      0D0,      -64D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=181,190 ) /
     :       -118D0,       0D0,     -1D0,       66D0,     0D0,     0D0,
     :        144D0,       0D0,     -1D0,      -61D0,     0D0,     0D0,
     :       -121D0,       0D0,      1D0,       60D0,     0D0,     0D0,
     :       -134D0,       0D0,      1D0,       56D0,     0D0,     1D0,
     :       -105D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       56D0,     0D0,     0D0,
     :        120D0,       0D0,      0D0,      -52D0,     0D0,     0D0,
     :        101D0,       0D0,      0D0,      -54D0,     0D0,     0D0,
     :       -113D0,       0D0,      0D0,       59D0,     0D0,     0D0,
     :       -106D0,       0D0,      0D0,       61D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=191,200 ) /
     :       -129D0,       0D0,      1D0,       55D0,     0D0,     0D0,
     :       -114D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :        113D0,       0D0,     -1D0,      -49D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,       51D0,     0D0,     0D0,
     :       -100D0,       0D0,     -1D0,       56D0,     0D0,     0D0,
     :         87D0,       0D0,      0D0,      -47D0,     0D0,     0D0,
     :        161D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         96D0,       0D0,      0D0,      -50D0,     0D0,     0D0,
     :        151D0,       0D0,     -1D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=201,210 ) /
     :       -104D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :       -110D0,       0D0,      0D0,       48D0,     0D0,     0D0,
     :       -100D0,       0D0,      1D0,       50D0,     0D0,     0D0,
     :         92D0,       0D0,     -5D0,       12D0,     0D0,    -2D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :        -78D0,       0D0,      0D0,       41D0,     0D0,     0D0,
     :        -77D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :          2D0,       0D0,      0D0,       54D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,      -40D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=211,220 ) /
     :        -93D0,       0D0,      0D0,       40D0,     0D0,     0D0,
     :        -83D0,       0D0,     10D0,       40D0,     0D0,    -2D0,
     :         83D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -91D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :        128D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -79D0,       0D0,      0D0,       34D0,     0D0,     0D0,
     :        -83D0,       0D0,      0D0,       47D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -44D0,     0D0,     0D0,
     :         83D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :         91D0,       0D0,      0D0,      -39D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=221,230 ) /
     :        -77D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         68D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :         71D0,       0D0,      0D0,      -31D0,     0D0,     0D0,
     :         62D0,       0D0,      0D0,      -34D0,     0D0,     0D0,
     :        -63D0,       0D0,      0D0,       33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=231,240 ) /
     :        -73D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :        115D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -103D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         63D0,       0D0,      0D0,      -28D0,     0D0,     0D0,
     :         74D0,       0D0,      0D0,      -32D0,     0D0,     0D0,
     :       -103D0,       0D0,     -3D0,        3D0,     0D0,    -1D0,
     :        -69D0,       0D0,      0D0,       30D0,     0D0,     0D0,
     :         57D0,       0D0,      0D0,      -29D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         64D0,       0D0,      0D0,      -33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=241,250 ) /
     :        -63D0,       0D0,      0D0,       26D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       23D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -24D0,     0D0,     0D0,
     :        -48D0,       0D0,      0D0,       25D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :         56D0,       0D0,      0D0,      -25D0,     0D0,     0D0,
     :         88D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :        -75D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=251,260 ) /
     :         85D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         49D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :        -74D0,       0D0,     -3D0,       -1D0,     0D0,    -1D0,
     :        -39D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         51D0,       0D0,      0D0,      -22D0,     0D0,     0D0,
     :        -40D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         41D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :        -42D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,       22D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=261,270 ) /
     :        -42D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         39D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :         46D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -53D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         81D0,       0D0,     -1D0,       -4D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         53D0,       0D0,      0D0,      -23D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=271,280 ) /
     :        -33D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -33D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :        -60D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         48D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,      -14D0,     0D0,     0D0,
     :         38D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=281,290 ) /
     :        -29D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :         44D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=291,300 ) /
     :        -60D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         36D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -37D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=301,310 ) /
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         65D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :        -30D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -31D0,       0D0,      0D0,       13D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=311,320 ) /
     :         49D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       18D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=321,330 ) /
     :        -21D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=331,340 ) /
     :         28D0,       0D0,      0D0,        0D0,     0D0,    -2D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=341,350 ) /
     :        -31D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :         20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=351,360 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         19D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        7D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=361,370 ) /
     :         13D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=371,380 ) /
     :        -26D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=381,390 ) /
     :         22D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         25D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=391,400 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=401,410 ) /
     :         15D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -25D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=411,420 ) /
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         21D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -8D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=421,430 ) /
     :         19D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=431,440 ) /
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :        -24D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=441,450 ) /
     :        -24D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :          0D0,       0D0,  -1988D0,        0D0,     0D0, -1679D0,
     :          0D0,       0D0,    -63D0,        0D0,     0D0,   -27D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     4D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          0D0,       0D0,    364D0,        0D0,     0D0,   176D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=451,460 ) /
     :          0D0,       0D0,  -1044D0,        0D0,     0D0,  -891D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    330D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=461,470 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=471,480 ) /
     :         -5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          0D0,       0D0,    -12D0,        0D0,     0D0,   -10D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=481,490 ) /
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=491,500 ) /
     :         -8D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=501,510 ) /
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=511,520 ) /
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=521,530 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=531,540 ) /
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         10D0,       0D0,     13D0,        6D0,     0D0,    -5D0,
     :          0D0,       0D0,     30D0,        0D0,     0D0,    14D0,
     :          0D0,       0D0,   -162D0,        0D0,     0D0,  -138D0,
     :          0D0,       0D0,     75D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=541,550 ) /
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=551,560 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,     -3D0,        3D0,     0D0,     1D0,
     :          0D0,       0D0,     -3D0,        0D0,     0D0,    -2D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=561,570 ) /
     :         -1D0,       0D0,      3D0,        3D0,     0D0,    -1D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    -13D0,        0D0,     0D0,   -11D0,
     :          3D0,       0D0,      6D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=571,580 ) /
     :          8D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=581,590 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=591,600 ) /
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          0D0,       0D0,    -26D0,        0D0,     0D0,   -11D0,
     :          0D0,       0D0,    -10D0,        0D0,     0D0,    -5D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=601,610 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=611,620 ) /
     :         13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=621,630 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,     -5D0,        0D0,     0D0,    -2D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=631,640 ) /
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=641,650 ) /
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=651,660 ) /
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=661,670 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=671,678 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /

*
*  Planetary argument multipliers
*    :         L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre

      DATA ( ( NAPL(I,J), I=1,14 ), J=  1, 10 ) /
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  2,  2,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  6, -3,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 11, 20 ) /
     :         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 21, 30 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  2,
     :         2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 31, 40 ) /
     :        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  2,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0, 18,-17,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 41, 50 ) /
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -4,  5,  0,  0,  0,
     :        -2,  0,  0,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -4,  3,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :        -1,  0,  1,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 51, 60 ) /
     :        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2, -2,  0,  0,  0,
     :        -2,  0,  2,  0,  2,  0,  0, -5,  9,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,
     :        -1,  0,  0,  1,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  2,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 61, 70 ) /
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 71, 80 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 81, 90 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -5,  5,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 91,100 ) /
     :        -2,  0,  1,  1,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -1, -5,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  1,  1,  1,  0,-20, 20,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=101,110 ) /
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  1,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=111,120 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -9, 17,  0,  0,  0,  0,  2,
     :         1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=121,130 ) /
     :         0,  0, -1,  1,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  1,  0, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=131,140 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0, -3,  7,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=141,150 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=151,160 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1,
     :        -1,  0,  0,  1,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=161,170 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,
     :         0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=171,180 ) /
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=181,190 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -7, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=191,200 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=201,210 ) /
     :         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=211,220 ) /
     :         0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=221,230 ) /
     :         0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=231,240 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 11,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=241,250 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=251,260 ) /
     :         0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  2,
     :         0,  0, -1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=261,270 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=271,280 ) /
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  2,
     :         0,  0, -2,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=281,290 ) /
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=291,300 ) /
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0, -2,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=301,310 ) /
     :         0,  0, -1,  1,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=311,320 ) /
     :        -2,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -7, 12,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=321,330 ) /
     :         0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=331,340 ) /
     :         0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=341,350 ) /
     :         0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=351,360 ) /
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=361,370 ) /
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=371,380 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 14,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=381,390 ) /
     :         0,  0,  0,  0,  0,  0,  0, -3,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=391,400 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=401,410 ) /
     :         0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=411,420 ) /
     :         0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=421,430 ) /
     :         0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=431,440 ) /
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=441,450 ) /
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=451,460 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=461,470 ) /
     :         0,  0,  0,  0,  0,  0,  0, -5, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=471,480 ) /
     :         0,  0,  0,  0,  0,  0, -3,  9, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=481,490 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=491,500 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=501,510 ) /
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=511,520 ) /
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=521,530 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=531,540 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=541,550 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 16,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5, 16, -4, -5,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=551,560 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=561,570 ) /
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=571,580 ) /
     :         0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=581,590 ) /
     :         0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=591,600 ) /
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -7,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=601,610 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=611,620 ) /
     :         0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=621,630 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,
     :         1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=631,640 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=641,650 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=651,660 ) /
     :         0,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0, 10, -3,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=661,670 ) /
     :         0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=671,680 ) /
     :         0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=681,687 ) /
     :         1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /

*
*  Planetary nutation coefficients, unit 1e-7 arcsec
*  longitude (sin, cos), obliquity (sin, cos)
*

      DATA ( ( ICPL(I,J), I=1,4 ), J=  1, 10 ) /
     :       1440,          0,          0,          0,
     :         56,       -117,        -42,        -40,
     :        125,        -43,          0,        -54,
     :          0,          5,          0,          0,
     :          3,         -7,         -3,          0,
     :          3,          0,          0,         -2,
     :       -114,          0,          0,         61,
     :       -219,         89,          0,          0,
     :         -3,          0,          0,          0,
     :       -462,       1604,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 11, 20 ) /
     :         99,          0,          0,        -53,
     :         -3,          0,          0,          2,
     :          0,          6,          2,          0,
     :          3,          0,          0,          0,
     :        -12,          0,          0,          0,
     :         14,       -218,        117,          8,
     :         31,       -481,       -257,        -17,
     :       -491,        128,          0,          0,
     :      -3084,       5123,       2735,       1647,
     :      -1444,       2409,      -1286,       -771 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 21, 30 ) /
     :         11,        -24,        -11,         -9,
     :         26,         -9,          0,          0,
     :        103,        -60,          0,          0,
     :          0,        -13,         -7,          0,
     :        -26,        -29,        -16,         14,
     :          9,        -27,        -14,         -5,
     :         12,          0,          0,         -6,
     :         -7,          0,          0,          0,
     :          0,         24,          0,          0,
     :        284,          0,          0,       -151 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 31, 40 ) /
     :        226,        101,          0,          0,
     :          0,         -8,         -2,          0,
     :          0,         -6,         -3,          0,
     :          5,          0,          0,         -3,
     :        -41,        175,         76,         17,
     :          0,         15,          6,          0,
     :        425,        212,       -133,        269,
     :       1200,        598,        319,       -641,
     :        235,        334,          0,          0,
     :         11,        -12,         -7,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 41, 50 ) /
     :          5,         -6,          3,          3,
     :         -5,          0,          0,          3,
     :          6,          0,          0,         -3,
     :         15,          0,          0,          0,
     :         13,          0,          0,         -7,
     :         -6,         -9,          0,          0,
     :        266,        -78,          0,          0,
     :       -460,       -435,       -232,        246,
     :          0,         15,          7,          0,
     :         -3,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 51, 60 ) /
     :          0,        131,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :        -17,        -19,        -10,          9,
     :         -9,        -11,          6,         -5,
     :         -6,          0,          0,          3,
     :        -16,          8,          0,          0,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 61, 70 ) /
     :         11,         24,         11,         -5,
     :         -3,         -4,         -2,          1,
     :          3,          0,          0,         -1,
     :          0,         -8,         -4,          0,
     :          0,          3,          0,          0,
     :          0,          5,          0,          0,
     :          0,          3,          2,          0,
     :         -6,          4,          2,          3,
     :         -3,         -5,          0,          0,
     :         -5,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 71, 80 ) /
     :          4,         24,         13,         -2,
     :        -42,         20,          0,          0,
     :        -10,        233,          0,          0,
     :         -3,          0,          0,          1,
     :         78,        -18,          0,          0,
     :          0,          3,          1,          0,
     :          0,         -3,         -1,          0,
     :          0,         -4,         -2,          1,
     :          0,         -8,         -4,         -1,
     :          0,         -5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 81, 90 ) /
     :         -7,          0,          0,          3,
     :        -14,          8,          3,          6,
     :          0,          8,         -4,          0,
     :          0,         19,         10,          0,
     :         45,        -22,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,          0,          0,
     :          0,          3,          0,          0,
     :          3,          5,          3,         -2,
     :         89,        -16,         -9,        -48 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 91,100 ) /
     :          0,          3,          0,          0,
     :         -3,          7,          4,          2,
     :       -349,        -62,          0,          0,
     :        -15,         22,          0,          0,
     :         -3,          0,          0,          0,
     :        -53,          0,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -8,          0,          0,
     :         15,         -7,         -4,         -8,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=101,110 ) /
     :        -21,        -78,          0,          0,
     :         20,        -70,        -37,        -11,
     :          0,          6,          3,          0,
     :          5,          3,          2,         -2,
     :        -17,         -4,         -2,          9,
     :          0,          6,          3,          0,
     :         32,         15,         -8,         17,
     :        174,         84,         45,        -93,
     :         11,         56,          0,          0,
     :        -66,        -12,         -6,         35 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=111,120 ) /
     :         47,          8,          4,        -25,
     :          0,          8,          4,          0,
     :         10,        -22,        -12,         -5,
     :         -3,          0,          0,          2,
     :        -24,         12,          0,          0,
     :          5,         -6,          0,          0,
     :          3,          0,          0,         -2,
     :          4,          3,          1,         -2,
     :          0,         29,         15,          0,
     :         -5,         -4,         -2,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=121,130 ) /
     :          8,         -3,         -1,         -5,
     :          0,         -3,          0,          0,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -5,          0,          0,          3,
     :         46,         66,         35,        -25,
     :        -14,          7,          0,          0,
     :          0,          3,          2,          0,
     :         -5,          0,          0,          0,
     :        -68,        -34,        -18,         36 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=131,140 ) /
     :          0,         14,          7,          0,
     :         10,         -6,         -3,         -5,
     :         -5,         -4,         -2,          3,
     :         -3,          5,          2,          1,
     :         76,         17,          9,        -41,
     :         84,        298,        159,        -45,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          2,
     :         -3,          0,          0,          1,
     :        -82,        292,        156,         44 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=141,150 ) /
     :        -73,         17,          9,         39,
     :         -9,        -16,          0,          0,
     :          3,          0,         -1,         -2,
     :         -3,          0,          0,          0,
     :         -9,         -5,         -3,          5,
     :       -439,          0,          0,          0,
     :         57,        -28,        -15,        -30,
     :          0,         -6,         -3,          0,
     :         -4,          0,          0,          2,
     :        -40,         57,         30,         21 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=151,160 ) /
     :         23,          7,          3,        -13,
     :        273,         80,         43,       -146,
     :       -449,        430,          0,          0,
     :         -8,        -47,        -25,          4,
     :          6,         47,         25,         -3,
     :          0,         23,         13,          0,
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -2,
     :        -48,       -110,        -59,         26,
     :         51,        114,         61,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=161,170 ) /
     :       -133,          0,          0,         57,
     :          0,          4,          0,          0,
     :        -21,         -6,         -3,         11,
     :          0,         -3,         -1,          0,
     :        -11,        -21,        -11,          6,
     :        -18,       -436,       -233,          9,
     :         35,         -7,          0,          0,
     :          0,          5,          3,          0,
     :         11,         -3,         -1,         -6,
     :         -5,         -3,         -1,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=171,180 ) /
     :        -53,         -9,         -5,         28,
     :          0,          3,          2,          1,
     :          4,          0,          0,         -2,
     :          0,         -4,          0,          0,
     :        -50,        194,        103,         27,
     :        -13,         52,         28,          7,
     :        -91,        248,          0,          0,
     :          6,         49,         26,         -3,
     :         -6,        -47,        -25,          3,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=181,190 ) /
     :         52,         23,         10,        -23,
     :         -3,          0,          0,          1,
     :          0,          5,          3,          0,
     :         -4,          0,          0,          0,
     :         -4,          8,          3,          2,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :          0,          8,          4,          0,
     :          0,          8,          4,          1,
     :         -4,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=191,200 ) /
     :         -4,          0,          0,          0,
     :         -8,          4,          2,          4,
     :          8,         -4,         -2,         -4,
     :          0,         15,          7,          0,
     :       -138,          0,          0,          0,
     :          0,         -7,         -3,          0,
     :          0,         -7,         -3,          0,
     :         54,          0,          0,        -29,
     :          0,         10,          4,          0,
     :         -7,          0,          0,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=201,210 ) /
     :        -37,         35,         19,         20,
     :          0,          4,          0,          0,
     :         -4,          9,          0,          0,
     :          8,          0,          0,         -4,
     :         -9,        -14,         -8,          5,
     :         -3,         -9,         -5,          3,
     :       -145,         47,          0,          0,
     :        -10,         40,         21,          5,
     :         11,        -49,        -26,         -7,
     :      -2150,          0,          0,        932 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=211,220 ) /
     :        -12,          0,          0,          5,
     :         85,          0,          0,        -37,
     :          4,          0,          0,         -2,
     :          3,          0,          0,         -2,
     :        -86,        153,          0,          0,
     :         -6,          9,          5,          3,
     :          9,        -13,         -7,         -5,
     :         -8,         12,          6,          4,
     :        -51,          0,          0,         22,
     :        -11,       -268,       -116,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=221,230 ) /
     :          0,         12,          5,          0,
     :          0,          7,          3,          0,
     :         31,          6,          3,        -17,
     :        140,         27,         14,        -75,
     :         57,         11,          6,        -30,
     :        -14,        -39,          0,          0,
     :          0,         -6,         -2,          0,
     :          4,         15,          8,         -2,
     :          0,          4,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=231,240 ) /
     :          0,         11,          5,          0,
     :          9,          6,          0,          0,
     :         -4,         10,          4,          2,
     :          5,          3,          0,          0,
     :         16,          0,          0,         -9,
     :         -3,          0,          0,          0,
     :          0,          3,          2,         -1,
     :          7,          0,          0,         -3,
     :        -25,         22,          0,          0,
     :         42,        223,        119,        -22 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=241,250 ) /
     :        -27,       -143,        -77,         14,
     :          9,         49,         26,         -5,
     :      -1166,          0,          0,        505,
     :         -5,          0,          0,          2,
     :         -6,          0,          0,          3,
     :         -8,          0,          1,          4,
     :          0,         -4,          0,          0,
     :        117,          0,          0,        -63,
     :         -4,          8,          4,          2,
     :          3,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=251,260 ) /
     :         -5,          0,          0,          2,
     :          0,         31,          0,          0,
     :         -5,          0,          1,          3,
     :          4,          0,          0,         -2,
     :         -4,          0,          0,          2,
     :        -24,        -13,         -6,         10,
     :          3,          0,          0,          0,
     :          0,        -32,        -17,          0,
     :          8,         12,          5,         -3,
     :          3,          0,          0,         -1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=261,270 ) /
     :          7,         13,          0,          0,
     :         -3,         16,          0,          0,
     :         50,          0,          0,        -27,
     :          0,         -5,         -3,          0,
     :         13,          0,          0,          0,
     :          0,          5,          3,          1,
     :         24,          5,          2,        -11,
     :          5,        -11,         -5,         -2,
     :         30,         -3,         -2,        -16,
     :         18,          0,          0,         -9 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=271,280 ) /
     :          8,        614,          0,          0,
     :          3,         -3,         -1,         -2,
     :          6,         17,          9,         -3,
     :         -3,         -9,         -5,          2,
     :          0,          6,          3,         -1,
     :       -127,         21,          9,         55,
     :          3,          5,          0,          0,
     :         -6,        -10,         -4,          3,
     :          5,          0,          0,          0,
     :         16,          9,          4,         -7 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=281,290 ) /
     :          3,          0,          0,         -2,
     :          0,         22,          0,          0,
     :          0,         19,         10,          0,
     :          7,          0,          0,         -4,
     :          0,         -5,         -2,          0,
     :          0,          3,          1,          0,
     :         -9,          3,          1,          4,
     :         17,          0,          0,         -7,
     :          0,         -3,         -2,         -1,
     :        -20,         34,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=291,300 ) /
     :        -10,          0,          1,          5,
     :         -4,          0,          0,          2,
     :         22,        -87,          0,          0,
     :         -4,          0,          0,          2,
     :         -3,         -6,         -2,          1,
     :        -16,         -3,         -1,          7,
     :          0,         -3,         -2,          0,
     :          4,          0,          0,          0,
     :        -68,         39,          0,          0,
     :         27,          0,          0,        -14 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=301,310 ) /
     :          0,         -4,          0,          0,
     :        -25,          0,          0,          0,
     :        -12,         -3,         -2,          6,
     :          3,          0,          0,         -1,
     :          3,         66,         29,         -1,
     :        490,          0,          0,       -213,
     :        -22,         93,         49,         12,
     :         -7,         28,         15,          4,
     :         -3,         13,          7,          2,
     :        -46,         14,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=311,320 ) /
     :         -5,          0,          0,          0,
     :          2,          1,          0,          0,
     :          0,         -3,          0,          0,
     :        -28,          0,          0,         15,
     :          5,          0,          0,         -2,
     :          0,          3,          0,          0,
     :        -11,          0,          0,          5,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          1,
     :         25,        106,         57,        -13 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=321,330 ) /
     :          5,         21,         11,         -3,
     :       1485,          0,          0,          0,
     :         -7,        -32,        -17,          4,
     :          0,          5,          3,          0,
     :         -6,         -3,         -2,          3,
     :         30,         -6,         -2,        -13,
     :         -4,          4,          0,          0,
     :        -19,          0,          0,         10,
     :          0,          4,          2,         -1,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=331,340 ) /
     :          4,          0,          0,         -2,
     :          0,         -3,         -1,          0,
     :         -3,          0,          0,          0,
     :          5,          3,          1,         -2,
     :          0,         11,          0,          0,
     :        118,          0,          0,        -52,
     :          0,         -5,         -3,          0,
     :        -28,         36,          0,          0,
     :          5,         -5,          0,          0,
     :         14,        -59,        -31,         -8 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=341,350 ) /
     :          0,          9,          5,          1,
     :       -458,          0,          0,        198,
     :          0,        -45,        -20,          0,
     :          9,          0,          0,         -5,
     :          0,         -3,          0,          0,
     :          0,         -4,         -2,         -1,
     :         11,          0,          0,         -6,
     :          6,          0,          0,         -2,
     :        -16,         23,          0,          0,
     :          0,         -4,         -2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=351,360 ) /
     :         -5,          0,          0,          2,
     :       -166,        269,          0,          0,
     :         15,          0,          0,         -8,
     :         10,          0,          0,         -4,
     :        -78,         45,          0,          0,
     :          0,         -5,         -2,          0,
     :          7,          0,          0,         -4,
     :         -5,        328,          0,          0,
     :          3,          0,          0,         -2,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=361,370 ) /
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -4,         -2,          0,
     :      -1223,        -26,          0,          0,
     :          0,          7,          3,          0,
     :          3,          0,          0,          0,
     :          0,          3,          2,          0,
     :         -6,         20,          0,          0,
     :       -368,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=371,380 ) /
     :        -75,          0,          0,          0,
     :         11,          0,          0,         -6,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :        -13,        -30,          0,          0,
     :         21,          3,          0,          0,
     :         -3,          0,          0,          1,
     :         -4,          0,          0,          2,
     :          8,        -27,          0,          0,
     :        -19,        -11,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=381,390 ) /
     :         -4,          0,          0,          2,
     :          0,          5,          2,          0,
     :         -6,          0,          0,          2,
     :         -8,          0,          0,          0,
     :         -1,          0,          0,          0,
     :        -14,          0,          0,          6,
     :          6,          0,          0,          0,
     :        -74,          0,          0,         32,
     :          0,         -3,         -1,          0,
     :          4,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=391,400 ) /
     :          8,         11,          0,          0,
     :          0,          3,          2,          0,
     :       -262,          0,          0,        114,
     :          0,         -4,          0,          0,
     :         -7,          0,          0,          4,
     :          0,        -27,        -12,          0,
     :        -19,         -8,         -4,          8,
     :        202,          0,          0,        -87,
     :         -8,         35,         19,          5,
     :          0,          4,          2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=401,410 ) /
     :         16,         -5,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -3,          0,          0,
     :          1,          0,          0,          0,
     :        -35,        -48,        -21,         15,
     :         -3,         -5,         -2,          1,
     :          6,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          0,         -5,          0,          0,
     :         12,         55,         29,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=411,420 ) /
     :          0,          5,          3,          0,
     :       -598,          0,          0,          0,
     :         -3,        -13,         -7,          1,
     :         -5,         -7,         -3,          2,
     :          3,          0,          0,         -1,
     :          5,         -7,          0,          0,
     :          4,          0,          0,         -2,
     :         16,         -6,          0,          0,
     :          8,         -3,          0,          0,
     :          8,        -31,        -16,         -4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=421,430 ) /
     :          0,          3,          1,          0,
     :        113,          0,          0,        -49,
     :          0,        -24,        -10,          0,
     :          4,          0,          0,         -2,
     :         27,          0,          0,          0,
     :         -3,          0,          0,          1,
     :          0,         -4,         -2,          0,
     :          5,          0,          0,         -2,
     :          0,         -3,          0,          0,
     :        -13,          0,          0,          6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=431,440 ) /
     :          5,          0,          0,         -2,
     :        -18,        -10,         -4,          8,
     :         -4,        -28,          0,          0,
     :         -5,          6,          3,          2,
     :         -3,          0,          0,          1,
     :         -5,         -9,         -4,          2,
     :         17,          0,          0,         -7,
     :         11,          4,          0,          0,
     :          0,         -6,         -2,          0,
     :         83,         15,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=441,450 ) /
     :         -4,          0,          0,          2,
     :          0,       -114,        -49,          0,
     :        117,          0,          0,        -51,
     :         -5,         19,         10,          2,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          2,
     :          0,         -3,         -1,          0,
     :          3,          0,          0,          0,
     :          0,         -6,         -2,          0,
     :        393,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=451,460 ) /
     :         -4,         21,         11,          2,
     :         -6,          0,         -1,          3,
     :         -3,          8,          4,          1,
     :          8,          0,          0,          0,
     :         18,        -29,        -13,         -8,
     :          8,         34,         18,         -4,
     :         89,          0,          0,          0,
     :          3,         12,          6,         -1,
     :         54,        -15,         -7,        -24,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=461,470 ) /
     :          3,          0,          0,         -1,
     :          0,         35,          0,          0,
     :       -154,        -30,        -13,         67,
     :         15,          0,          0,          0,
     :          0,          4,          2,          0,
     :          0,          9,          0,          0,
     :         80,        -71,        -31,        -35,
     :          0,        -20,         -9,          0,
     :         11,          5,          2,         -5,
     :         61,        -96,        -42,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=471,480 ) /
     :         14,          9,          4,         -6,
     :        -11,         -6,         -3,          5,
     :          0,         -3,         -1,          0,
     :        123,       -415,       -180,        -53,
     :          0,          0,          0,        -35,
     :         -5,          0,          0,          0,
     :          7,        -32,        -17,         -4,
     :          0,         -9,         -5,          0,
     :          0,         -4,          2,          0,
     :        -89,          0,          0,         38 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=481,490 ) /
     :          0,        -86,        -19,         -6,
     :          0,          0,        -19,          6,
     :       -123,       -416,       -180,         53,
     :          0,         -3,         -1,          0,
     :         12,         -6,         -3,         -5,
     :        -13,          9,          4,          6,
     :          0,        -15,         -7,          0,
     :          3,          0,          0,         -1,
     :        -62,        -97,        -42,         27,
     :        -11,          5,          2,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=491,500 ) /
     :          0,        -19,         -8,          0,
     :         -3,          0,          0,          1,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :        -85,        -70,        -31,         37,
     :        163,        -12,         -5,        -72,
     :        -63,        -16,         -7,         28,
     :        -21,        -32,        -14,          9,
     :          0,         -3,         -1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=501,510 ) /
     :          3,          0,          0,         -2,
     :          0,          8,          0,          0,
     :          3,         10,          4,         -1,
     :          3,          0,          0,         -1,
     :          0,         -7,         -3,          0,
     :          0,         -4,         -2,          0,
     :          6,         19,          0,          0,
     :          5,       -173,        -75,         -2,
     :          0,         -7,         -3,          0,
     :          7,        -12,         -5,         -3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=511,520 ) /
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -1,
     :         74,          0,          0,        -32,
     :         -3,         12,          6,          2,
     :         26,        -14,         -6,        -11,
     :         19,          0,          0,         -8,
     :          6,         24,         13,         -3,
     :         83,          0,          0,          0,
     :          0,        -10,         -5,          0,
     :         11,         -3,         -1,         -5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=521,530 ) /
     :          3,          0,          1,         -1,
     :          3,          0,          0,         -1,
     :         -4,          0,          0,          0,
     :          5,        -23,        -12,         -3,
     :       -339,          0,          0,        147,
     :          0,        -10,         -5,          0,
     :          5,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          0,         -4,         -2,          0,
     :         18,         -3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=531,540 ) /
     :          9,        -11,         -5,         -4,
     :         -8,          0,          0,          4,
     :          3,          0,          0,         -1,
     :          0,          9,          0,          0,
     :          6,         -9,         -4,         -2,
     :         -4,        -12,          0,          0,
     :         67,        -91,        -39,        -29,
     :         30,        -18,         -8,        -13,
     :          0,          0,          0,          0,
     :          0,       -114,        -50,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=541,550 ) /
     :          0,          0,          0,         23,
     :        517,         16,          7,       -224,
     :          0,         -7,         -3,          0,
     :        143,         -3,         -1,        -62,
     :         29,          0,          0,        -13,
     :         -4,          0,          0,          2,
     :         -6,          0,          0,          3,
     :          5,         12,          5,         -2,
     :        -25,          0,          0,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=551,560 ) /
     :          0,          4,          2,          0,
     :        -22,         12,          5,         10,
     :         50,          0,          0,        -22,
     :          0,          7,          4,          0,
     :          0,          3,          1,          0,
     :         -4,          4,          2,          2,
     :         -5,        -11,         -5,          2,
     :          0,          4,          2,          0,
     :          4,         17,          9,         -2,
     :         59,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=561,570 ) /
     :          0,         -4,         -2,          0,
     :         -8,          0,          0,          4,
     :         -3,          0,          0,          0,
     :          4,        -15,         -8,         -2,
     :        370,         -8,          0,       -160,
     :          0,          0,         -3,          0,
     :          0,          3,          1,          0,
     :         -6,          3,          1,          3,
     :          0,          6,          0,          0,
     :        -10,          0,          0,          4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=571,580 ) /
     :          0,          9,          4,          0,
     :          4,         17,          7,         -2,
     :         34,          0,          0,        -15,
     :          0,          5,          3,          0,
     :         -5,          0,          0,          2,
     :        -37,         -7,         -3,         16,
     :          3,         13,          7,         -2,
     :         40,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :       -184,         -3,         -1,         80 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=581,590 ) /
     :         -3,          0,          0,          1,
     :         -3,          0,          0,          0,
     :          0,        -10,         -6,         -1,
     :         31,         -6,          0,        -13,
     :         -3,        -32,        -14,          1,
     :         -7,          0,          0,          3,
     :          0,         -8,         -4,          0,
     :          3,         -4,          0,          0,
     :          0,          4,          0,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=591,600 ) /
     :         19,        -23,        -10,          2,
     :          0,          0,          0,        -10,
     :          0,          3,          2,          0,
     :          0,          9,          5,         -1,
     :         28,          0,          0,          0,
     :          0,         -7,         -4,          0,
     :          8,         -4,          0,         -4,
     :          0,          0,         -2,          0,
     :          0,          3,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=601,610 ) /
     :         -9,          0,          1,          4,
     :          3,         12,          5,         -1,
     :         17,         -3,         -1,          0,
     :          0,          7,          4,          0,
     :         19,          0,          0,          0,
     :          0,         -5,         -3,          0,
     :         14,         -3,          0,         -1,
     :          0,          0,         -1,          0,
     :          0,          0,          0,         -5,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=611,620 ) /
     :         13,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :          2,          9,          4,          3,
     :          0,          0,          0,         -4,
     :          8,          0,          0,          0,
     :          0,          4,          2,          0,
     :          6,          0,          0,         -3,
     :          6,          0,          0,          0,
     :          0,          3,          1,          0,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=621,630 ) /
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          0,
     :          6,          0,          0,          0,
     :          7,          0,          0,          0,
     :         -4,          0,          0,          0,
     :          4,          0,          0,          0,
     :          6,          0,          0,          0,
     :          0,         -4,          0,          0,
     :          0,         -4,          0,          0,
     :          5,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=631,640 ) /
     :         -3,          0,          0,          0,
     :          4,          0,          0,          0,
     :         -5,          0,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :         13,          0,          0,          0,
     :         21,         11,          0,          0,
     :          0,         -5,          0,          0,
     :          0,         -5,         -2,          0,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=641,650 ) /
     :          0,         -5,          0,          0,
     :         -3,          0,          0,          2,
     :         20,         10,          0,          0,
     :        -34,          0,          0,          0,
     :        -19,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :         -6,          0,          0,          3,
     :         -4,          0,          0,          0,
     :          3,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=651,660 ) /
     :          3,          0,          0,          0,
     :          4,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          6,          0,          0,         -3,
     :         -8,          0,          0,          3,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :        126,        -63,        -27,        -55,
     :         -5,          0,          1,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=661,670 ) /
     :         -3,         28,         15,          2,
     :          5,          0,          1,         -2,
     :          0,          9,          4,          1,
     :          0,          9,          4,         -1,
     :       -126,        -63,        -27,         55,
     :          3,          0,          0,         -1,
     :         21,        -11,         -6,        -11,
     :          0,         -4,          0,          0,
     :        -21,        -11,         -6,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=671,680 ) /
     :          0,          3,          1,          0,
     :          8,          0,          0,         -4,
     :         -6,          0,          0,          3,
     :         -3,          0,          0,          1,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          1,
     :         -5,          0,          0,          2,
     :         24,        -12,         -5,        -11,
     :          0,          3,          1,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=681,687 ) /
     :          0,          3,          2,          0,
     :        -24,        -12,         -5,         10,
     :          4,          0,         -1,         -2,
     :         13,          0,          0,         -6,
     :          7,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          3,          0,          0,         -1 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental date J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  -------------------
*  LUNI-SOLAR NUTATION
*  -------------------

*
*  Fundamental (Delaunay) arguments
*

*  Mean anomaly of the Moon (IERS 2003).
      EL = iau_FAL03 ( T )

*  Mean anomaly of the Sun (MHB2000).
      ELP = MOD (       1287104.79305D0 +
     :            T*( 129596581.0481D0 +
     :            T*(       - 0.5532D0 +
     :            T*(         0.000136D0 +
     :            T*(       - 0.00001149D0 )))), TURNAS ) * DAS2R

*  Mean longitude of the Moon minus that of the ascending node
*  (IERS 2003.
      F = iau_FAF03 ( T )

*  Mean elongation of the Moon from the Sun (MHB2000).
      D = MOD (        1072260.70369D0 +
     :          T*( 1602961601.2090D0 +
     :          T*(        - 6.3706D0 +
     :          T*(          0.006593D0 +
     :          T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R

*  Mean longitude of the ascending node of the Moon (IERS 2003).
      OM = iau_FAOM03 ( T )

*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0

*  Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NALS(1,I) ) * EL  +
     :               DBLE ( NALS(2,I) ) * ELP +
     :               DBLE ( NALS(3,I) ) * F   +
     :               DBLE ( NALS(4,I) ) * D   +
     :               DBLE ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG

 100  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R

*  ------------------
*  PLANETARY NUTATION
*  ------------------

*  n.b.  The MHB2000 code computes the luni-solar and planetary nutation
*        in different routines, using slightly different Delaunay
*        arguments in the two cases.  This behaviour is faithfully
*        reproduced here.  Use of the IERS 2003 expressions for both
*        cases leads to negligible changes, well below
*        0.1 microarcsecond.

*  Mean anomaly of the Moon (MHB2000).
      AL = MOD ( 2.35555598D0 + 8328.6914269554D0 * T, D2PI )

*  Mean anomaly of the Sun (MHB2000).
      ALSU = MOD ( 6.24006013D0 + 628.301955D0 * T, D2PI )

*  Mean longitude of the Moon minus that of the ascending node
* (MHB2000).
      AF = MOD ( 1.627905234D0 + 8433.466158131D0 * T, D2PI )

*  Mean elongation of the Moon from the Sun (MHB2000).
      AD = MOD ( 5.198466741D0 + 7771.3771468121D0 * T, D2PI )

*  Mean longitude of the ascending node of the Moon (MHB2000).
      AOM = MOD ( 2.18243920D0 - 33.757045D0 * T, D2PI )

*  General accumulated precession in longitude (IERS 2003).
      APA = iau_FAPA03 ( T )

*  Planetary longitudes, Mercury through Uranus (IERS 2003).
      ALME = iau_FAME03 ( T )
      ALVE = iau_FAVE03 ( T )
      ALEA = iau_FAE03 ( T )
      ALMA = iau_FAMA03 ( T )
      ALJU = iau_FAJU03 ( T )
      ALSA = iau_FASA03 ( T )
      ALUR = iau_FAUR03 ( T )

*  Neptune longitude (MHB2000).
      ALNE = MOD ( 5.321159000D0 + 3.8127774000D0 * T, D2PI )

*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0

*  Summation of planetary nutation series (in reverse order).
      DO 200 I = NPL, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NAPL( 1,I) ) * AL   +
     :               DBLE ( NAPL( 2,I) ) * ALSU +
     :               DBLE ( NAPL( 3,I) ) * AF   +
     :               DBLE ( NAPL( 4,I) ) * AD   +
     :               DBLE ( NAPL( 5,I) ) * AOM  +
     :               DBLE ( NAPL( 6,I) ) * ALME +
     :               DBLE ( NAPL( 7,I) ) * ALVE +
     :               DBLE ( NAPL( 8,I) ) * ALEA +
     :               DBLE ( NAPL( 9,I) ) * ALMA +
     :               DBLE ( NAPL(10,I) ) * ALJU +
     :               DBLE ( NAPL(11,I) ) * ALSA +
     :               DBLE ( NAPL(12,I) ) * ALUR +
     :               DBLE ( NAPL(13,I) ) * ALNE +
     :               DBLE ( NAPL(14,I) ) * APA, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + DBLE( ICPL(1,I)) * SARG + DBLE( ICPL(2,I)) * CARG
         DE = DE + DBLE( ICPL(3,I)) * SARG + DBLE( ICPL(4,I)) * CARG

 200  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSIPL = DP * U2R
      DEPSPL = DE * U2R

*  -------
*  RESULTS
*  -------

*  Add luni-solar and planetary components.
      DPSI = DPSILS + DPSIPL
      DEPS = DEPSLS + DEPSPL

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      SUBROUTINE iau_BPN2XY ( RBPN, X, Y )
*+
*  - - - - - - - - - - -
*   i a u _ B P N 2 X Y
*  - - - - - - - - - - -
*
*  Extract from the bias-precession-nutation matrix the X,Y coordinates
*  of the Celestial Intermediate Pole.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     RBPN      d(3,3)    celestial-to-true matrix (Note 1)
*
*  Returned:
*     X,Y         d       Celestial Intermediate Pole (Note 2)
*
*  Notes:
*
*  1) The matrix RBPN transforms vectors from GCRS to true equator (and
*     CIO or equinox) of date, and therefore the Celestial Intermediate
*     Pole unit vector is the bottom row of the matrix.
*
*  2) X,Y are components of the Celestial Intermediate Pole unit vector
*     in the Geocentric Celestial Reference System.
*
*  Reference:
*
*     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
*     "Expressions for the Celestial Intermediate Pole and Celestial
*     Ephemeris Origin consistent with the IAU 2000A precession-nutation
*     model", Astron.Astrophys. 400, 1145-1154 (2003)
*
*     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
*          intermediate origin" (CIO) by IAU 2006 Resolution 2.
*
*  This revision:  2010 January 18
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RBPN(3,3), X, Y

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Extract the X,Y coordinates.
      X = RBPN(3,1)
      Y = RBPN(3,2)

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_S00 ( DATE1, DATE2, X, Y )
*+
*  - - - - - - - -
*   i a u _ S 0 0
*  - - - - - - - -
*
*  The CIO locator s, positioning the Celestial Intermediate Origin on
*  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
*  coordinates.  Compatible with IAU 2000A precession-nutation.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2    d      TT as a 2-part Julian Date (Note 1)
*     X,Y            d      CIP coordinates (Note 3)
*
*  Returned:
*     iau_S00        d      the CIO locator s in radians (Note 2)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The CIO locator s is the difference between the right ascensions
*     of the same point in two systems:  the two systems are the GCRS
*     and the CIP,CIO, and the point is the ascending node of the
*     CIP equator.  The quantity s remains below 0.1 arcsecond
*     throughout 1900-2100.
*
*  3) The series used to compute s is in fact for s+XY/2, where X and Y
*     are the x and y components of the CIP unit vector;  this series is
*     more compact than a direct series for s would be.  This routine
*     requires X,Y to be supplied by the caller, who is responsible for
*     providing values that are consistent with the supplied date.
*
*  4) The model is consistent with the IAU 2000A precession-nutation.
*
*  Called:
*     iau_FAL03    mean anomaly of the Moon
*     iau_FALP03   mean anomaly of the Sun
*     iau_FAF03    mean argument of the latitude of the Moon
*     iau_FAD03    mean elongation of the Moon from the Sun
*     iau_FAOM03   mean longitude of the Moon's ascending node
*     iau_FAVE03   mean longitude of Venus
*     iau_FAE03    mean longitude of Earth
*     iau_FAPA03   general accumulated precession in longitude
*
*  References:
*
*     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
*     "Expressions for the Celestial Intermediate Pole and Celestial
*     Ephemeris Origin consistent with the IAU 2000A precession-nutation
*     model", Astron.Astrophys. 400, 1145-1154 (2003)
*
*     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
*          intermediate origin" (CIO) by IAU 2006 Resolution 2.
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*  This revision:  2010 January 18
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, X, Y

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Time since J2000.0, in Julian centuries
      DOUBLE PRECISION T

*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION A, S0, S1, S2, S3, S4, S5
      DOUBLE PRECISION iau_FAL03, iau_FALP03, iau_FAF03,
     :                 iau_FAD03, iau_FAOM03, iau_FAVE03, iau_FAE03,
     :                 iau_FAPA03

*  Fundamental arguments
      DOUBLE PRECISION FA(8)

*  ---------------------
*  The series for s+XY/2
*  ---------------------

*  Number of terms in the series
      INTEGER NSP, NS0, NS1, NS2, NS3, NS4
      PARAMETER ( NSP=6, NS0=33, NS1=3, NS2=25, NS3=4, NS4=1 )

*  Polynomial coefficients
      DOUBLE PRECISION SP ( NSP )

*  Coefficients of l,l',F,D,Om,LVe,LE,pA
      INTEGER KS0 ( 8, NS0 ),
     :        KS1 ( 8, NS1 ),
     :        KS2 ( 8, NS2 ),
     :        KS3 ( 8, NS3 ),
     :        KS4 ( 8, NS4 )

*  Sine and cosine coefficients
      DOUBLE PRECISION SS0 ( 2, NS0 ),
     :                 SS1 ( 2, NS1 ),
     :                 SS2 ( 2, NS2 ),
     :                 SS3 ( 2, NS3 ),
     :                 SS4 ( 2, NS4 )

*  Polynomial coefficients
      DATA SP /    94    D-6,
     :           3808.35 D-6,
     :           -119.94 D-6,
     :         -72574.09 D-6,
     :             27.70 D-6,
     :             15.61 D-6 /

*  Argument coefficients for t^0
      DATA ( ( KS0(I,J), I=1,8), J=1,10 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,
     :  0,  0,  0,  0,  2,  0,  0,  0,
     :  0,  0,  2, -2,  3,  0,  0,  0,
     :  0,  0,  2, -2,  1,  0,  0,  0,
     :  0,  0,  2, -2,  2,  0,  0,  0,
     :  0,  0,  2,  0,  3,  0,  0,  0,
     :  0,  0,  2,  0,  1,  0,  0,  0,
     :  0,  0,  0,  0,  3,  0,  0,  0,
     :  0,  1,  0,  0,  1,  0,  0,  0,
     :  0,  1,  0,  0, -1,  0,  0,  0 /
      DATA ( ( KS0(I,J), I=1,8), J=11,20 ) /
     :  1,  0,  0,  0, -1,  0,  0,  0,
     :  1,  0,  0,  0,  1,  0,  0,  0,
     :  0,  1,  2, -2,  3,  0,  0,  0,
     :  0,  1,  2, -2,  1,  0,  0,  0,
     :  0,  0,  4, -4,  4,  0,  0,  0,
     :  0,  0,  1, -1,  1, -8, 12,  0,
     :  0,  0,  2,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  2,  0,  0,  0,
     :  1,  0,  2,  0,  3,  0,  0,  0,
     :  1,  0,  2,  0,  1,  0,  0,  0 /
      DATA ( ( KS0(I,J), I=1,8), J=21,30 ) /
     :  0,  0,  2, -2,  0,  0,  0,  0,
     :  0,  1, -2,  2, -3,  0,  0,  0,
     :  0,  1, -2,  2, -1,  0,  0,  0,
     :  0,  0,  0,  0,  0,  8,-13, -1,
     :  0,  0,  0,  2,  0,  0,  0,  0,
     :  2,  0, -2,  0, -1,  0,  0,  0,
     :  0,  1,  2, -2,  2,  0,  0,  0,
     :  1,  0,  0, -2,  1,  0,  0,  0,
     :  1,  0,  0, -2, -1,  0,  0,  0,
     :  0,  0,  4, -2,  4,  0,  0,  0 /
      DATA ( ( KS0(I,J), I=1,8), J=31,NS0 ) /
     :  0,  0,  2, -2,  4,  0,  0,  0,
     :  1,  0, -2,  0, -3,  0,  0,  0,
     :  1,  0, -2,  0, -1,  0,  0,  0 /

*  Argument coefficients for t^1
      DATA ( ( KS1(I,J), I=1,8), J=1,NS1 ) /
     :  0,  0,  0,  0,  2,  0,  0,  0,
     :  0,  0,  0,  0,  1,  0,  0,  0,
     :  0,  0,  2, -2,  3,  0,  0,  0 /

*  Argument coefficients for t^2
      DATA ( ( KS2(I,J), I=1,8), J=1,10 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,
     :  0,  0,  2, -2,  2,  0,  0,  0,
     :  0,  0,  2,  0,  2,  0,  0,  0,
     :  0,  0,  0,  0,  2,  0,  0,  0,
     :  0,  1,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  2,  0,  0,  0,
     :  0,  0,  2,  0,  1,  0,  0,  0,
     :  1,  0,  2,  0,  2,  0,  0,  0,
     :  0,  1, -2,  2, -2,  0,  0,  0 /
      DATA ( ( KS2(I,J), I=1,8), J=11,20 ) /
     :  1,  0,  0, -2,  0,  0,  0,  0,
     :  0,  0,  2, -2,  1,  0,  0,  0,
     :  1,  0, -2,  0, -2,  0,  0,  0,
     :  0,  0,  0,  2,  0,  0,  0,  0,
     :  1,  0,  0,  0,  1,  0,  0,  0,
     :  1,  0, -2, -2, -2,  0,  0,  0,
     :  1,  0,  0,  0, -1,  0,  0,  0,
     :  1,  0,  2,  0,  1,  0,  0,  0,
     :  2,  0,  0, -2,  0,  0,  0,  0,
     :  2,  0, -2,  0, -1,  0,  0,  0 /
      DATA ( ( KS2(I,J), I=1,8), J=21,NS2 ) /
     :  0,  0,  2,  2,  2,  0,  0,  0,
     :  2,  0,  2,  0,  2,  0,  0,  0,
     :  2,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2, -2,  2,  0,  0,  0,
     :  0,  0,  2,  0,  0,  0,  0,  0 /

*  Argument coefficients for t^3
      DATA ( ( KS3(I,J), I=1,8), J=1,NS3 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,
     :  0,  0,  2, -2,  2,  0,  0,  0,
     :  0,  0,  2,  0,  2,  0,  0,  0,
     :  0,  0,  0,  0,  2,  0,  0,  0 /

*  Argument coefficients for t^4
      DATA ( ( KS4(I,J), I=1,8), J=1,NS4 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0 /

*  Sine and cosine coefficients for t^0
      DATA ( ( SS0(I,J), I=1,2), J=1,10 ) /
     :            -2640.73D-6,          +0.39D-6,
     :              -63.53D-6,          +0.02D-6,
     :              -11.75D-6,          -0.01D-6,
     :              -11.21D-6,          -0.01D-6,
     :               +4.57D-6,           0.00D-6,
     :               -2.02D-6,           0.00D-6,
     :               -1.98D-6,           0.00D-6,
     :               +1.72D-6,           0.00D-6,
     :               +1.41D-6,          +0.01D-6,
     :               +1.26D-6,          +0.01D-6 /
      DATA ( ( SS0(I,J), I=1,2), J=11,20 ) /
     :               +0.63D-6,           0.00D-6,
     :               +0.63D-6,           0.00D-6,
     :               -0.46D-6,           0.00D-6,
     :               -0.45D-6,           0.00D-6,
     :               -0.36D-6,           0.00D-6,
     :               +0.24D-6,          +0.12D-6,
     :               -0.32D-6,           0.00D-6,
     :               -0.28D-6,           0.00D-6,
     :               -0.27D-6,           0.00D-6,
     :               -0.26D-6,           0.00D-6 /
      DATA ( ( SS0(I,J), I=1,2), J=21,30 ) /
     :               +0.21D-6,           0.00D-6,
     :               -0.19D-6,           0.00D-6,
     :               -0.18D-6,           0.00D-6,
     :               +0.10D-6,          -0.05D-6,
     :               -0.15D-6,           0.00D-6,
     :               +0.14D-6,           0.00D-6,
     :               +0.14D-6,           0.00D-6,
     :               -0.14D-6,           0.00D-6,
     :               -0.14D-6,           0.00D-6,
     :               -0.13D-6,           0.00D-6 /
      DATA ( ( SS0(I,J), I=1,2), J=31,NS0 ) /
     :               +0.11D-6,           0.00D-6,
     :               -0.11D-6,           0.00D-6,
     :               -0.11D-6,           0.00D-6 /

*  Sine and cosine coefficients for t^1
      DATA ( ( SS1(I,J), I=1,2), J=1,NS1 ) /
     :               -0.07D-6,          +3.57D-6,
     :               +1.71D-6,          -0.03D-6,
     :                0.00D-6,          +0.48D-6 /

*  Sine and cosine coefficients for t^2
      DATA ( ( SS2(I,J), I=1,2), J=1,10 ) /
     :             +743.53D-6,          -0.17D-6,
     :              +56.91D-6,          +0.06D-6,
     :               +9.84D-6,          -0.01D-6,
     :               -8.85D-6,          +0.01D-6,
     :               -6.38D-6,          -0.05D-6,
     :               -3.07D-6,           0.00D-6,
     :               +2.23D-6,           0.00D-6,
     :               +1.67D-6,           0.00D-6,
     :               +1.30D-6,           0.00D-6,
     :               +0.93D-6,           0.00D-6 /
      DATA ( ( SS2(I,J), I=1,2), J=11,20 ) /
     :               +0.68D-6,           0.00D-6,
     :               -0.55D-6,           0.00D-6,
     :               +0.53D-6,           0.00D-6,
     :               -0.27D-6,           0.00D-6,
     :               -0.27D-6,           0.00D-6,
     :               -0.26D-6,           0.00D-6,
     :               -0.25D-6,           0.00D-6,
     :               +0.22D-6,           0.00D-6,
     :               -0.21D-6,           0.00D-6,
     :               +0.20D-6,           0.00D-6 /
      DATA ( ( SS2(I,J), I=1,2), J=21,NS2 ) /
     :               +0.17D-6,           0.00D-6,
     :               +0.13D-6,           0.00D-6,
     :               -0.13D-6,           0.00D-6,
     :               -0.12D-6,           0.00D-6,
     :               -0.11D-6,           0.00D-6 /

*  Sine and cosine coefficients for t^3
      DATA ( ( SS3(I,J), I=1,2), J=1,NS3 ) /
     :               +0.30D-6,         -23.51D-6,
     :               -0.03D-6,          -1.39D-6,
     :               -0.01D-6,          -0.24D-6,
     :                0.00D-6,          +0.22D-6 /

*  Sine and cosine coefficients for t^4
      DATA ( ( SS4(I,J), I=1,2), J=1,NS4 ) /
     :               -0.26D-6,          -0.01D-6 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  Fundamental Arguments (from IERS Conventions 2003)

*  Mean anomaly of the Moon.
      FA(1) = iau_FAL03 ( T )

*  Mean anomaly of the Sun.
      FA(2) = iau_FALP03 ( T )

*  Mean longitude of the Moon minus that of the ascending node.
      FA(3) = iau_FAF03 ( T )

*  Mean elongation of the Moon from the Sun.
      FA(4) = iau_FAD03 ( T )

*  Mean longitude of the ascending node of the Moon.
      FA(5) = iau_FAOM03 ( T )

*  Mean longitude of Venus.
      FA(6) = iau_FAVE03 ( T )

*  Mean longitude of Earth.
      FA(7) = iau_FAE03 ( T )

*  General precession in longitude.
      FA(8) = iau_FAPA03 ( T )

*  Evaluate s.
      S0 = SP(1)
      S1 = SP(2)
      S2 = SP(3)
      S3 = SP(4)
      S4 = SP(5)
      S5 = SP(6)

      DO 2 I = NS0,1,-1
         A = 0D0
         DO 1 J=1,8
            A = A + DBLE(KS0(J,I))*FA(J)
 1       CONTINUE
         S0 = S0 + ( SS0(1,I)*SIN(A) + SS0(2,I)*COS(A) )
 2    CONTINUE

      DO 4 I = NS1,1,-1
         A = 0D0
         DO 3 J=1,8
            A = A + DBLE(KS1(J,I))*FA(J)
 3       CONTINUE
         S1 = S1 + ( SS1(1,I)*SIN(A) + SS1(2,I)*COS(A) )
 4    CONTINUE

      DO 6 I = NS2,1,-1
         A = 0D0
         DO 5 J=1,8
            A = A + DBLE(KS2(J,I))*FA(J)
 5       CONTINUE
         S2 = S2 + ( SS2(1,I)*SIN(A) + SS2(2,I)*COS(A) )
 6    CONTINUE

      DO 8 I = NS3,1,-1
         A = 0D0
         DO 7 J=1,8
            A = A + DBLE(KS3(J,I))*FA(J)
 7       CONTINUE
         S3 = S3 + ( SS3(1,I)*SIN(A) + SS3(2,I)*COS(A) )
 8    CONTINUE

      DO 10 I = NS4,1,-1
         A = 0D0
         DO 9 J=1,8
            A = A + DBLE(KS4(J,I))*FA(J)
 9       CONTINUE
         S4 = S4 + ( SS4(1,I)*SIN(A) + SS4(2,I)*COS(A) )
 10   CONTINUE

      iau_S00 = ( S0 +
     :          ( S1 +
     :          ( S2 +
     :          ( S3 +
     :          ( S4 +
     :            S5 * T ) * T ) * T ) * T ) * T ) * DAS2R - X*Y/2D0

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAD03 ( T )
*+
*  - - - - - - - - - -
*   i a u _ F A D 0 3
*  - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean elongation of the Moon from the Sun.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAD03   d    D, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     is from Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Arcseconds to radians.
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle.
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean elongation of the Moon from the Sun (IERS Conventions 2003).
      iau_FAD03 = MOD (      1072260.703692D0 +
     :                  T*( 1602961601.2090D0 +
     :                  T*(        - 6.3706D0 +
     :                  T*(          0.006593D0 +
     :                  T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAE03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A E 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of Earth.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAE03   d    mean longitude of Earth, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     comes from Souchay et al. (1999) after Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  2Pi.
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of Earth (IERS Conventions 2003).
      iau_FAE03= MOD ( 1.753470314D0 + 628.3075849991D0 * T, D2PI )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAF03 ( T )
*+
*  - - - - - - - - - -
*   i a u _ F A F 0 3
*  - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of the Moon minus mean longitude of the ascending
*  node.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAF03   d    F, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     is from Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Arcseconds to radians.
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle.
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of the Moon minus that of the ascending node
*  (IERS Conventions 2003).
      iau_FAF03 = MOD (       335779.526232D0 +
     :                  T*( 1739527262.8478D0 +
     :                  T*(       - 12.7512D0 +
     :                  T*(       -  0.001037D0 +
     :                  T*(          0.00000417D0 )))), TURNAS ) * DAS2R

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAJU03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A J U 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of Jupiter.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAJU03  d    mean longitude of Jupiter, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     comes from Souchay et al. (1999) after Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  2Pi.
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of Jupiter (IERS Conventions 2003).
      iau_FAJU03= MOD ( 0.599546497D0 + 52.9690962641D0 * T, D2PI )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAL03 ( T )
*+
*  - - - - - - - - - -
*   i a u _ F A L 0 3
*  - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean anomaly of the Moon.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAL03   d    l, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     is from Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Arcseconds to radians.
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle.
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean anomaly of the Moon (IERS Conventions 2003).
      iau_FAL03 = MOD (       485868.249036D0 +
     :                  T*( 1717915923.2178D0 +
     :                  T*(         31.8792D0 +
     :                  T*(          0.051635D0 +
     :                  T*(        - 0.00024470D0 )))), TURNAS ) * DAS2R

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FALP03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A L P 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean anomaly of the Sun.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FALP03  d    l', radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     is from Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Arcseconds to radians.
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle.
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean anomaly of the Sun (IERS Conventions 2003).
      iau_FALP03 = MOD (     1287104.793048D0 +
     :                   T*( 129596581.0481D0 +
     :                   T*(       - 0.5532D0 +
     :                   T*(         0.000136D0 +
     :                   T*(       - 0.00001149D0 )))), TURNAS ) * DAS2R

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAMA03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A M A 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of Mars.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAMA03  d    mean longitude of Mars, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     comes from Souchay et al. (1999) after Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  2Pi.
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of Mars (IERS Conventions 2003).
      iau_FAMA03= MOD ( 6.203480913D0 + 334.0612426700D0 * T, D2PI )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAME03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A M E 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of Mercury.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAME03  d    mean longitude of Mercury, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     comes from Souchay et al. (1999) after Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  2Pi.
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of Mercury (IERS Conventions 2003).
      iau_FAME03 = MOD ( 4.402608842D0 + 2608.7903141574D0 * T, D2PI )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAOM03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A O M 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of the Moon's ascending node.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAOM03  d    Omega, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     is from Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Arcseconds to radians.
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle.
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of the Moon's ascending node (IERS Conventions 2003).
      iau_FAOM03 = MOD (      450160.398036D0 +
     :                   T*( - 6962890.5431D0 +
     :                   T*(         7.4722D0 +
     :                   T*(         0.007702D0 +
     :                   T*(       - 0.00005939D0 )))), TURNAS ) * DAS2R

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAPA03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A P A 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  general accumulated precession in longitude.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAPA03  d    general precession in longitude, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003).  It
*     is taken from Kinoshita & Souchay (1990) and comes originally from
*     Lieske et al. (1977).
*
*  References:
*
*     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
*     48, 187
*
*     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
*     Astron.Astrophys. 58, 1-16
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  General accumulated precession in longitude.
      iau_FAPA03= ( 0.024381750D0 + 0.00000538691D0 * T ) * T

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FASA03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A S A 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of Saturn.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FASA03  d    mean longitude of Saturn, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     comes from Souchay et al. (1999) after Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  2Pi.
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of Saturn (IERS Conventions 2003).
      iau_FASA03= MOD ( 0.874016757D0 + 21.3299104960D0 * T, D2PI )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAUR03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A U R 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of Uranus.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAUR03  d    mean longitude of Uranus, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     is adapted from Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  2Pi.
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of Uranus (IERS Conventions 2003).
      iau_FAUR03= MOD ( 5.481293872D0 + 7.4781598567D0 * T, D2PI )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_FAVE03 ( T )
*+
*  - - - - - - - - - - -
*   i a u _ F A V E 0 3
*  - - - - - - - - - - -
*
*  Fundamental argument, IERS Conventions (2003):
*  mean longitude of Venus.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     T           d    TDB, Julian centuries since J2000.0 (Note 1)
*
*  Returned:
*     iau_FAVE03  d    mean longitude of Venus, radians (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use TT,
*     which makes no significant difference.
*
*  2) The expression used is as adopted in IERS Conventions (2003) and
*     comes from Souchay et al. (1999) after Simon et al. (1994).
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  2Pi.
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Mean longitude of Venus (IERS Conventions 2003).
      iau_FAVE03= MOD ( 3.176146697D0 + 1021.3285546211D0 * T, D2PI )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      SUBROUTINE iau_PN00 ( DATE1, DATE2, DPSI, DEPS,
     :                      EPSA, RB, RP, RBP, RN, RBPN )
*+
*  - - - - - - - - -
*   i a u _ P N 0 0
*  - - - - - - - - -
*
*  Precession-nutation, IAU 2000 model:  a multi-purpose routine,
*  supporting classical (equinox-based) use directly and CIO-based
*  use indirectly.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1,DATE2   d       TT as a 2-part Julian Date (Note 1)
*     DPSI,DEPS     d       nutation (Note 2)
*
*  Returned:
*     EPSA          d       mean obliquity (Note 3)
*     RB          d(3,3)    frame bias matrix (Note 4)
*     RP          d(3,3)    precession matrix (Note 5)
*     RBP         d(3,3)    bias-precession matrix (Note 6)
*     RN          d(3,3)    nutation matrix (Note 7)
*     RBPN        d(3,3)    GCRS-to-true matrix (Note 8)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The caller is responsible for providing the nutation components;
*     they are in longitude and obliquity, in radians and are with
*     respect to the equinox and ecliptic of date.  For high-accuracy
*     applications, free core nutation should be included as well as
*     any other relevant corrections to the position of the CIP.
*
*  3) The returned mean obliquity is consistent with the IAU 2000
*     precession-nutation models.
*
*  4) The matrix RB transforms vectors from GCRS to J2000.0 mean equator
*     and equinox by applying frame bias.
*
*  5) The matrix RP transforms vectors from J2000.0 mean equator and
*     equinox to mean equator and equinox of date by applying
*     precession.
*
*  6) The matrix RBP transforms vectors from GCRS to mean equator and
*     equinox of date by applying frame bias then precession.  It is the
*     product RP x RB.
*
*  7) The matrix RN transforms vectors from mean equator and equinox of
*     date to true equator and equinox of date by applying the nutation
*     (luni-solar + planetary).
*
*  8) The matrix RBPN transforms vectors from GCRS to true equator and
*     equinox of date.  It is the product RN x RBP, applying frame bias,
*     precession and nutation in that order.
*
*  Called:
*     iau_PR00     IAU 2000 precession adjustments
*     iau_OBL80    mean obliquity, IAU 1980
*     iau_BP00     frame bias and precession matrices, IAU 2000
*     iau_NUMAT    form nutation matrix
*     iau_RXR      product of two r-matrices
*
*  Reference:
*
*     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
*     "Expressions for the Celestial Intermediate Pole and Celestial
*     Ephemeris Origin consistent with the IAU 2000A precession-nutation
*     model", Astron.Astrophys. 400, 1145-1154 (2003)
*
*     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
*          intermediate origin" (CIO) by IAU 2006 Resolution 2.
*
*  This revision:  2010 January 18
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS,
     :                 EPSA, RB(3,3), RP(3,3), RBP(3,3),
     :                 RN(3,3), RBPN(3,3)

      DOUBLE PRECISION DPSIPR, DEPSPR

      DOUBLE PRECISION iau_OBL80

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  IAU 2000 precession-rate adjustments.
      CALL iau_PR00 ( DATE1, DATE2, DPSIPR, DEPSPR )

*  Mean obliquity, consistent with IAU 2000 precession-nutation.
      EPSA = iau_OBL80 ( DATE1, DATE2 ) + DEPSPR

*  Frame bias and precession matrices and their product.
      CALL iau_BP00 ( DATE1, DATE2, RB, RP, RBP )

*  Nutation matrix.
      CALL iau_NUMAT ( EPSA, DPSI, DEPS, RN )

*  Bias-precession-nutation matrix (classical).
      CALL iau_RXR ( RN, RBP, RBPN )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      DOUBLE PRECISION FUNCTION iau_OBL80 ( DATE1, DATE2 )
*+
*  - - - - - - - - - -
*   i a u _ O B L 8 0
*  - - - - - - - - - -
*
*  Mean obliquity of the ecliptic, IAU 1980 model.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2     d      TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     iau_OBL80       d      obliquity of the ecliptic (radians, Note 2)
*
*  Notes:
*
*  1) The date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TDB)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*             DATE1         DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The result is the angle between the ecliptic and mean equator of
*     date DATE1+DATE2.
*
*  Reference:
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P. Kenneth Seidelmann (ed), University Science Books (1992),
*     Expression 3.222-1 (p114).
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  Mean obliquity of date.
      iau_OBL80 = DAS2R * ( 84381.448D0 +
     :                      ( -46.8150D0 +
     :                       ( -0.00059D0 +
     :                          0.001813D0 * T ) * T ) * T )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      SUBROUTINE iau_RXR ( A, B, ATB )
*+
*  - - - - - - - -
*   i a u _ R X R
*  - - - - - - - -
*
*  Multiply two r-matrices.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     A        d(3,3)    first r-matrix
*     B        d(3,3)    second r-matrix
*
*  Returned:
*     ATB      d(3,3)    A * B
*
*  Called:
*     iau_CR       copy r-matrix
*
*  This revision:  2000 November 25
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION A(3,3), B(3,3), ATB(3,3)

      INTEGER I, J, K
      DOUBLE PRECISION W, WM(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 3 I=1,3
         DO 2 J=1,3
            W = 0D0
            DO 1 K=1,3
               W = W + A(I,K)*B(K,J)
 1          CONTINUE
            WM(I,J) = W
 2       CONTINUE
 3    CONTINUE
      CALL iau_CR ( WM, ATB )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      SUBROUTINE iau_CR ( R, C )
*+
*  - - - - - - -
*   i a u _ C R
*  - - - - - - -
*
*  Copy an r-matrix.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     R        d(3,3)    r-matrix to be copied
*
*  Returned:
*     C        d(3,3)    copy
*
*  Called:
*     iau_CP       copy p-vector
*
*  This revision:  2000 November 25
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION R(3,3), C(3,3)

      INTEGER I

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 1 I=1,3
         CALL iau_CP ( R(1,I), C(1,I) )
 1    CONTINUE

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_CP ( P, C )
*+
*  - - - - - - -
*   i a u _ C P
*  - - - - - - -
*
*  Copy a p-vector.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     P        d(3)     p-vector to be copied
*
*  Returned:
*     C        d(3)     copy
*
*  This revision:  2000 November 25
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION P(3), C(3)

      INTEGER I

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 1 I=1,3
         C(I) = P(I)
 1    CONTINUE

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_NUMAT ( EPSA, DPSI, DEPS, RMATN )
*+
*  - - - - - - - - - -
*   i a u _ N U M A T
*  - - - - - - - - - -
*
*  Form the matrix of nutation.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     EPSA          d       mean obliquity of date (Note 1)
*     DPSI,DEPS     d       nutation (Note 2)
*
*  Returned:
*     RMATN       d(3,3)    nutation matrix (Note 3)
*
*  Notes:
*
*
*  1) The supplied mean obliquity EPSA, must be consistent with the
*     precession-nutation models from which DPSI and DEPS were obtained.
*
*  2) The caller is responsible for providing the nutation components;
*     they are in longitude and obliquity, in radians and are with
*     respect to the equinox and ecliptic of date.
*
*  3) The matrix operates in the sense V(true) = RMATN * V(mean),
*     where the p-vector V(true) is with respect to the true
*     equatorial triad of date and the p-vector V(mean) is with
*     respect to the mean equatorial triad of date.
*
*  Called:
*     iau_IR       initialize r-matrix to identity
*     iau_RX       rotate around X-axis
*     iau_RZ       rotate around Z-axis
*
*  Reference:
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P. Kenneth Seidelmann (ed), University Science Books (1992),
*     Section 3.222-3 (p114).
*
*  This revision:  2006 November 13
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION EPSA, DPSI, DEPS, RMATN(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Build the rotation matrix.
      CALL iau_IR ( RMATN )
      CALL iau_RX ( EPSA, RMATN )
      CALL iau_RZ ( -DPSI, RMATN )
      CALL iau_RX ( -(EPSA+DEPS), RMATN )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_IR ( R )
*+
*  - - - - - - -
*   i a u _ I R
*  - - - - - - -
*
*  Initialize an r-matrix to the identity matrix.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Returned:
*     R        d(3,3)    r-matrix
*
*  Called:
*     iau_ZR       zero r-matrix
*
*  This revision:  2000 November 25
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION R(3,3)

      INTEGER I

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      CALL iau_ZR ( R )
      DO 1 I=1,3
         R(I,I) = 1D0
 1    CONTINUE

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_ZR ( R )
*+
*  - - - - - - -
*   i a u _ Z R
*  - - - - - - -
*
*  Initialize an r-matrix to the null matrix.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Returned:
*     R        d(3,3)    r-matrix
*
*  This revision:  2000 November 25
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION R(3,3)

      INTEGER I, J

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 2 J=1,3
         DO 1 I=1,3
            R(I,J) = 0D0
 1       CONTINUE
 2    CONTINUE

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_RX ( PHI, R )
*+
*  - - - - - - -
*   i a u _ R X
*  - - - - - - -
*
*  Rotate an r-matrix about the x-axis.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     PHI      d         angle (radians)
*
*  Given and returned:
*     R        d(3,3)    r-matrix
*
*  Sign convention:  The matrix can be used to rotate the
*  reference frame of a vector.  Calling this routine with
*  positive PHI incorporates in the matrix an additional
*  rotation, about the x-axis, anticlockwise as seen looking
*  towards the origin from positive x.
*
*  Called:
*     iau_IR       initialize r-matrix to identity
*     iau_RXR      product of two r-matrices
*     iau_CR       copy r-matrix
*
*  This revision:  2006 November 13
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION PHI, R(3,3)

      DOUBLE PRECISION S, C, A(3,3), W(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Matrix representing new rotation.
      S = SIN(PHI)
      C = COS(PHI)
      CALL iau_IR ( A )
      A(2,2) = C
      A(3,2) = -S
      A(2,3) = S
      A(3,3) = C

*  Rotate.
      CALL iau_RXR ( A, R, W )

*  Return result.
      CALL iau_CR ( W, R )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_RY ( THETA, R )
*+
*  - - - - - - -
*   i a u _ R Y
*  - - - - - - -
*
*  Rotate an r-matrix about the y-axis.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     THETA    d         angle (radians)
*
*  Given and returned:
*     R        d(3,3)    r-matrix
*
*  Sign convention:  The matrix can be used to rotate the
*  reference frame of a vector.  Calling this routine with
*  positive THETA incorporates in the matrix an additional
*  rotation, about the y-axis, anticlockwise as seen looking
*  towards the origin from positive y.
*
*  Called:
*     iau_IR       initialize r-matrix to identity
*     iau_RXR      product of two r-matrices
*     iau_CR       copy r-matrix
*
*  This revision:  2006 November 13
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION THETA, R(3,3)

      DOUBLE PRECISION S, C, A(3,3), W(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Matrix representing new rotation.
      S = SIN(THETA)
      C = COS(THETA)
      CALL iau_IR ( A )
      A(1,1) = C
      A(3,1) = S
      A(1,3) = -S
      A(3,3) = C

*  Rotate.
      CALL iau_RXR ( A, R, W )

*  Return result.
      CALL iau_CR ( W, R )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_RZ ( PSI, R )
*+
*  - - - - - - -
*   i a u _ R Z
*  - - - - - - -
*
*  Rotate an r-matrix about the z-axis.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     PSI      d         angle (radians)
*
*  Given and returned:
*     R        d(3,3)    r-matrix, rotated
*
*  Sign convention:  The matrix can be used to rotate the
*  reference frame of a vector.  Calling this routine with
*  positive PSI incorporates in the matrix an additional
*  rotation, about the z-axis, anticlockwise as seen looking
*  towards the origin from positive z.
*
*  Called:
*     iau_IR       initialize r-matrix to identity
*     iau_RXR      product of two r-matrices
*     iau_CR       copy r-matrix
*
*  This revision:  2006 November 13
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION PSI, R(3,3)

      DOUBLE PRECISION S, C, A(3,3), W(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Matrix representing new rotation.
      S = SIN(PSI)
      C = COS(PSI)
      CALL iau_IR ( A )
      A(1,1) = C
      A(2,1) = -S
      A(1,2) = S
      A(2,2) = C

*  Rotate.
      CALL iau_RXR ( A, R, W )

*  Return result.
      CALL iau_CR ( W, R )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_BP00 ( DATE1, DATE2, RB, RP, RBP )
*+
*  - - - - - - - - -
*   i a u _ B P 0 0
*  - - - - - - - - -
*
*  Frame bias and precession, IAU 2000.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2    d       TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     RB           d(3,3)    frame bias matrix (Note 2)
*     RP           d(3,3)    precession matrix (Note 3)
*     RBP          d(3,3)    bias-precession matrix (Note 4)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The matrix RB transforms vectors from GCRS to mean J2000.0 by
*     applying frame bias.
*
*  3) The matrix RP transforms vectors from J2000.0 mean equator and
*     equinox to mean equator and equinox of date by applying
*     precession.
*
*  4) The matrix RBP transforms vectors from GCRS to mean equator and
*     equinox of date by applying frame bias then precession.  It is the
*     product RP x RB.
*
*  Called:
*     iau_BI00     frame bias components, IAU 2000
*     iau_PR00     IAU 2000 precession adjustments
*     iau_IR       initialize r-matrix to identity
*     iau_RX       rotate around X-axis
*     iau_RY       rotate around Y-axis
*     iau_RZ       rotate around Z-axis
*     iau_RXR      product of two r-matrices
*
*  Reference:
*
*     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
*     "Expressions for the Celestial Intermediate Pole and Celestial
*     Ephemeris Origin consistent with the IAU 2000A precession-nutation
*     model", Astron.Astrophys. 400, 1145-1154 (2003)
*
*     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
*          intermediate origin" (CIO) by IAU 2006 Resolution 2.
*
*  This revision:  2010 January 18
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, RB(3,3), RP(3,3), RBP(3,3)

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  J2000.0 obliquity (Lieske et al. 1977)
      DOUBLE PRECISION EPS0
      PARAMETER ( EPS0 = 84381.448D0 * DAS2R )

      DOUBLE PRECISION T, DPSIBI, DEPSBI, DRA0,
     :                 PSIA77, OMA77, CHIA, DPSIPR, DEPSPR, PSIA, OMA

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  Frame bias.
      CALL iau_BI00 ( DPSIBI, DEPSBI, DRA0 )

*  Precession angles (Lieske et al. 1977)
      PSIA77 =        ( 5038.7784D0 +
     :                (   -1.07259D0 +
     :                (   -0.001147D0 ) * T ) * T ) * T * DAS2R
      OMA77  = EPS0 + (
     :                (    0.05127D0 +
     :                (   -0.007726D0 ) * T ) * T ) * T * DAS2R
      CHIA   =        (   10.5526D0 +
     :                (   -2.38064D0 +
     :                (   -0.001125D0 ) * T ) * T ) * T * DAS2R

*  Apply IAU 2000 precession corrections.
      CALL iau_PR00 ( DATE1, DATE2, DPSIPR, DEPSPR )
      PSIA = PSIA77 + DPSIPR
      OMA  = OMA77  + DEPSPR

*  Frame bias matrix: GCRS to J2000.0.
      CALL iau_IR ( RB )
      CALL iau_RZ ( DRA0, RB )
      CALL iau_RY ( DPSIBI*SIN(EPS0), RB )
      CALL iau_RX ( -DEPSBI, RB )

*  Precession matrix: J2000.0 to mean of date.
      CALL iau_IR ( RP )
      CALL iau_RX ( EPS0, RP )
      CALL iau_RZ ( -PSIA, RP )
      CALL iau_RX ( -OMA, RP )
      CALL iau_RZ ( CHIA, RP )

*  Bias-precession matrix: GCRS to mean of date.
      CALL iau_RXR ( RP, RB, RBP )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_BI00 ( DPSIBI, DEPSBI, DRA )
*+
*  - - - - - - - - -
*   i a u _ B I 0 0
*  - - - - - - - - -
*
*  Frame bias components of IAU 2000 precession-nutation models (part of
*  MHB2000 with additions).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Returned:
*     DPSIBI,DEPSBI  d   longitude and obliquity corrections
*     DRA            d   the ICRS RA of the J2000.0 mean equinox
*
*  Notes:
*
*  1) The frame bias corrections in longitude and obliquity (radians)
*     are required in order to correct for the offset between the GCRS
*     pole and the J2000.0 mean pole.  They define, with respect to the
*     GCRS frame, a J2000.0 mean pole that is consistent with the rest
*     of the IAU 2000A precession-nutation model.
*
*  2) In addition to the displacement of the pole, the complete
*     description of the frame bias requires also an offset in right
*     ascension.  This is not part of the IAU 2000A model, and is from
*     Chapront et al. (2002).  It is returned in radians.
*
*  3) This is a supplemented implementation of one aspect of the IAU
*     2000A nutation model, formally adopted by the IAU General Assembly
*     in 2000, namely MHB2000 (Mathews et al. 2002).
*
*  References:
*
*     Chapront, J., Chapront-Touze, M. & Francou, G., Astron.Astrophys.,
*     387, 700, 2002.
*
*     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
*     and precession   New nutation series for nonrigid Earth and
*     insights into the Earth's interior", J.Geophys.Res., 107, B4,
*     2002.  The MHB2000 code itself was obtained on 9th September 2002
*     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DPSIBI, DEPSBI, DRA

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  The frame bias corrections in longitude and obliquity
      DOUBLE PRECISION DPBIAS, DEBIAS
      PARAMETER ( DPBIAS = -0.041775D0 * DAS2R,
     :            DEBIAS = -0.0068192D0 * DAS2R )

*  The ICRS RA of the J2000.0 equinox (Chapront et al., 2002)
      DOUBLE PRECISION DRA0
      PARAMETER ( DRA0 = -0.0146D0 * DAS2R )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Return the results (which are fixed).
      DPSIBI = DPBIAS
      DEPSBI = DEBIAS
      DRA = DRA0

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END

      SUBROUTINE iau_PR00 ( DATE1, DATE2, DPSIPR, DEPSPR )
*+
*  - - - - - - - - -
*   i a u _ P R 0 0
*  - - - - - - - - -
*
*  Precession-rate part of the IAU 2000 precession-nutation models
*  (part of MHB2000).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2    d   TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     DPSIPR,DEPSPR  d   precession corrections (Notes 2,3)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The precession adjustments are expressed as "nutation components",
*     corrections in longitude and obliquity with respect to the J2000.0
*     equinox and ecliptic.
*
*  3) Although the precession adjustments are stated to be with respect
*     to Lieske et al. (1977), the MHB2000 model does not specify which
*     set of Euler angles are to be used and how the adjustments are to
*     be applied.  The most literal and straightforward procedure is to
*     adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
*     to add DPSIPR to psi_A and DEPSPR to both omega_A and eps_A.
*
*  4) This is an implementation of one aspect of the IAU 2000A nutation
*     model, formally adopted by the IAU General Assembly in 2000,
*     namely MHB2000 (Mathews et al. 2002).
*
*  References:
*
*     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
*     for the precession quantities based upon the IAU (1976) System of
*     Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)
*
*     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
*     and precession   New nutation series for nonrigid Earth and
*     insights into the Earth's interior", J.Geophys.Res., 107, B4,
*     2002.  The MHB2000 code itself was obtained on 9th September 2002
*     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
*
*     Wallace, P.T., "Software for Implementing the IAU 2000
*     Resolutions", in IERS Workshop 5.1 (2002).
*
*  This revision:  2009 December 15
*
*  SOFA release 2010-12-01
*
*  Copyright (C) 2010 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSIPR, DEPSPR

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

      DOUBLE PRECISION T

*  ------------------------------------
*  Precession and obliquity corrections (radians per century)
*  ------------------------------------

      DOUBLE PRECISION PRECOR, OBLCOR
      PARAMETER ( PRECOR = -0.29965D0 * DAS2R,
     :            OBLCOR = -0.02524D0 * DAS2R )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  Precession rate contributions with respect to IAU 1976/80.
      DPSIPR = PRECOR * T
      DEPSPR = OBLCOR * T

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2010
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The name(s) of all routine(s) in your derived work shall not
*        include the prefix "iau".
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall
*     acknowledge that the SOFA software was used in obtaining those
*     results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  6. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*----------------------------------------------------------------------

      END


      subroutine admint2 (AMP,F,P,NOUT, AMPIN2,IDTIN,PHIN2,NIN,ITIN)
*f2py intent(in) AMPIN2,IDTIN,PHIN2,NIN,ITIN
*f2py intent(out) AMP,F,P,NOUT
*+
*  - - - - - - - - - - -
*   A D M I N T
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine returns the ocean loading displacement amplitude,
*  frequency, and phase of a set of tidal constituents generated by
*  the Bos-Scherneck website at http://www.oso.chalmers.se/~loading/.
*  The variable nin is input as the number wanted, and the variable 
*  nout is returned as the number provided.  The constituents used
*  are stored in the arrays idd (Doodson number) and tamp
*  (Cartwright-Edden amplitude).  The actual amp and phase of each
*  of these are determined by spline interpolation of the real and
*  imaginary part of the admittance, as specified at a subset of the
*  constituents.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 1 model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as
*     a Class 1, 2, or 3 model.
*
*  Given:
*     AMPIN       d      Cartwright-Edden amplitude of tidal constituents
*     IDTIN       i      Doodson number of tidal constituents
*     PHIN        d      Phase of tidal constituents
*     NIN         i      Number of harmonics used
*
*  Returned:
*     AMP         d      Amplitude due to ocean loading
*     F           d      Frequency due to ocean loading
*     P           d      Phase due to ocean loading
*     NOUT        i      Number of harmonics returned
*
*  Notes:
*
*  1) The phase is determined for a time set in COMMON block /date/ in
*     the subroutine TDFRPH.
*  
*  2) The arrays F and P must be specified as double precision. 
*
*  Called:
*     TDFRPH             Returns frequency and phase of a tidal
*                        constituent with given Doodson number            
*     SPLINE             Sets up array for cubic spline interpolation
*     EVAL               Performs cubic spline interpolation 
*     SHELLS             Sorts an array using Shell Sort
*
*  Test case:
*     Test cases are provided in the main program HARDISP.F.
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*  Revisions:  
*  2009 June 17 B.E. Stetzler  Initial changes to header
*  2009 June 18 B.E. Stetzler  Used IMPLICIT NONE, declared more variables,
*                              and added D0 to DOUBLE PRECISION variables 
*  2009 August 19 B.E.Stetzler Capitalized all variables for FORTRAN 77
*                              compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE
      integer I,J,K,IT,KEY,LL,NCON,NDI,NIN,NLP,NSD,NT,II,KK,ITIN
      integer*8 IDD,IDTIN,NOUT

*+----------------------------------------------------------------------
*  The parameters below set the number of harmonics used in the prediction
*  (nt; This must also be set in the main program) and the number of
*  constituents whose amp and phase may be specified (ncon)
*-----------------------------------------------------------------------
      PARAMETER (NT=342)
      PARAMETER (NCON=20)

      REAL AIM,DI,DR,DTR,RF,RL,SCR,SDI,SDR,TAMP,ZDI,ZDR,
     .     EVAL,AM,RE,SF, AMPIN,PHIN
      DOUBLE PRECISION AMPIN2,PHIN2
      DOUBLE PRECISION AMP,F,FR,P,PR

C      DIMENSION AMPIN(*),IDTIN(6,*),PHIN(*)
C      DIMENSION AMP(*),F(*),P(*)
      DIMENSION AMPIN(11),IDTIN(6,11),PHIN(11),AMPIN2(11),PHIN2(11)
      DIMENSION AMP(342),F(342),P(342)
      DIMENSION ITIN(5)

*  Arrays containing information about all stored constituents
      DIMENSION IDD(6,NT),TAMP(NT)
      
      COMMON/DATE/IT(5)

*  Arrays containing information about the subset whose amp and phase may
*  be specified, and scratch arrays for the spline routines for which
*  at most ncon constituents may be specified.

      DIMENSION RL(NCON),AIM(NCON),RF(NCON),KEY(NCON),SCR(NCON),
     . ZDI(NCON),ZDR(NCON),DI(NCON),DR(NCON),SDI(NCON),SDR(NCON)
      DATA DTR/.01745329252/
      DATA RL/NCON*0.0/,AIM/NCON*0.0/,RF/NCON*0.0/
      DATA ZDI/NCON*0.0/,ZDR/NCON*0.0/,DI/NCON*0.0/,DR/NCON*0.0/
      DATA SDI/NCON*0.0/,SDR/NCON*0.0/
      DATA TAMP/
     .  .632208, .294107, .121046, .079915, .023818,-.023589, .022994,
     .  .019333,-.017871, .017192, .016018, .004671,-.004662,-.004519,
     .  .004470, .004467, .002589,-.002455,-.002172, .001972, .001947,
     .  .001914,-.001898, .001802, .001304, .001170, .001130, .001061,
     . -.001022,-.001017, .001014, .000901,-.000857, .000855, .000855,
     .  .000772, .000741, .000741,-.000721, .000698, .000658, .000654,
     . -.000653, .000633, .000626,-.000598, .000590, .000544, .000479,
     . -.000464, .000413,-.000390, .000373, .000366, .000366,-.000360,
     . -.000355, .000354, .000329, .000328, .000319, .000302, .000279,
     . -.000274,-.000272, .000248,-.000225, .000224,-.000223,-.000216,
     .  .000211, .000209, .000194, .000185,-.000174,-.000171, .000159,
     .  .000131, .000127, .000120, .000118, .000117, .000108, .000107,
     .  .000105,-.000102, .000102, .000099,-.000096, .000095,-.000089,
     . -.000085,-.000084,-.000081,-.000077,-.000072,-.000067, .000066,
     .  .000064, .000063, .000063, .000063, .000062, .000062,-.000060,
     .  .000056, .000053, .000051, .000050, .368645,-.262232,-.121995,
     . -.050208, .050031,-.049470, .020620, .020613, .011279,-.009530,
     . -.009469,-.008012, .007414,-.007300, .007227,-.007131,-.006644,
     .  .005249, .004137, .004087, .003944, .003943, .003420, .003418,
     .  .002885, .002884, .002160,-.001936, .001934,-.001798, .001690,
     .  .001689, .001516, .001514,-.001511, .001383, .001372, .001371,
     . -.001253,-.001075, .001020, .000901, .000865,-.000794, .000788,
     .  .000782,-.000747,-.000745, .000670,-.000603,-.000597, .000542,
     .  .000542,-.000541,-.000469,-.000440, .000438, .000422, .000410,
     . -.000374,-.000365, .000345, .000335,-.000321,-.000319, .000307,
     .  .000291, .000290,-.000289, .000286, .000275, .000271, .000263,
     . -.000245, .000225, .000225, .000221,-.000202,-.000200,-.000199,
     .  .000192, .000183, .000183, .000183,-.000170, .000169, .000168,
     .  .000162, .000149,-.000147,-.000141, .000138, .000136, .000136,
     .  .000127, .000127,-.000126,-.000121,-.000121, .000117,-.000116,
     . -.000114,-.000114,-.000114, .000114, .000113, .000109, .000108,
     .  .000106,-.000106,-.000106, .000105, .000104,-.000103,-.000100,
     . -.000100,-.000100, .000099,-.000098, .000093, .000093, .000090,
     . -.000088, .000083,-.000083,-.000082,-.000081,-.000079,-.000077,
     . -.000075,-.000075,-.000075, .000071, .000071,-.000071, .000068,
     .  .000068, .000065, .000065, .000064, .000064, .000064,-.000064,
     . -.000060, .000056, .000056, .000053, .000053, .000053,-.000053,
     .  .000053, .000053, .000052, .000050,-.066607,-.035184,-.030988,
     .  .027929,-.027616,-.012753,-.006728,-.005837,-.005286,-.004921,
     . -.002884,-.002583,-.002422, .002310, .002283,-.002037, .001883,
     . -.001811,-.001687,-.001004,-.000925,-.000844, .000766, .000766,
     . -.000700,-.000495,-.000492, .000491, .000483, .000437,-.000416,
     . -.000384, .000374,-.000312,-.000288,-.000273, .000259, .000245,
     . -.000232, .000229,-.000216, .000206,-.000204,-.000202, .000200,
     .  .000195,-.000190, .000187, .000180,-.000179, .000170, .000153,
     . -.000137,-.000119,-.000119,-.000112,-.000110,-.000110, .000107,
     . -.000095,-.000095,-.000091,-.000090,-.000081,-.000079,-.000079,
     .  .000077,-.000073, .000069,-.000067,-.000066, .000065, .000064,
     . -.000062, .000060, .000059,-.000056, .000055,-.000051/
      DATA IDD/
     .  2, 0, 0, 0, 0, 0,   2, 2,-2, 0, 0, 0,   2,-1, 0, 1, 0, 0,  
     .  2, 2, 0, 0, 0, 0,   2, 2, 0, 0, 1, 0,   2, 0, 0, 0,-1, 0,  
     .  2,-1, 2,-1, 0, 0,   2,-2, 2, 0, 0, 0,   2, 1, 0,-1, 0, 0,  
     .  2, 2,-3, 0, 0, 1,   2,-2, 0, 2, 0, 0,   2,-3, 2, 1, 0, 0,  
     .  2, 1,-2, 1, 0, 0,   2,-1, 0, 1,-1, 0,   2, 3, 0,-1, 0, 0,  
     .  2, 1, 0, 1, 0, 0,   2, 2, 0, 0, 2, 0,   2, 2,-1, 0, 0,-1,  
     .  2, 0,-1, 0, 0, 1,   2, 1, 0, 1, 1, 0,   2, 3, 0,-1, 1, 0,  
     .  2, 0, 1, 0, 0,-1,   2, 0,-2, 2, 0, 0,   2,-3, 0, 3, 0, 0,  
     .  2,-2, 3, 0, 0,-1,   2, 4, 0, 0, 0, 0,   2,-1, 1, 1, 0,-1,  
     .  2,-1, 3,-1, 0,-1,   2, 2, 0, 0,-1, 0,   2,-1,-1, 1, 0, 1,  
     .  2, 4, 0, 0, 1, 0,   2,-3, 4,-1, 0, 0,   2,-1, 2,-1,-1, 0,  
     .  2, 3,-2, 1, 0, 0,   2, 1, 2,-1, 0, 0,   2,-4, 2, 2, 0, 0,  
     .  2, 4,-2, 0, 0, 0,   2, 0, 2, 0, 0, 0,   2,-2, 2, 0,-1, 0,  
     .  2, 2,-4, 0, 0, 2,   2, 2,-2, 0,-1, 0,   2, 1, 0,-1,-1, 0,  
     .  2,-1, 1, 0, 0, 0,   2, 2,-1, 0, 0, 1,   2, 2, 1, 0, 0,-1,  
     .  2,-2, 0, 2,-1, 0,   2,-2, 4,-2, 0, 0,   2, 2, 2, 0, 0, 0,  
     .  2,-4, 4, 0, 0, 0,   2,-1, 0,-1,-2, 0,   2, 1, 2,-1, 1, 0,  
     .  2,-1,-2, 3, 0, 0,   2, 3,-2, 1, 1, 0,   2, 4, 0,-2, 0, 0,  
     .  2, 0, 0, 2, 0, 0,   2, 0, 2,-2, 0, 0,   2, 0, 2, 0, 1, 0,  
     .  2,-3, 3, 1, 0,-1,   2, 0, 0, 0,-2, 0,   2, 4, 0, 0, 2, 0,  
     .  2, 4,-2, 0, 1, 0,   2, 0, 0, 0, 0, 2,   2, 1, 0, 1, 2, 0,  
     .  2, 0,-2, 0,-2, 0,   2,-2, 1, 0, 0, 1,   2,-2, 1, 2, 0,-1,  
     .  2,-1, 1,-1, 0, 1,   2, 5, 0,-1, 0, 0,   2, 1,-3, 1, 0, 1,  
     .  2,-2,-1, 2, 0, 1,   2, 3, 0,-1, 2, 0,   2, 1,-2, 1,-1, 0,  
     .  2, 5, 0,-1, 1, 0,   2,-4, 0, 4, 0, 0,   2,-3, 2, 1,-1, 0,  
     .  2,-2, 1, 1, 0, 0,   2, 4, 0,-2, 1, 0,   2, 0, 0, 2, 1, 0,  
     .  2,-5, 4, 1, 0, 0,   2, 0, 2, 0, 2, 0,   2,-1, 2, 1, 0, 0,  
     .  2, 5,-2,-1, 0, 0,   2, 1,-1, 0, 0, 0,   2, 2,-2, 0, 0, 2,  
     .  2,-5, 2, 3, 0, 0,   2,-1,-2, 1,-2, 0,   2,-3, 5,-1, 0,-1,  
     .  2,-1, 0, 0, 0, 1,   2,-2, 0, 0,-2, 0,   2, 0,-1, 1, 0, 0,  
     .  2,-3, 1, 1, 0, 1,   2, 3, 0,-1,-1, 0,   2, 1, 0, 1,-1, 0,  
     .  2,-1, 2, 1, 1, 0,   2, 0,-3, 2, 0, 1,   2, 1,-1,-1, 0, 1,  
     .  2,-3, 0, 3,-1, 0,   2, 0,-2, 2,-1, 0,   2,-4, 3, 2, 0,-1,  
     .  2,-1, 0, 1,-2, 0,   2, 5, 0,-1, 2, 0,   2,-4, 5, 0, 0,-1,  
     .  2,-2, 4, 0, 0,-2,   2,-1, 0, 1, 0, 2,   2,-2,-2, 4, 0, 0,  
     .  2, 3,-2,-1,-1, 0,   2,-2, 5,-2, 0,-1,   2, 0,-1, 0,-1, 1,  
     .  2, 5,-2,-1, 1, 0,   1, 1, 0, 0, 0, 0,   1,-1, 0, 0, 0, 0,  
     .  1, 1,-2, 0, 0, 0,   1,-2, 0, 1, 0, 0,   1, 1, 0, 0, 1, 0,  
     .  1,-1, 0, 0,-1, 0,   1, 2, 0,-1, 0, 0,   1, 0, 0, 1, 0, 0,  
     .  1, 3, 0, 0, 0, 0,   1,-2, 2,-1, 0, 0,   1,-2, 0, 1,-1, 0,  
     .  1,-3, 2, 0, 0, 0,   1, 0, 0,-1, 0, 0,   1, 1, 0, 0,-1, 0,  
     .  1, 3, 0, 0, 1, 0,   1, 1,-3, 0, 0, 1,   1,-3, 0, 2, 0, 0,  
     .  1, 1, 2, 0, 0, 0,   1, 0, 0, 1, 1, 0,   1, 2, 0,-1, 1, 0,  
     .  1, 0, 2,-1, 0, 0,   1, 2,-2, 1, 0, 0,   1, 3,-2, 0, 0, 0,  
     .  1,-1, 2, 0, 0, 0,   1, 1, 1, 0, 0,-1,   1, 1,-1, 0, 0, 1,  
     .  1, 4, 0,-1, 0, 0,   1,-4, 2, 1, 0, 0,   1, 0,-2, 1, 0, 0,  
     .  1,-2, 2,-1,-1, 0,   1, 3, 0,-2, 0, 0,   1,-1, 0, 2, 0, 0,  
     .  1,-1, 0, 0,-2, 0,   1, 3, 0, 0, 2, 0,   1,-3, 2, 0,-1, 0,  
     .  1, 4, 0,-1, 1, 0,   1, 0, 0,-1,-1, 0,   1, 1,-2, 0,-1, 0,  
     .  1,-3, 0, 2,-1, 0,   1, 1, 0, 0, 2, 0,   1, 1,-1, 0, 0,-1,  
     .  1,-1,-1, 0, 0, 1,   1, 0, 2,-1, 1, 0,   1,-1, 1, 0, 0,-1,  
     .  1,-1,-2, 2, 0, 0,   1, 2,-2, 1, 1, 0,   1,-4, 0, 3, 0, 0,  
     .  1,-1, 2, 0, 1, 0,   1, 3,-2, 0, 1, 0,   1, 2, 0,-1,-1, 0,  
     .  1, 0, 0, 1,-1, 0,   1,-2, 2, 1, 0, 0,   1, 4,-2,-1, 0, 0,  
     .  1,-3, 3, 0, 0,-1,   1,-2, 1, 1, 0,-1,   1,-2, 3,-1, 0,-1,  
     .  1, 0,-2, 1,-1, 0,   1,-2,-1, 1, 0, 1,   1, 4,-2, 1, 0, 0,  
     .  1,-4, 4,-1, 0, 0,   1,-4, 2, 1,-1, 0,   1, 5,-2, 0, 0, 0,  
     .  1, 3, 0,-2, 1, 0,   1,-5, 2, 2, 0, 0,   1, 2, 0, 1, 0, 0,  
     .  1, 1, 3, 0, 0,-1,   1,-2, 0, 1,-2, 0,   1, 4, 0,-1, 2, 0,  
     .  1, 1,-4, 0, 0, 2,   1, 5, 0,-2, 0, 0,   1,-1, 0, 2, 1, 0,  
     .  1,-2, 1, 0, 0, 0,   1, 4,-2, 1, 1, 0,   1,-3, 4,-2, 0, 0,  
     .  1,-1, 3, 0, 0,-1,   1, 3,-3, 0, 0, 1,   1, 5,-2, 0, 1, 0,  
     .  1, 1, 2, 0, 1, 0,   1, 2, 0, 1, 1, 0,   1,-5, 4, 0, 0, 0,  
     .  1,-2, 0,-1,-2, 0,   1, 5, 0,-2, 1, 0,   1, 1, 2,-2, 0, 0,  
     .  1, 1,-2, 2, 0, 0,   1,-2, 2, 1, 1, 0,   1, 0, 3,-1, 0,-1,  
     .  1, 2,-3, 1, 0, 1,   1,-2,-2, 3, 0, 0,   1,-1, 2,-2, 0, 0,  
     .  1,-4, 3, 1, 0,-1,   1,-4, 0, 3,-1, 0,   1,-1,-2, 2,-1, 0,  
     .  1,-2, 0, 3, 0, 0,   1, 4, 0,-3, 0, 0,   1, 0, 1, 1, 0,-1,  
     .  1, 2,-1,-1, 0, 1,   1, 2,-2, 1,-1, 0,   1, 0, 0,-1,-2, 0,  
     .  1, 2, 0, 1, 2, 0,   1, 2,-2,-1,-1, 0,   1, 0, 0, 1, 2, 0,  
     .  1, 0, 1, 0, 0, 0,   1, 2,-1, 0, 0, 0,   1, 0, 2,-1,-1, 0,  
     .  1,-1,-2, 0,-2, 0,   1,-3, 1, 0, 0, 1,   1, 3,-2, 0,-1, 0,  
     .  1,-1,-1, 0,-1, 1,   1, 4,-2,-1, 1, 0,   1, 2, 1,-1, 0,-1,  
     .  1, 0,-1, 1, 0, 1,   1,-2, 4,-1, 0, 0,   1, 4,-4, 1, 0, 0,  
     .  1,-3, 1, 2, 0,-1,   1,-3, 3, 0,-1,-1,   1, 1, 2, 0, 2, 0,  
     .  1, 1,-2, 0,-2, 0,   1, 3, 0, 0, 3, 0,   1,-1, 2, 0,-1, 0,  
     .  1,-2, 1,-1, 0, 1,   1, 0,-3, 1, 0, 1,   1,-3,-1, 2, 0, 1,  
     .  1, 2, 0,-1, 2, 0,   1, 6,-2,-1, 0, 0,   1, 2, 2,-1, 0, 0,  
     .  1,-1, 1, 0,-1,-1,   1,-2, 3,-1,-1,-1,   1,-1, 0, 0, 0, 2,  
     .  1,-5, 0, 4, 0, 0,   1, 1, 0, 0, 0,-2,   1,-2, 1, 1,-1,-1,  
     .  1, 1,-1, 0, 1, 1,   1, 1, 2, 0, 0,-2,   1,-3, 1, 1, 0, 0,  
     .  1,-4, 4,-1,-1, 0,   1, 1, 0,-2,-1, 0,   1,-2,-1, 1,-1, 1,  
     .  1,-3, 2, 2, 0, 0,   1, 5,-2,-2, 0, 0,   1, 3,-4, 2, 0, 0,  
     .  1, 1,-2, 0, 0, 2,   1,-1, 4,-2, 0, 0,   1, 2, 2,-1, 1, 0,  
     .  1,-5, 2, 2,-1, 0,   1, 1,-3, 0,-1, 1,   1, 1, 1, 0, 1,-1,  
     .  1, 6,-2,-1, 1, 0,   1,-2, 2,-1,-2, 0,   1, 4,-2, 1, 2, 0,  
     .  1,-6, 4, 1, 0, 0,   1, 5,-4, 0, 0, 0,   1,-3, 4, 0, 0, 0,  
     .  1, 1, 2,-2, 1, 0,   1,-2, 1, 0,-1, 0,   0, 2, 0, 0, 0, 0,  
     .  0, 1, 0,-1, 0, 0,   0, 0, 2, 0, 0, 0,   0, 0, 0, 0, 1, 0,  
     .  0, 2, 0, 0, 1, 0,   0, 3, 0,-1, 0, 0,   0, 1,-2, 1, 0, 0,  
     .  0, 2,-2, 0, 0, 0,   0, 3, 0,-1, 1, 0,   0, 0, 1, 0, 0,-1,  
     .  0, 2, 0,-2, 0, 0,   0, 2, 0, 0, 2, 0,   0, 3,-2, 1, 0, 0,  
     .  0, 1, 0,-1,-1, 0,   0, 1, 0,-1, 1, 0,   0, 4,-2, 0, 0, 0,  
     .  0, 1, 0, 1, 0, 0,   0, 0, 3, 0, 0,-1,   0, 4, 0,-2, 0, 0,  
     .  0, 3,-2, 1, 1, 0,   0, 3,-2,-1, 0, 0,   0, 4,-2, 0, 1, 0,  
     .  0, 0, 2, 0, 1, 0,   0, 1, 0, 1, 1, 0,   0, 4, 0,-2, 1, 0,  
     .  0, 3, 0,-1, 2, 0,   0, 5,-2,-1, 0, 0,   0, 1, 2,-1, 0, 0,  
     .  0, 1,-2, 1,-1, 0,   0, 1,-2, 1, 1, 0,   0, 2,-2, 0,-1, 0,  
     .  0, 2,-3, 0, 0, 1,   0, 2,-2, 0, 1, 0,   0, 0, 2,-2, 0, 0,  
     .  0, 1,-3, 1, 0, 1,   0, 0, 0, 0, 2, 0,   0, 0, 1, 0, 0, 1,  
     .  0, 1, 2,-1, 1, 0,   0, 3, 0,-3, 0, 0,   0, 2, 1, 0, 0,-1,  
     .  0, 1,-1,-1, 0, 1,   0, 1, 0, 1, 2, 0,   0, 5,-2,-1, 1, 0,  
     .  0, 2,-1, 0, 0, 1,   0, 2, 2,-2, 0, 0,   0, 1,-1, 0, 0, 0,  
     .  0, 5, 0,-3, 0, 0,   0, 2, 0,-2, 1, 0,   0, 1, 1,-1, 0,-1,  
     .  0, 3,-4, 1, 0, 0,   0, 0, 2, 0, 2, 0,   0, 2, 0,-2,-1, 0,  
     .  0, 4,-3, 0, 0, 1,   0, 3,-1,-1, 0, 1,   0, 0, 2, 0, 0,-2,  
     .  0, 3,-3, 1, 0, 1,   0, 2,-4, 2, 0, 0,   0, 4,-2,-2, 0, 0,  
     .  0, 3, 1,-1, 0,-1,   0, 5,-4, 1, 0, 0,   0, 3,-2,-1,-1, 0,  
     .  0, 3,-2, 1, 2, 0,   0, 4,-4, 0, 0, 0,   0, 6,-2,-2, 0, 0,  
     .  0, 5, 0,-3, 1, 0,   0, 4,-2, 0, 2, 0,   0, 2, 2,-2, 1, 0,  
     .  0, 0, 4, 0, 0,-2,   0, 3,-1, 0, 0, 0,   0, 3,-3,-1, 0, 1,  
     .  0, 4, 0,-2, 2, 0,   0, 1,-2,-1,-1, 0,   0, 2,-1, 0, 0,-1,  
     .  0, 4,-4, 2, 0, 0,   0, 2, 1, 0, 1,-1,   0, 3,-2,-1, 1, 0,  
     .  0, 4,-3, 0, 1, 1,   0, 2, 0, 0, 3, 0,   0, 6,-4, 0, 0, 0/

*  Initialize variables.
      K   = 0
      NLP = 0
      NDI = 0
      NSD = 0
      
      IT = int(ITIN)
            
C      print *, IDTIN(1:6,11)
      AMPIN = real(AMPIN2)
      PHIN = real(PHIN2)

      DO LL=1,NIN
*  See if Doodson numbers match
         DO KK=1,NT
            II = 0
            DO I=1,6
               II = II + abs(IDD(I,KK)-IDTIN(I,LL))
*      print *, 'I=', I, 'KK=', KK, 'LL=', LL
*      print *, IDD(I,KK),IDTIN(I,LL),abs(IDD(I,KK)-IDTIN(I,LL))
            ENDDO
            IF(II.EQ.0) GO TO 5
         ENDDO
*  If you have a match, put line into array
 5       IF(II.EQ.0.AND.K.LT.NCON) THEN
            K = K + 1
            RL(K) = AMPIN(LL)*COS(DTR*PHIN(LL))/ABS(TAMP(KK))
            AIM(K)= AMPIN(LL)*SIN(DTR*PHIN(LL))/ABS(TAMP(KK))
*+---------------------------------------------------------------------
*  Now have real and imaginary parts of admittance, scaled by Cartwright-
*  Edden amplitude. Admittance phase is whatever was used in the original
*  expression. (Usually phase is given relative to some reference,
*  but amplitude is in absolute units). Next get frequency.
*----------------------------------------------------------------------
            CALL TDFRPH(IDD(1,KK),FR,PR)
            RF(K) = FR
         ENDIF
      ENDDO

*      pause
*+---------------------------------------------------------------------
*  Done going through constituents; there are k of them.
*  Have specified admittance at a number of points. Sort these by frequency
*  and separate diurnal and semidiurnal, recopying admittances to get them
*  in order using Shell Sort.
*----------------------------------------------------------------------

      CALL SHELLS(RF,KEY,K)
      DO I=1,K
         IF(RF(I).LT.0.5) NLP = NLP + 1
         IF(RF(I).LT.1.5.AND.RF(I).GT.0.5) NDI = NDI + 1
         IF(RF(I).LT.2.5.AND.RF(I).GT.1.5) NSD = NSD + 1
         SCR(I) = RL(KEY(I))
      ENDDO
*      print *, SCR
*      pause
      DO I=1,K
         RL(I) = SCR(I)
         SCR(I) = AIM(KEY(I))
      ENDDO
      DO I=1,K
         AIM(I) = SCR(I)
      ENDDO
*+---------------------------------------------------------------------
*  now set up splines (8 cases - four species, each real and imaginary)
*  We have to allow for the case when there are no constituent amplitudes
*  for the long-period tides.
*----------------------------------------------------------------------
*      print *, RF
      IF(NLP.NE.0) CALL SPLINE(NLP,RF,RL,ZDR,SCR)
      IF(NLP.NE.0) CALL SPLINE(NLP,RF,AIM,ZDI,SCR)
      CALL SPLINE(NDI,RF(NLP+1),RL(NLP+1),DR,SCR)
      CALL SPLINE(NDI,RF(NLP+1),AIM(NLP+1),DI,SCR)
      CALL SPLINE(NSD,RF(NLP+NDI+1),RL(NLP+NDI+1),SDR,SCR)
      CALL SPLINE(NSD,RF(NLP+NDI+1),AIM(NLP+NDI+1),SDI,SCR)
*  Evaluate all harmonics using the interpolated admittance
      J = 1
      DO I=1,NT
         IF(IDD(1,I).EQ.0.AND.NLP.EQ.0) GO TO 11
         CALL TDFRPH(IDD(1,I),F(J),P(J))
*  Compute phase corrections to equilibrium tide using function EVAL
         IF(IDD(1,I).EQ.0) P(J) = P(J) + 180.
         IF(IDD(1,I).EQ.1) P(J) = P(J) + 90.
         SF = F(J)
         IF(IDD(1,I).EQ.0) RE = EVAL(SF,NLP,RF,RL,ZDR)
         IF(IDD(1,I).EQ.0) AM = EVAL(SF,NLP,RF,AIM,ZDI)
         IF(IDD(1,I).EQ.1) RE = EVAL(SF,NDI,RF(NLP+1),RL(NLP+1),DR)
         IF(IDD(1,I).EQ.1) AM = EVAL(SF,NDI,RF(NLP+1),AIM(NLP+1),DI)
         IF(IDD(1,I).EQ.2) RE =
     .      EVAL(SF,NSD,RF(NLP+NDI+1),RL(NLP+NDI+1),SDR)
         IF(IDD(1,I).EQ.2) AM =
     .      EVAL(SF,NSD,RF(NLP+NDI+1),AIM(NLP+NDI+1),SDI)
*      print *, RF(NLP+NDI+1:NLP+NDI+NSD), RL(NLP+NDI+1:NLP+NDI+NSD)
*      pause
*      print *, RE, AM
         AMP(J) = TAMP(I)*SQRT(RE**2+AM**2)
         P(J) = P(J) + ATAN2(AM,RE)/DTR
         IF(P(J).GT.180) P(J)=P(J)-360.
*         print *, P(J)
         J = J + 1
 11      CONTINUE
      ENDDO
      NOUT = J - 1
*      open(unit=23,file="logt.log",access="append",status="old")
*      write(23,*), NOUT
*      close(23)
      RETURN

*  Finished.


      END
      
      SUBROUTINE TDFRPH (IDOOD,FREQ,PHASE)
*+
*  - - - - - - - - - - -
*   T D F R P H 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine returns the frequency and phase of a tidal
*  constituent when its Doodson number is given as input. 
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 1 model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as
*     a Class 1, 2, or 3 model.
*
*  Given:
*     idood       i      Doodson number of a tidal constituent
*
*  Returned:
*     freq        d      Frequency of a tidal constituent
*     phase       d      Phase of a tidal constituent (Note 1)
*
*  Notes:
*
*  1) The phases must be decreased by 90 degrees if the sum of the order 
*     and the species number is odd (as for the 2nd degree diurnals, and 
*     3rd degree low frequency and semidiurnals).
*     
*     These phases may need further adjustment to allow for the spherical
*     harmonic normalization used; e.g. for that used for the potential
*     by Cartwright and Tayler, 180 degrees must be added for (species,
*     order) = (1,2), (1,3), or (3,3). 
*
*  Called:
*     TOYMD     Converts year-day of year to year-month-day format 
*     LEAP      Returns true if a given year is a leap year
*     JULDAT    Converts Gregorian date to Julian date 
*     ETUTC     Returns difference of Epheremis Time (ET) and 
*               Coordinated Universal Time (UTC) 
*
*  Test case:
*     given input: For June 25, 2009 0 Hr 0 Min, M2 tide
*                  DATA IDOOD = / 2, 0, 0, 0, 0, 0 /  
*
*     expected output: FREQ = 1.93227361605688D0
*                      PHASE = 132.8193176853237674D0
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:  
*  2009 June   15 B.E.Stetzler  Initial changes to code 
*  2009 August 19 B.E.Stetzler  Capitalized all variables for FORTRAN
*                               77 compatibility
*  2010 March  19 B.E.Stetzler  Provided test case
*-----------------------------------------------------------------------

      IMPLICIT NONE
C      SAVE ITMSAVE,D,DD
      INTEGER I,INITIAL,ITM,ITM2,ITMSAVE,JD,JULDAT,LEAP
      INTEGER*8 IDOOD
      DOUBLE PRECISION YEAR,DELTA,FREQ,PHASE,
     .                 D,DAYFR,DD,DJD,F1,F2,F3,F4,F5,
     .                 FD1,FD2,FD3,FD4,FD5,T
      DIMENSION IDOOD(6),ITM2(6),ITMSAVE(5),D(6),DD(6)

* Common block 'date' stores time information in Universal Time (UT)

      COMMON/DATE/ITM(5)
C      DATA ITMSAVE/5*0/
*------------------------------------------------------------------------
*  Test to see if time has changed; if so, set the phases and frequencies
*  for each of the Doodson arguments
*------------------------------------------------------------------------
      INITIAL=0

C      DO I=1,5
C        IF(ITM(I).NE.ITMSAVE(I)) INITIAL=1
C      ENDDO

C      IF(INITIAL.EQ.1) THEN
C        DO I=1,5
C           ITMSAVE(I) = ITM(I)
C        ENDDO

* Convert times to Julian days (UT) then to Julian centuries from J2000.0
*   (ET)
        CALL TOYMD(ITM,ITM2)
        JD = JULDAT(ITM2)
        DAYFR=  ITM(3)/24.0D0 + ITM(4)/1440.0D0 + ITM(5)/84600.0D0
        YEAR=ITM(1)+(ITM(2)+DAYFR)/(365.0D0+LEAP(ITM(1)))
        CALL ETUTC(YEAR,DELTA)
        DJD= JD - 0.5D0 + DAYFR
        T = (DJD - 2451545.0D0 + DELTA/86400.0D0)/36525.0D0

* IERS expressions for the Delaunay arguments, in degrees

        F1 =     134.9634025100D0 +
     .    T*( 477198.8675605000D0 +
     .    T*(      0.0088553333D0 +
     .    T*(      0.0000143431D0 +
     .    T*(     -0.0000000680D0 ))))
        F2 =     357.5291091806D0 +
     .    T*(  35999.0502911389D0 +
     .    T*(     -0.0001536667D0 +
     .    T*(      0.0000000378D0 +
     .    T*(     -0.0000000032D0 ))))
        F3 =      93.2720906200D0 +
     .    T*( 483202.0174577222D0 +
     .    T*(     -0.0035420000D0 +
     .    T*(     -0.0000002881D0 +
     .    T*(      0.0000000012D0 ))))
        F4 =     297.8501954694D0 +
     .    T*( 445267.1114469445D0 +
     .    T*(     -0.0017696111D0 +
     .    T*(      0.0000018314D0 +
     .    T*(     -0.0000000088D0 ))))
        F5 =     125.0445550100D0 +
     .    T*(  -1934.1362619722D0 +
     .    T*(      0.0020756111D0 +
     .    T*(      0.0000021394D0 +
     .    T*(     -0.0000000165D0 ))))

*  Convert to Doodson (Darwin) variables

        D(1) = 360.0D0*DAYFR - F4
        D(2) = F3 + F5
        D(3) = D(2) - F4
        D(4) = D(2) - F1
        D(5) = -F5
        D(6) = D(3) - F2

*  Find frequencies of Delauney variables (in cycles/day), and from these
*  the same for the Doodson arguments

        FD1 =  0.0362916471D0 + 0.0000000013D0*T
        FD2 =  0.0027377786D0
        FD3 =  0.0367481951D0 - 0.0000000005D0*T
        FD4 =  0.0338631920D0 - 0.0000000003D0*T
        FD5 = -0.0001470938D0 + 0.0000000003D0*T
        DD(1) = 1.0D0 - FD4
        DD(2) = FD3 + FD5
        DD(3) = DD(2) - FD4
        DD(4) = DD(2) - FD1
        DD(5) = -FD5
        DD(6) = DD(3) - FD2
C      ENDIF

*  End of intialization (likely to be called only once)

*  Compute phase and frequency of the given tidal constituent

      FREQ=0.0D0
      PHASE=0.0D0
      DO I = 1,6
         FREQ =  FREQ + IDOOD(I)*DD(I)
         PHASE = PHASE + IDOOD(I)*D(I)
      ENDDO

*      print *, FREQ
*      print *, PHASE
*      print *,'\n'
* Adjust phases so that they fall in the positive range 0 to 360
      PHASE = DMOD(PHASE,360.0D0)
      IF(PHASE.LT.0.0D0) PHASE = PHASE + 360.0D0

      RETURN

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END

      SUBROUTINE TOYMD (IT1,IT2)
*+
*  - - - - - - - - - - -
*   T O Y M D
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine converts times given in it1 expressed in year and
*  day of year to year-month-day in it2.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     it1           i(2)      time given in year and day of year (Note 1)
*
*  Returned:
*     it2           i(3)      time given in year-month-day format
*
*  Notes:
*
*  1) The time is split into a year, given as it1(1) and the day of the
*     year, given as it1(2). 
*
*  Called:
*    LEAP 
*
*  Test case:
*    Given input:  it1(1) = 2008
*                  it1(2) = 120
*
*    Expected output: it2(1) = 2008
*                     it2(2) = 4
*                     it2(3) = 29
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 April   20 B.E.Stetzler Initial standardization of subroutine 
*  2009 April   23 B.E.Stetzler Provided test case
*  2009 August  19 B.E.Stetzler Capitalized all variables for FORTRAN
*                               77 compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER IDN,IT1,IT2,JJ,M,MON,LEAP
*      DIMENSION IT1(*),IT2(*)
      DIMENSION IT1(2),IT2(3)

      IDN(M) = MOD((367*(M-2-12*((M-14)/12)))/12+29,365)
      MON(JJ,M) = (12*(JJ-29-M))/367 + 2 + (JJ-200)/169
      IT2(1) = IT1(1)
      IT2(2) = MON(IT1(2),LEAP(IT1(1)))
      IT2(3) = IT1(2) - IDN(IT2(2)) - LEAP(IT2(1))*((9+IT2(2))/12)
      RETURN

* Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END

      SUBROUTINE SPLINE (NN,X,U,S,A)
*+
*  - - - - - - - - -
*   S P L I N E
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to find an array s for the spline
*  interpolator function EVAL.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model	
* 
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given: This is a support routine of the main program HARDISP.F.
*     nn             i      number of data points supplied, which may be
*                           negative (Note 1)
*     x              d      array containing x-coordinates where function
*                           is sampled (Note 2)
*     u              d      array containing sample values that are to be
*                           interpolated 
*     a              d      working space array of dimension at least nn
*
*  Returned:
*     s              d      output array of 2nd derivative at sample points 
*
*  Notes:
*
*  1) If the user wishes to force the derivatives at the ends of the series
*     to assume specified values, he or she should put du(1)/dx and du(n)/dx
*     in the variables s1 and s2 and call the subroutine with nn = -(number
*     of terms in the series).  Normally a parabola is fitted through the 
*     1st and last 3 points to find the slopes.  If less than 4 points are
*     supplied, straight lines are fitted.
* 
*  2) The sequence xx(1), xx(2), ... xx(nn) must be strictly increasing.
*
*  Called:
*     None
*
*  Test case:
*     Not provided for this subroutine.
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 June   08 B.E.Stetzler    Added header and copyright
*  2009 August 19 B.E.Stetzler    Capitalized all variables for FORTRAN
*                                 77 compatibility
*  2009 August 26 B.E.Stetzler    Used IMPLICIT NONE and defined all variables
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I,J,N,N1,NN,NMAX
      PARAMETER (NMAX = 20)
      REAL A,C,Q,Q1,QN,S,U,X,U1,U2,X1,X2
      DIMENSION X(NMAX),U(NMAX),S(NMAX),A(NMAX)

      Q(U1,X1,U2,X2)=(U1/X1**2-U2/X2**2)/(1.0/X1-1.0/X2)
      

      N = IABS(NN)

      IF (N.LE.3) THEN

*  series too short for cubic spline - use straight lines.
         DO I=1,N
            S(I)=0.0
         ENDDO
         RETURN
      ENDIF

      Q1=Q(U(2)-U(1),X(2)-X(1),U(3)-U(1),X(3)-X(1))
      QN=Q(U(N-1)-U(N),X(N-1)-X(N),U(N-2)-U(N),X(N-2)-X(N))

      IF (NN.LE.0) THEN
         Q1=S(1)
         QN=S(2)
      ENDIF

      S(1)=6.0*((U(2)-U(1))/(X(2)-X(1)) - Q1)
      N1= N - 1

      DO I=2,N1
         S(I)= (U(I-1)/(X(I)-X(I-1)) - U(I)*(1.0/(X(I)-X(I-1))
     .   + 1.0/(X(I+1)-X(I))) + U(I+1)/(X(I+1)-X(I)))*6.0
      ENDDO

      S(N)=6.0*(QN + (U(N1)-U(N))/(X(N)-X(N1)))
      A(1)=2.0*(X(2)-X(1))
      A(2)=1.5*(X(2)-X(1)) + 2.0*(X(3)-X(2))
      S(2)=S(2) - 0.5*S(1)

      DO I=3,N1
         C=(X(I)-X(I-1))/A(I-1)
         A(I)=2.0*(X(I+1)-X(I-1)) - C*(X(I)-X(I-1))
         S(I)=S(I) - C*S(I-1)
      ENDDO

      C=(X(N)-X(N1))/A(N1)
      A(N)=(2.0-C)*(X(N)-X(N1))
      S(N)=S(N) - C*S(N1)

*  Back substitute
      S(N)= S(N)/A(N)

      DO J=1,N1
         I=N-J
         S(I) =(S(I) - (X(I+1)-X(I))*S(I+1))/A(I)
      ENDDO
      RETURN

* Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END

      SUBROUTINE SHELLS (X,K,N)
*+
*  - - - - - - - - -
*   S H E L L S
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The subroutine sorts an array x, of length n, sorting upward,
*  and returns an array k which may be used to key another array
*  to the sorted pattern (i.e., if we had an array f to which x
*  corresponded before sorting, then after calling SHELLS,
*  f(k(1)) will be the element of f corresponding to the
*  smallest x, f(k(2)) the next smallest, and so on).
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model	
* 
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     x              d      array to be sorted (Note 1)
*     n              i      length of the input array x
*
*  Returned:
*     k              i      sorted array that may be used to key another 
*                           array
*  Notes:
*
*  1) See the subroutine ADMINT.F header comments for detailed information.
* 
*  Called:
*     None
*
*  Test case:
*     Not provided for this subroutine.
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  1982 December 29              Revised so that array k is sorted in turn
*                                
*  2009 June   05 B.E. Stetzler    Added header and copyright
*  2009 August 19 B.E. Stetzler    Capitalized all variables for FORTRAN
*                                  77 compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER I,IGAP,IEX,IK,IMAX,IPL,J,K,L,N
      REAL SV,X
      DIMENSION X(N),K(N)

      IGAP = N

      DO 1 I = 1,N
 1    K(I) = I
 5    IF(IGAP.LE.1) GO TO 25

      IGAP = IGAP/2
      IMAX = N - IGAP
 10   IEX = 0
      DO 20 I = 1,IMAX
      IPL = I + IGAP
      IF(X(I).LE.X(IPL)) GO TO 20
      SV = X(I)
      IK = K(I)
      X(I) = X(IPL)
      K(I) = K(IPL)
      X(IPL) = SV
      K(IPL) = IK
      IEX = IEX + 1
 20   CONTINUE

      IF(IEX.GT.0) GO TO 10
      GO TO 5

*  Now sort k's (for identical values of x, if any)

 25   J = 1
 30   IF(J.GE.N) RETURN
      IF(X(J).EQ.X(J+1)) GO TO 33
      J = J + 1
      GO TO 30
*  Have at least two x's with the same value. See how long this is true
 33   L = J
 35   IF(X(L).NE.X(L+1)) GO TO 38
      L = L + 1
      IF(L.LT.N) GO TO 35
*  j and l are the indices within which x(i) does not change - sort k
 38   IGAP = L - J + 1
 40   IF(IGAP.LE.1) J = L + 1
      IF(IGAP.LE.1) GO TO 30

      IGAP = IGAP/2
      IMAX = L-J+1 - IGAP
 45   IEX = 0

      DO 50 I=1,IMAX
      IPL = I + IGAP + J - 1
      IF(K(I+J-1).LE.K(IPL)) GO TO 50
      IK = K(I+J-1)
      K(I+J-1) = K(IPL)
      K(IPL) = IK
      IEX = IEX + 1
 50   CONTINUE
      IF(IEX.GT.0) GO TO 45
      GO TO 40

* Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
      
      INTEGER FUNCTION MDAY (IY, M)
*+
*  - - - - - - - - - - -
*   M D A Y 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This function finds the day number of days before start of month m,
*  of year iy, in Gregorian intercalation.  
* 
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     iy           i      a year
*     m            i      a month
*
*  Returned:
*     mday         i      day number of day before start of a month
*
*  Notes:
*
*  1)  This function needs to test for a leap year. 
*
*  Called:
*     None
*
*  Test case:  This is a support function of the main program HARDISP.F.
*     given input: iy = 2009
*                   m = 5 
*     expected output: mday = 120 
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 April 22 B.E. Stetzler Initial standardization of function
*                              and provided a test case 
*  2009 July  29 B.E. Stetzler Capitalized all variables for FORTRAN 77
*                              compatibility 
*-----------------------------------------------------------------------
      
      IMPLICIT NONE

      INTEGER IY, LEAP, M

      LEAP = 1 - (MOD(IY,4)+3)/4
      IF(MOD(IY,100).EQ.0.AND.MOD(IY,400).NE.0) LEAP=0
      MDAY = MOD((367*(M-2-12*((M-14)/12)))/12+29,365) + LEAP*((9+M)/12)

      RETURN

* Finished.
  
*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END

      INTEGER FUNCTION LEAP (IY)
*+
*  - - - - - - - - - - -
*   L E A P 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This function determines whether a given integer year is a leap
*  year.
* 
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     iy           i      a year (Note 1)
*
*  Returned:
*     0            i      if year is not a leap year
*     1            i      if year is a leap year
*
*  Notes:
*
*  1) The year is a Gregorian year. 
*
*  Called:
*     None
*
*  Test case:  This is a support function of the main program HARDISP.F.
*     given input: IY = 2009
*     
*     expected output: 0
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 April  20 B.E.Stetzler  Initial standardization of function
*                               and provided a test case 
*  2009 August 19 B.E.Stetzler  Capitalized all variables for FORTRAN
*                               77 compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER IY

      LEAP = 1 - (MOD(IY,4)+3)/4

      IF(MOD(IY,100).EQ.0.AND.MOD(IY,400).NE.0) LEAP=0

      RETURN

* Finished.
  
*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END

      INTEGER FUNCTION JULDAT (IT)
*+
*  - - - - - - - - - - -
*   J U L D A T
*  - - - - - - - - - - -
*
*  This function is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This function converts a Gregorian date to a Julian date.
* 
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     it           i      a Gregorian date (Note 1)
*
*  Returned:
*     juldat       i      a Julian date (Note 2) 
*
*  Notes:
*
*  1)  The format of the Gregorian date should be yyyy-mm-dd. 
*  2)  The date is valid for all positive values.
*
*  Called:
*     None
*
*  Test case:  This is a support routine of the main program HARDISP.F.
*     given input: it(1) = 2008
*                  it(2) = 12
*                  it(3) = 12
*     expected output: juldat = 2454813
*
*  References:
*
*     Explanatory Supplement American Ephemeris & Nautical Almanac
*     (cf Comm CACM, 11, 657 (1968) and 15, 918 (1972)), p. 604
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 April  22 B.E.Stetzler  Initial standardization of function
*                               and provided a test case 
*  2009 August 19 B.E.Stetzler  Capitalized all variables for FORTRAN
*                               77 compatibility
*-----------------------------------------------------------------------
      
      IMPLICIT NONE

      INTEGER IT
      DIMENSION IT(*)

      JULDAT = (1461*(IT(1)+4800+(IT(2)-14)/12))/4
     .       + (367*(IT(2)-2-12*((IT(2)-14)/12)))/12
     .       - (3*((IT(1)+4900+(IT(2)-14)/12)/100))/4+IT(3)-32075

      RETURN

* Finished.
  
*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END

      REAL FUNCTION EVAL (Y,NN,X,U,S)
*+
*  - - - - - - - - - - -
*   E V A L 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This function performs cubic spline interpolation of a given function
*  sampled at unequally spaced intervals.  The subroutine SPLINE needs
*  to be called beforehand to set up the array s.
* 
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     y            d     the coordinate at which a function value is
*                        desired (Note 1)
*     nn           i     number of samples of the original function
*     x            d     array containing sample coordinates x(1),x(2),...
*                        x(nn) (Note 2)
*     s            d     array containing the 2nd derivatives at the sample
*                        points (Note 3)
*
*  Returned:
*     u            d     array containing samples of a function at the
*                        coordinates x(1),x(2),...x(nn)
*
*  Notes:
*
*  1) If y falls outside the range (x(1),x(nn)), the value at the nearest
*     endpoint of the series is used.
*
*  2) The sequence x(1),x(2),...x(nn) must be strictly increasing.
*
*  3) This array is found by the subroutine SPLINE, which must be called
*     once before beginning this interpolation.
*
*  Called:
*     None
*
*  Test case:
*     
*  Not provided for this function.  This is a support routine of the main
*  program HARDISP.F.
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 June   03 B.E.Stetzler Initial standardization of function
*  2009 August 19 B.E.Stetzler Capitalized all variables for FORTRAN 77
*                               compatibility
*-----------------------------------------------------------------------
      
      IMPLICIT NONE
      INTEGER K,K1,K2,NN
      REAL DELI,DK,DY,DY1,F1,F2,F3,FF1,FF2,S,U,X,Y
      DIMENSION X(*),U(*),S(*)

      NN = IABS(NN)

*     If y is out of range, substitute endpoint values

      IF (Y.LE.X(1)) THEN
         EVAL=U(1)
         RETURN
      ENDIF

      IF (Y.GE.X(NN)) THEN
         EVAL=U(NN)
         RETURN
      ENDIF

*    Locate interval (x(k1),x(k2)) which contains y
      DO 100 K=2,NN
         IF(X(K-1).LT.Y.AND.X(K).GE.Y) THEN
           K1=K-1
           K2=K
         ENDIF
100   CONTINUE

*    Evaluate and then interpolate
      DY=X(K2)-Y
      DY1=Y-X(K1)
      DK=X(K2)-X(K1)
      DELI=1.0D0/(6.0D0*DK)
      FF1=S(K1)*DY*DY*DY
      FF2=S(K2)*DY1*DY1*DY1
      F1=(FF1+FF2)*DELI
      F2=DY1*((U(K2)/DK)-(S(K2)*DK)/6.0D0)
      F3= DY*((U(K1)/DK)-(S(K1)*DK)/6.0D0)
      EVAL=F1+F2+F3
      RETURN

* Finished.
  
*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END

      SUBROUTINE ETUTC (YEAR, DELTA)
*+
*  - - - - - - - - -
*   E T U T C
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to compute the difference, delta,
*  between Epheremis Time (ET) and Coordinated Universal Time (UTC).
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model	
* 
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     year           d      decimal year (Note 1)
*
*  Returned:
*     delta          d      ET - UTC (Note 2)
*
*     :------------------------------------------:
*     :                                          :
*     :                 IMPORTANT                :
*     :                                          :
*     :  A new version of this routine must be   :
*     :  produced whenever a new leap second is  :
*     :  announced.  There are three items to    :
*     :  change on each such occasion:           :
*     :                                          :
*     :  1) Update the nstep variable            :
*     :  2) Update the arrays st and si          :                              
*     :  3) Change date of latest leap second    :
*     :                                          :
*     :  Latest leap second:  2012 June 31   :
*     :                                          :
*     :__________________________________________:
*
*  Notes:
*
*  1) This subroutine is valid only from 1700.-until next leap second.
*     Currently, this is up to 2009.5.
* 
*  2) The expression used in given in seconds.
* 
*  3) Leap second table in GAMIT UTC (and UT) is the time most 
*     often used (e.g. in time signals)
*
*  Test case:
*     given input: year = 2007.0 
*
*     expected output: delta = 65.1840000000000 seconds
*
*  References:
*
*     Broucke, R. A., Explanatory Supplement American Ephemeris &
*     Nautical Almanac (cf Comm CACM, 11, 657 (1968) and 15, 918 (1972)),
*     p. 90 (Tables)
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 April 22 B.E. Stetzler    Added header and copyright and
*                                 provided test case
*  2009 August 19 B.E. Stetzler   Capitalized all variables for FORTRAN
*                                 77 compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER I,N,NSTEP
      PARAMETER (NSTEP=25)
      REAL*8 D,DELTA,FRAC,ST,SI,TX,TY,YEAR
      DIMENSION D(142),TX(39),TY(39),ST(NSTEP),SI(NSTEP)

* si gives amount of step, at the times given in st

      DATA SI/25*1./
      DATA ST/1972.5,1973.0,1974.0,1975.0,1976.0,1977.0,1978.0,
     .        1979.0,1980.0,1981.5,1982.5,1983.5,1985.5,1988.0,
     .        1990.0,1991.0,1992.5,1993.5,1994.5,1996.0,1997.5,
     .        1999.0,2006.0,2009.0,2012.5/

      DATA D/ 5.15, 4.64, 5.36, 3.49, 3.27, 2.45, 4.03, 1.76, 3.30,
     .  1.00, 2.42, 0.94, 2.31, 2.27,-0.22, 0.03,-0.05,-0.06,-0.57,
     .  0.03,-0.47, 0.98,-0.86, 2.45, 0.22, 0.37, 2.79, 1.20, 3.52,
     .  1.17, 2.67, 3.06, 2.66, 2.97, 3.28, 3.31, 3.33, 3.23, 3.60,
     .  3.52, 4.27, 2.68, 2.75, 2.67, 1.94, 1.39, 1.66, 0.88, 0.33,
     . -0.17,-1.88,-3.43,-4.05,-5.77,-7.06,-7.36,-7.67,-7.64,-7.93,
     . -7.82,-8.35,-7.91,-8.03,-9.14,-8.18,-7.88,-7.62,-7.17,-8.14,
     . -7.59,-7.17,-7.94,-8.23,-7.88,-7.68,-6.94,-6.89,-7.11,-5.87,
     . -5.04,-3.90,-2.87,-0.58, 0.71, 1.80, 3.08, 4.63, 5.86, 7.21,
     .  8.58,10.50,12.10,12.49,14.41,15.59,15.81,17.52,19.01,18.39,
     . 19.55,20.36,21.01,21.81,21.76,22.35,22.68,22.94,22.93,22.69,
     . 22.94,23.20,23.31,23.63,23.47,23.68,23.62,23.53,23.59,23.99,
     . 23.80,24.20,24.99,24.97,25.72,26.21,26.37,26.89,27.68,28.13,
     . 28.94,29.42,29.66,30.29,30.96,31.09,31.59,32.06,31.82,32.69,
     . 33.05,33.16,33.59/
      DATA TX/61.5,
     .62.     ,62.5     ,63.      ,63.5     ,64.      ,64.5     ,65.   ,
     .65.5    ,66.      ,66.5     ,67.      ,67.5     ,68.      ,68.25 ,
     .68.5    ,68.75    ,69.      ,69.25    ,69.5     ,69.75    ,70.   ,
     .70.25   ,70.5     ,70.75    ,71.      ,71.085   ,71.162   ,71.247,
     .71.329  ,71.414   ,71.496   ,71.581   ,71.666   ,71.748   ,71.833,
     .71.915  ,71.999   ,72.0/
      DATA TY/33.59,
     .34.032  ,34.235   ,34.441   ,34.644   ,34.95    ,35.286   ,35.725,
     .36.16   ,36.498   ,36.968   ,37.444   ,37.913   ,38.39    ,38.526,
     .38.76   ,39.      ,39.238   ,39.472   ,39.707   ,39.946   ,40.185,
     .40.42   ,40.654   ,40.892   ,41.131   ,41.211   ,41.284   ,41.364,
     .41.442  ,41.522   ,41.600   ,41.680   ,41.761   ,41.838   ,41.919,
     .41.996  ,42.184   ,42.184/

*  For the oldest epochs, use approximations

      IF(YEAR.LT.1700.0D0) THEN
        DELTA = 0.0D0
        RETURN
      ENDIF
      IF(YEAR.LT.1785.0D0) THEN
        DELTA = (YEAR-1750.0D0)/5.0D0
        RETURN
      ENDIF
      IF(YEAR.LT.1820.5D0) THEN
        DELTA = 6.0D0
        RETURN
      ENDIF

*  For 1820.5 to 1961.5, data is spaced at yearly intervals

      IF(YEAR.LT.1961.5D0) THEN
         N = YEAR - 1819.5
         FRAC = YEAR - (1819.5 + N)
         DELTA = (D(N+1) - D(N))*FRAC + D(N)
         RETURN
      ENDIF

*  For 1961.5 to 1972.0, interpolate between unequispaced data

      IF(YEAR.LT.1972.0D0) THEN
        DO 150 I = 1,38
           IF(YEAR-1900.0D0.EQ.TX(I)) THEN
              DELTA = TY(I)
              RETURN
           ENDIF
           IF(YEAR-1900.0D0.LT.TX(I)) THEN
              DELTA=TY(I-1) + (TY(I)-TY(I-1))*
     .                ((YEAR-1900.0D0-TX(I-1))/(TX(I)-TX(I-1)))
              RETURN
           ENDIF
150     CONTINUE
      ENDIF

*--------------------------------------------------------------------------*
*   after 1972 et-utc has only step offsets. st is the array of step times,*
*   and si is the step sizes (an added second is +1.)                      *
*--------------------------------------------------------------------------*
      DELTA = 42.184D0
      DO 250 I = 1,NSTEP
         IF(YEAR.GE.ST(I)) DELTA = DELTA + SI(I)
         IF(YEAR.LT.ST(I)) RETURN
250   CONTINUE
      RETURN

* Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
      

C     ******************************************************************
C     **                 SUB-ROUTINES FROM JPL                        **
C     ******************************************************************

      subroutine pleph ( ET, NTARG, NCENT, NAMFIL, RRD )
Cf2py intent(in) ET, NTARG, NCENT, NAMFIL
Cf2py intent(out) RRD
C
C++++++++++++++++++++++++++
C  NOTE : Over the years, different versions of PLEPH have had a fifth argument:
C  sometimes, an error return statement number; sometimes, a logical denoting
C  whether or not the requested date is covered by the ephemeris.  We apologize
C  for this inconsistency; in this present version, we use only the four necessary 
C  arguments and do the testing outside of the subroutine.
C
C
C
C     THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS
C     AND GIVES THE POSITION AND VELOCITY OF THE POINT 'NTARG'
C     WITH RESPECT TO 'NCENT'.
C
C     CALLING SEQUENCE PARAMETERS:
C
C       ET = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION
C            IS WANTED.
C
C       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **
C          THE REASON FOR THIS OPTION IS DISCUSSED IN THE 
C          SUBROUTINE STATE
C
C     NTARG = INTEGER NUMBER OF 'TARGET' POINT.
C
C     NCENT = INTEGER NUMBER OF CENTER POINT.
C
C            THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
C
C                1 = MERCURY           8 = NEPTUNE
C                2 = VENUS             9 = PLUTO
C                3 = EARTH            10 = MOON
C                4 = MARS             11 = SUN
C                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
C                6 = SATURN           13 = EARTH-MOON BARYCENTER
C                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
C                            15 = LIBRATIONS, IF ON EPH FILE
C
C             (IF NUTATIONS ARE WANTED, SET NTARG = 14. FOR LIBRATIONS,
C              SET NTARG = 15. SET NCENT=0.)
C
C      RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY
C            OF POINT 'NTARG' RELATIVE TO 'NCENT'. THE UNITS ARE AU AND
C            AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS
C            PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF
C            RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF
C            RADIANS AND RADIANS/DAY.
C
C            The option is available to have the units in km and km/sec.
C            For this, set km=.true. in the STCOMX common block.
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER NTARG, NCENT
      DIMENSION RRD(6),ET2Z(2),ET2(2),PV(6,13)
      DIMENSION SS(3),CVAL(400),PVSUN(6),zips(2)
      data zips/2*0.d0/

      LOGICAL BSAVE,KM,BARY
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      INTEGER LIST(12),IPT(39),DENUM
      
      CHARACTER*180 NAMFIL


      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT

      COMMON/STCOMX/KM,BARY,PVSUN
      KM = .TRUE.


C     INITIALIZE ET2 FOR 'STATE' AND SET UP COMPONENT COUNT
C
      ET2(1)=ET
      ET2(2)=0.D0
      GO TO 11

C     ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT 
C          (SEE THE DISCUSSION IN THE SUBROUTINE STATE)

      ENTRY DPLEPH(ET2Z,NTARG,NCENT, NAMFIL,RRD)

      ET2(1)=ET2Z(1)
      ET2(2)=ET2Z(2)

  11  IF(FIRST) CALL STATE(zips,list,pv,pnut, NAMFIL)
      FIRST=.FALSE.

  96  IF(NTARG .EQ. NCENT) RETURN

      DO I=1,12
      LIST(I)=0
      ENDDO

C     CHECK FOR NUTATION CALL

      IF(NTARG.NE.14) GO TO 97
        IF(IPT(35).GT.0) THEN
          LIST(11)=2
          CALL STATE(ET2,LIST,PV,RRD, NAMFIL)
          RETURN
        ELSE
          do i=1,4
          rrd(i)=0.d0
          enddo
          WRITE(6,297)
  297     FORMAT(' *****  NO NUTATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

C     CHECK FOR LIBRATIONS

  97  do i=1,6
      rrd(i)=0.d0
      enddo

      IF(NTARG.NE.15) GO TO 98
        IF(IPT(38).GT.0) THEN
          LIST(12)=2
          CALL STATE(ET2,LIST,PV,RRD, NAMFIL)
          DO I=1,6
          RRD(I)=PV(I,11)
          ENDDO
          RETURN
        ELSE
          WRITE(6,298)
  298     FORMAT(' *****  NO LIBRATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

C       FORCE BARYCENTRIC OUTPUT BY 'STATE'

  98  BSAVE=BARY
      BARY=.TRUE.

C       SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL

      DO I=1,2
      K=NTARG
      IF(I .EQ. 2) K=NCENT
      IF(K .LE. 10) LIST(K)=2
      IF(K .EQ. 10) LIST(3)=2
      IF(K .EQ. 3) LIST(10)=2
      IF(K .EQ. 13) LIST(3)=2
      ENDDO

C       MAKE CALL TO STATE

      CALL STATE(ET2,LIST,PV,RRD, NAMFIL)

      IF(NTARG .EQ. 11 .OR. NCENT .EQ. 11) THEN
      DO I=1,6
      PV(I,11)=PVSUN(I)
      ENDDO
      ENDIF

      IF(NTARG .EQ. 12 .OR. NCENT .EQ. 12) THEN
      DO I=1,6
      PV(I,12)=0.D0
      ENDDO
      ENDIF

      IF(NTARG .EQ. 13 .OR. NCENT .EQ. 13) THEN
      DO I=1,6
      PV(I,13)=PV(I,3)
      ENDDO
      ENDIF

      IF(NTARG*NCENT .EQ. 30 .AND. NTARG+NCENT .EQ. 13) THEN
      DO I=1,6
      PV(I,3)=0.D0
      ENDDO
      GO TO 99
      ENDIF

      IF(LIST(3) .EQ. 2) THEN
      DO I=1,6
      PV(I,3)=PV(I,3)-PV(I,10)/(1.D0+EMRAT)
      ENDDO
      ENDIF

      IF(LIST(10) .EQ. 2) THEN
      DO I=1,6
      PV(I,10)=PV(I,3)+PV(I,10)
      ENDDO
      ENDIF

  99  DO I=1,6
      RRD(I)=PV(I,NTARG)-PV(I,NCENT)
      ENDDO

      BARY=BSAVE

      RETURN
      END

C++++++++++++++++++++++++
C
      SUBROUTINE FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)
C
C++++++++++++++++++++++++
C
C  THE SUBROUTINE SETS THE VALUES OF  NRECL, KSIZE, NRFILE, AND NAMFIL.

      SAVE

      CHARACTER*180 NAMFIL

C  *****************************************************************
C  *****************************************************************
C
C  THE PARAMETERS NRECL, NRFILE, AND NAMFIL ARE TO BE SET BY THE USER

C  *****************************************************************

C  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
C  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES

      NRECL=4

C  *****************************************************************

C  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE (DEFAULT: 12)

      NRFILE=12

C  *****************************************************************

C  NAMFIL IS THE EXTERNAL NAME OF THE BINARY EPHEMERIS FILE

C      NAMFIL= 'jpl_eph/JPLEPH.421'

C  *****************************************************************

C  KSIZE must be set by the user according to the ephemeris to be read

C  For  de200, set KSIZE to 1652
C  For  de405, set KSIZE to 2036
C  For  de406, set KSIZE to 1456
C  For  de414, set KSIZE to 2036
C  For  inpop13c, set KSIZE to 1876

      KSIZE = 2036
C      KSIZE = 1876

C  *******************************************************************

      RETURN

      END
      
C+++++++++++++++++++++++++++++++++
C
      SUBROUTINE INTERP(BUF,T,NCF,NCM,NA,IFL,PV)
C
C+++++++++++++++++++++++++++++++++
C
C     THIS SUBROUTINE DIFFERENTIATES AND INTERPOLATES A
C     SET OF CHEBYSHEV COEFFICIENTS TO GIVE POSITION AND VELOCITY
C
C     CALLING SEQUENCE PARAMETERS:
C
C       INPUT:
C
C         BUF   1ST LOCATION OF ARRAY OF D.P. CHEBYSHEV COEFFICIENTS OF POSITION
C
C           T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY
C               COEFFICIENTS AT WHICH INTERPOLATION IS WANTED
C               (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE
C               INTERVAL IN INPUT TIME UNITS.
C
C         NCF   # OF COEFFICIENTS PER COMPONENT
C
C         NCM   # OF COMPONENTS PER SET OF COEFFICIENTS
C
C          NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY
C               (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)
C
C          IFL  INTEGER FLAG: =1 FOR POSITIONS ONLY
C                             =2 FOR POS AND VEL
C
C
C       OUTPUT:
C
C         PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION
C               EXPECTED IS PV(NCM,IFL), DP.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      SAVE
C
      DOUBLE PRECISION BUF(NCF,NCM,*),T(2),PV(NCM,*),PC(18),VC(18)

C
      DATA NP/2/
      DATA NV/3/
      DATA TWOT/0.D0/
      DATA PC(1),PC(2)/1.D0,0.D0/
      DATA VC(2)/1.D0/
C
C       ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET
C       OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME
C       WITHIN THAT SUBINTERVAL.
C
      DNA=DBLE(NA)
      DT1=DINT(T(1))
      TEMP=DNA*T(1)
      L=IDINT(TEMP-DT1)+1

C         TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1)

      TC=2.D0*(DMOD(TEMP,1.D0)+DT1)-1.D0

C       CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,
C       AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.
C       (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE
C       CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)

      IF(TC.NE.PC(2)) THEN
        NP=2
        NV=3
        PC(2)=TC
        TWOT=TC+TC
      ENDIF
C
C       BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED
C       AND ARE STORED IN THE ARRAY 'PC'.
C
      IF(NP.LT.NCF) THEN
        DO 1 I=NP+1,NCF
        PC(I)=TWOT*PC(I-1)-PC(I-2)
    1   CONTINUE
        NP=NCF
      ENDIF
C
C       INTERPOLATE TO GET POSITION FOR EACH COMPONENT
C
      DO 2 I=1,NCM
      PV(I,1)=0.D0
      DO 3 J=NCF,1,-1
      PV(I,1)=PV(I,1)+PC(J)*BUF(J,I,L)
    3 CONTINUE
    2 CONTINUE
      IF(IFL.LE.1) RETURN
C
C       IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH
C       DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.
C
      VFAC=(DNA+DNA)/T(2)
      VC(3)=TWOT+TWOT
      IF(NV.LT.NCF) THEN
        DO 4 I=NV+1,NCF
        VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)
    4   CONTINUE
        NV=NCF
      ENDIF
C
C       INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT
C
      DO 5 I=1,NCM
      PV(I,2)=0.D0
      DO 6 J=NCF,2,-1
      PV(I,2)=PV(I,2)+VC(J)*BUF(J,I,L)
    6 CONTINUE
      PV(I,2)=PV(I,2)*VFAC
    5 CONTINUE
C
      RETURN
C
      END
      
C+++++++++++++++++++++++++
C
      SUBROUTINE SPLIT(TT,FR)
C
C+++++++++++++++++++++++++
C
C     THIS SUBROUTINE BREAKS A D.P. NUMBER INTO A D.P. INTEGER
C     AND A D.P. FRACTIONAL PART.
C
C     CALLING SEQUENCE PARAMETERS:
C
C       TT = D.P. INPUT NUMBER
C
C       FR = D.P. 2-WORD OUTPUT ARRAY.
C            FR(1) CONTAINS INTEGER PART
C            FR(2) CONTAINS FRACTIONAL PART
C
C            FOR NEGATIVE INPUT NUMBERS, FR(1) CONTAINS THE NEXT
C            MORE NEGATIVE INTEGER; FR(2) CONTAINS A POSITIVE FRACTION.
C
C       CALLING SEQUENCE DECLARATIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION FR(2)

C       MAIN ENTRY -- GET INTEGER AND FRACTIONAL PARTS

      FR(1)=DINT(TT)
      FR(2)=TT-FR(1)

      IF(TT.GE.0.D0 .OR. FR(2).EQ.0.D0) RETURN

C       MAKE ADJUSTMENTS FOR NEGATIVE INPUT NUMBER

      FR(1)=FR(1)-1.D0
      FR(2)=FR(2)+1.D0

      RETURN

      END


C++++++++++++++++++++++++++++++++
C
      SUBROUTINE STATE(ET2,LIST,PV,PNUT, NAMFIL)
C
C++++++++++++++++++++++++++++++++
C
C THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS FILE
C
C     CALLING SEQUENCE PARAMETERS:
C
C     INPUT:
C
C         ET2   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION
C               IS WANTED.  ANY COMBINATION OF ET2(1)+ET2(2) WHICH FALLS
C               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.
C
C                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE
C                   ENTIRE EPOCH IN ET2(1) AND SET ET2(2)=0.
C
C                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET ET2(1) =
C                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION
C                   EPOCH AND SET ET2(2) = FRACTIONAL PART OF A DAY
C                   ELAPSED BETWEEN ET2(1) AND EPOCH.
C
C                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET
C                   ET2(1) = SOME FIXED EPOCH, SUCH AS START OF INTEGRATION,
C                   AND ET2(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH.
C
C        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION
C               IS WANTED FOR EACH OF THE BODIES ON THE FILE.
C
C                         LIST(I)=0, NO INTERPOLATION FOR BODY I
C                                =1, POSITION ONLY
C                                =2, POSITION AND VELOCITY
C
C               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:
C
C                         I = 1: MERCURY
C                           = 2: VENUS
C                           = 3: EARTH-MOON BARYCENTER
C                           = 4: MARS
C                           = 5: JUPITER
C                           = 6: SATURN
C                           = 7: URANUS
C                           = 8: NEPTUNE
C                           = 9: PLUTO
C                           =10: GEOCENTRIC MOON
C                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY
C                           =12: LUNAR LIBRATIONS (IF ON FILE)
C
C
C     OUTPUT:
C
C          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED
C               QUANTITIES.  THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS
C               STATE IN THE ARRAY STARTING AT PV(1,I).  (ON ANY GIVEN
C               CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE AFFECTED BY THE
C               FIRST 10 'LIST' ENTRIES (AND BY LIST(12) IF LIBRATIONS ARE
C               ON THE FILE) ARE SET.  THE REST OF THE 'PV' ARRAY
C               IS UNTOUCHED.)  THE ORDER OF COMPONENTS STARTING IN
C               PV(1,I) IS: X,Y,Z,DX,DY,DZ.
C
C               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN
C               EQUATOR AND EQUINOX OF J2000 IF THE DE NUMBER IS 200 OR
C               GREATER; OF B1950 IF THE DE NUMBER IS LESS THAN 200. 
C
C               THE MOON STATE IS ALWAYS GEOCENTRIC; THE OTHER NINE STATES 
C               ARE EITHER HELIOCENTRIC OR SOLAR-SYSTEM BARYCENTRIC, 
C               DEPENDING ON THE SETTING OF COMMON FLAGS (SEE BELOW).
C
C               LUNAR LIBRATIONS, IF ON FILE, ARE PUT INTO PV(K,11) IF
C               LIST(12) IS 1 OR 2.
C
C         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,
C               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF
C               QUANTITIES IN NUT IS:
C
C                        D PSI  (NUTATION IN LONGITUDE)
C                        D EPSILON (NUTATION IN OBLIQUITY)
C                        D PSI DOT
C                        D EPSILON DOT
C
C           *   STATEMENT # FOR ERROR RETURN, IN CASE OF EPOCH OUT OF
C               RANGE OR I/O ERRORS.
C
C
C     COMMON AREA STCOMX:
C
C          KM   LOGICAL FLAG DEFINING PHYSICAL UNITS OF THE OUTPUT
C               STATES. KM = .TRUE., KM AND KM/SEC
C                          = .FALSE., AU AND AU/DAY
C               DEFAULT VALUE = .FALSE.  (KM DETERMINES TIME UNIT
C               FOR NUTATIONS AND LIBRATIONS.  ANGLE UNIT IS ALWAYS RADIANS.)
C
C        BARY   LOGICAL FLAG DEFINING OUTPUT CENTER.
C               ONLY THE 9 PLANETS ARE AFFECTED.
C                        BARY = .TRUE. =\ CENTER IS SOLAR-SYSTEM BARYCENTER
C                             = .FALSE. =\ CENTER IS SUN
C               DEFAULT VALUE = .FALSE.
C
C       PVSUN   DP 6-WORD ARRAY CONTAINING THE BARYCENTRIC POSITION AND
C               VELOCITY OF THE SUN.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      SAVE

      DIMENSION ET2(2),PV(6,12),PNUT(4),T(2),PJD(4),BUF(1500),
     . SS(3),CVAL(400),PVSUN(6)

      INTEGER LIST(12),IPT(3,13)

      LOGICAL FIRST
      DATA FIRST/.TRUE./

      CHARACTER*6 TTL(14,3),CNAM(400)
      CHARACTER*180 NAMFIL

      LOGICAL KM,BARY

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,NUMDE,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL
      COMMON/STCOMX/KM,BARY,PVSUN

C
C       ENTRY POINT - 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE
C
      IF(FIRST) THEN
        FIRST=.FALSE.

C ************************************************************************
C ************************************************************************

C THE USER MUST SELECT ONE OF THE FOLLOWING BY DELETING THE 'C' IN COLUMN 1

C ************************************************************************

C        CALL FSIZER1(NRECL,KSIZE,NRFILE,NAMFIL)
C        CALL FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL)
        CALL FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)

      IF(NRECL .EQ. 0) WRITE(*,*)'  ***** FSIZER IS NOT WORKING *****'

C ************************************************************************
C ************************************************************************

      IRECSZ=NRECL*KSIZE
      NCOEFFS=KSIZE/2

        OPEN(NRFILE,
     *       FILE=NAMFIL,
     *       ACCESS='DIRECT',
     *       FORM='UNFORMATTED',
     *       RECL=IRECSZ,
     *       STATUS='OLD')

      READ(NRFILE,REC=1)TTL,CNAM,SS,NCON,AU,EMRAT,
     . ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3)

      READ(NRFILE,REC=2)CVAL

      NRL=0

      ENDIF


C       ********** MAIN ENTRY POINT **********


      IF(ET2(1) .EQ. 0.D0) RETURN

      S=ET2(1)-.5D0
      CALL SPLIT(S,PJD(1))
      CALL SPLIT(ET2(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)+.5D0
      PJD(2)=PJD(2)+PJD(4)
      CALL SPLIT(PJD(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)

C       ERROR RETURN FOR EPOCH OUT OF RANGE

      IF(PJD(1)+PJD(4).LT.SS(1) .OR. PJD(1)+PJD(4).GT.SS(2)) GO TO 98

C       CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL

      NR=IDINT((PJD(1)-SS(1))/SS(3))+3
      IF(PJD(1).EQ.SS(2)) NR=NR-1

        tmp1 = DBLE(NR-3)*SS(3) + SS(1)
        tmp2 = PJD(1) - tmp1
        T(1) = (tmp2 + PJD(4))/SS(3)

C       READ CORRECT RECORD IF NOT IN CORE

      IF(NR.NE.NRL) THEN
        NRL=NR
        READ(NRFILE,REC=NR,ERR=99)(BUF(K),K=1,NCOEFFS)
      ENDIF

      IF(KM) THEN
      T(2)=SS(3)*86400.D0
      AUFAC=1.D0
      ELSE
      T(2)=SS(3)
      AUFAC=1.D0/AU
      ENDIF

C   INTERPOLATE SSBARY SUN

      CALL INTERP(BUF(IPT(1,11)),T,IPT(2,11),3,IPT(3,11),2,PVSUN)

      DO I=1,6
      PVSUN(I)=PVSUN(I)*AUFAC
      ENDDO

C   CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED

      DO 4 I=1,10
      IF(LIST(I).EQ.0) GO TO 4

      CALL INTERP(BUF(IPT(1,I)),T,IPT(2,I),3,IPT(3,I),
     & LIST(I),PV(1,I))

      DO J=1,6
       IF(I.LE.9 .AND. .NOT.BARY) THEN
       PV(J,I)=PV(J,I)*AUFAC-PVSUN(J)
       ELSE
       PV(J,I)=PV(J,I)*AUFAC
       ENDIF
      ENDDO

   4  CONTINUE

C       DO NUTATIONS IF REQUESTED (AND IF ON FILE)

      IF(LIST(11).GT.0 .AND. IPT(2,12).GT.0)
     * CALL INTERP(BUF(IPT(1,12)),T,IPT(2,12),2,IPT(3,12),
     * LIST(11),PNUT)

C       GET LIBRATIONS IF REQUESTED (AND IF ON FILE)

      IF(LIST(12).GT.0 .AND. IPT(2,13).GT.0)
     * CALL INTERP(BUF(IPT(1,13)),T,IPT(2,13),3,IPT(3,13),
     * LIST(12),PV(1,11))

      RETURN

  98  WRITE(*,198)ET2(1)+ET2(2),SS(1),SS(2)
 198  format(' ***  Requested JED,',f12.2,
     * ' not within ephemeris limits,',2f12.2,'  ***')

      stop

   99 WRITE(*,'(2F12.2,A80)')ET2,'ERROR RETURN IN STATE'

      STOP

      END
C+++++++++++++++++++++++++++++
C
      SUBROUTINE CONST(NAM,VAL,SSS,N,NAMFIL)
C
C+++++++++++++++++++++++++++++
C
C     THIS ENTRY OBTAINS THE CONSTANTS FROM THE EPHEMERIS FILE
C
C     CALLING SEQEUNCE PARAMETERS (ALL OUTPUT):
C
C       NAM = CHARACTER*6 ARRAY OF CONSTANT NAMES
C
C       VAL = D.P. ARRAY OF VALUES OF CONSTANTS
C
C       SSS = D.P. JD START, JD STOP, STEP OF EPHEMERIS
C
C         N = INTEGER NUMBER OF ENTRIES IN 'NAM' AND 'VAL' ARRAYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      SAVE

      CHARACTER*6 NAM(*),TTL(14,3),CNAM(400)
      CHARACTER*180 NAMFIL

      DOUBLE PRECISION VAL(*),SSS(3),SS(3),CVAL(400),zips(2)
      DOUBLE PRECISION xx(99)
      data zips/2*0.d0/

      INTEGER IPT(3,13),DENUM,list(11)
      logical first
      data first/.true./

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL

C  CALL STATE TO INITIALIZE THE EPHEMERIS AND READ IN THE CONSTANTS

      IF(FIRST) CALL STATE(zips,list,xx,xx,NAMFIL)
      first=.false.

      N=NCON

      DO I=1,3
      SSS(I)=SS(I)
      ENDDO

      DO I=1,N
      NAM(I)=CNAM(I)
      VAL(I)=CVAL(I)
      ENDDO

      RETURN

      END
