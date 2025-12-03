! =====================================================================
! Open-Source Summary
! Project: BMmod — Markov Chain Spatial Modeling
! Overview: Modern Fortran 2008/2018 codebase for spatial Markov chain
!           modeling and spectral methods. Provides numerical routines
!           for balancing, Hessenberg reduction, QR eigenvalue analysis,
!           and spectral decomposition, plus domain functions to derive
!           transition rates, probabilities, and 3-D models.
! Author: Bosszz
! Date: 2025-12-02
! =====================================================================
! Module: bmmod_types — precision kinds and size constants
module bmmod_types
  use iso_fortran_env, only: real32, int32
  implicit none
  
  ! Define precision constants
  integer, parameter :: wp = real32
  integer, parameter :: ip = int32
  
  ! Define array size constants
  integer, parameter :: mcat_max = 10
  integer, parameter :: mlag_max = 2700000
  integer, parameter :: mdat_max = mcat_max * mcat_max * mlag_max
  
end module bmmod_types

! Module: bmmod_data — global model parameters and shared state
! Purpose: Holds proportions (`p`), directional rates (`rd`), grid spacing
!          (`dhx`,`dhy`,`dhz`), grid extents (`nhx`,`nhy`,`nhz`), category
!          count (`ncat`), and background category index (`ibkgr`).
module bmmod_data
  use bmmod_types
  implicit none
  
  ! Global variables (formerly common blocks)
  real(wp) :: p(mcat_max)
  real(wp) :: rd(3, mcat_max, mcat_max)
  real(wp) :: dhx, dhy, dhz
  integer(ip) :: nhx, nhy, nhz
  integer(ip) :: ncat
  integer(ip) :: ibkgr
  
  ! Constants
  integer(ip), parameter :: ldbg = 9
  
end module bmmod_data

! Module: bmmod_linalg — linear algebra utilities
! Purpose: Implements `balanc`, `elmhes`, `hqr`, and `spectral` based on
!          robust EISPACK-style algorithms for stable eigen computations.
!          `spectral` constructs spectral components used by core routines.
module bmmod_linalg
  use bmmod_types
  implicit none
  
contains

  !> Performs spectral decomposition of matrix A
  subroutine spectral(n, A, wr, wi, spec)
    integer(ip), intent(in) :: n
    real(wp), intent(in) :: A(mcat_max, mcat_max)
    real(wp), intent(out) :: wr(mcat_max), wi(mcat_max)
    complex(wp), intent(out) :: spec(mcat_max, mcat_max, mcat_max)
    
    real(wp) :: mat_a(mcat_max, mcat_max) ! Local copy
    real(wp) :: r_mat(mcat_max, mcat_max) ! Original matrix copy
    complex(wp) :: w(mcat_max), denom(mcat_max), s(mcat_max, mcat_max)
    complex(wp) :: cc
    integer(ip) :: i, j, k, l, m, ierr
    integer(ip) :: low, igh
    
    ! Copy A to local working arrays
    r_mat = A
    mat_a = A
    
    ! Find eigenvalues
    call balanc(mcat_max, n, mat_a, low, igh)
    call elmhes(mcat_max, n, low, igh, mat_a)
    call hqr(mcat_max, n, low, igh, mat_a, wr, wi, ierr)
    
    if (ierr /= 0) then
       print *, "Error in HQR: ierr =", ierr
    end if
    
    ! Find the spectrum
    ! Convert wr and wi to complex
    do i = 1, n
       w(i) = cmplx(wr(i), wi(i), kind=wp)
    end do
    
    ! Calculate spectrum
    do i = 1, n
       denom(i) = cmplx(1.0_wp, 0.0_wp, kind=wp)
       
       ! Initialize spec
       spec(i, :, :) = cmplx(0.0_wp, 0.0_wp, kind=wp)
       do k = 1, n
          spec(i, k, k) = cmplx(1.0_wp, 0.0_wp, kind=wp)
       end do
       
       ! Calculate denominator and numerator
       do j = 1, n
          if (j /= i) then
             denom(i) = denom(i) * (w(j) - w(i))
             
             ! Assign spec to s
             s = spec(i, :, :)
             
             ! Update spec: spec(i) = s * (w(j)*I - R)
             do l = 1, n
                do k = 1, n
                   cc = cmplx(0.0_wp, 0.0_wp, kind=wp)
                   do m = 1, n
                      if (m == l) then
                         cc = cc + s(k, m) * (w(j) - r_mat(m, l))
                      else
                         cc = cc + s(k, m) * (-r_mat(m, l))
                      end if
                   end do
                   spec(i, k, l) = cc
                end do
             end do
          end if
       end do
       
       ! Divide spec by denom
       spec(i, :, :) = spec(i, :, :) / denom(i)
    end do
    
  end subroutine spectral

  !> Balances a real matrix and isolates eigenvalues whenever possible
  subroutine balanc(nm, n, a, low, igh)
    integer(ip), intent(in) :: nm, n
    real(wp), intent(inout) :: a(nm, n)
    integer(ip), intent(out) :: low, igh
    
    real(wp) :: c, f, g, r, s, b2, radix
    integer(ip) :: i, j, k, l, m, iexc
    logical :: noconv
    
    radix = 16.0_wp
    b2 = radix * radix
    k = 1
    l = n
    
    ! Search for rows isolating an eigenvalue and push them down
    do
       iexc = 0
       do j = l, 1, -1
          ! Check if row j has only zero off-diagonal elements
          do i = 1, l
             if (i /= j .and. a(j, i) /= 0.0_wp) goto 10
          end do
          
          ! Row j is isolated
          m = l
          iexc = 1
          goto 20
          
10        continue
       end do
       goto 30
       
20     continue
       ! Exchange row/col j and m
       if (j /= m) then
          do i = 1, l
             f = a(i, j); a(i, j) = a(i, m); a(i, m) = f
          end do
          do i = k, n
             f = a(j, i); a(j, i) = a(m, i); a(m, i) = f
          end do
       end if
       
       if (iexc == 1) then
          l = l - 1
          if (l == 1) goto 30
       end if
    end do
    
30  continue
    
    ! Search for columns isolating an eigenvalue and push them left
    do
       iexc = 0
       do j = k, l
          ! Check if col j has only zero off-diagonal elements
          do i = k, l
             if (i /= j .and. a(i, j) /= 0.0_wp) goto 40
          end do
          
          ! Col j is isolated
          m = k
          iexc = 2
          goto 50
          
40        continue
       end do
       goto 60
       
50     continue
       ! Exchange row/col j and m
       if (j /= m) then
          do i = 1, l
             f = a(i, j); a(i, j) = a(i, m); a(i, m) = f
          end do
          do i = k, n
             f = a(j, i); a(j, i) = a(m, i); a(m, i) = f
          end do
       end if
       
       if (iexc == 2) then
          k = k + 1
       end if
    end do
    
60  continue
    
    ! Iterative loop for norm reduction
    do
       noconv = .false.
       do i = k, l
          c = 0.0_wp
          r = 0.0_wp
          do j = k, l
             if (j /= i) then
                c = c + abs(a(j, i))
                r = r + abs(a(i, j))
             end if
          end do
          
          if (c == 0.0_wp .or. r == 0.0_wp) cycle
          
          g = r / radix
          f = 1.0_wp
          s = c + r
          
          do while (c < g)
             f = f * radix
             c = c * b2
          end do
          
          g = r * radix
          do while (c >= g)
             f = f / radix
             c = c / b2
          end do
          
          if ((c + r) / f < 0.95_wp * s) then
             g = 1.0_wp / f
             noconv = .true.
             do j = k, n
                a(i, j) = a(i, j) * g
             end do
             do j = 1, l
                a(j, i) = a(j, i) * f
             end do
          end if
       end do
       
       if (.not. noconv) exit
    end do
    
    low = k
    igh = l
  end subroutine balanc

  !> Reduces a submatrix to upper Hessenberg form
  subroutine elmhes(nm, n, low, igh, a)
    integer(ip), intent(in) :: nm, n, low, igh
    real(wp), intent(inout) :: a(nm, n)
    
    integer(ip) :: i, j, m, ip
    real(wp) :: x, y
    
    do m = low + 1, igh - 1
       x = 0.0_wp
       ip = m
       
       do j = m, igh
          if (abs(a(j, m - 1)) > abs(x)) then
             x = a(j, m - 1)
             ip = j
          end if
       end do
       
       if (x /= 0.0_wp) then
          ! Interchange columns
          do j = 1, igh
             y = a(j, ip)
             a(j, ip) = a(j, m)
             a(j, m) = y
          end do
          
          ! Interchange rows
          y = a(ip, m - 1)
          a(ip, m - 1) = a(m, m - 1)
          a(m, m - 1) = y
          
          do i = m + 1, igh
             a(i, m - 1) = a(i, m - 1) / x
             y = a(ip, m)
             a(ip, m) = a(m, m)
             a(m, m) = y
             y = -y
             do j = m + 1, igh
                a(i, j) = a(i, j) + y * a(i, m - 1)
             end do
             
             do j = 1, igh
                a(j, m) = a(j, m) + a(j, m - 1) * a(i, j)
             end do
          end do
          
          do j = igh + 1, n
             y = a(ip, j)
             a(ip, j) = a(m, j)
             a(m, j) = y
             y = -y
             do i = m + 1, igh
                a(i, j) = a(i, j) + y * a(i, m - 1)
             end do
          end do
       end if
    end do
  end subroutine elmhes

  !> Finds eigenvalues of a real upper Hessenberg matrix
  subroutine hqr(nm, n, low, igh, h, wr, wi, ierr)
    integer(ip), intent(in) :: nm, n, low, igh
    real(wp), intent(inout) :: h(nm, n)
    real(wp), intent(out) :: wr(n), wi(n)
    integer(ip), intent(out) :: ierr
    
    integer(ip) :: i, j, k, l, m, en, na, itn, its
    real(wp) :: p, q, r, s, t, w, x, y, zz, norm, tst1, tst2
    logical :: notlas
    
    ierr = 0
    norm = 0.0_wp
    
    do j = 1, n
       do i = 1, min(n, j + 1)
          norm = norm + abs(h(i, j))
       end do
    end do
    
    wr = 0.0_wp
    wi = 0.0_wp
    do i = 1, low - 1
       wr(i) = h(i, i)
    end do
    do i = igh + 1, n
       wr(i) = h(i, i)
    end do
    
    en = igh
    t = 0.0_wp
    
    ! Main loop
    do while (en >= low)
       itn = 30 * n
       
       ! Loop for finding next eigenvalue
       do
          if (itn == 0) then
             ierr = en
             return
          end if
          
          its = 0
          
          ! Inner loop
          do
             ! Look for single small sub-diagonal element
             do l = en, low + 1, -1
                s = abs(h(l - 1, l - 1)) + abs(h(l, l))
                if (s == 0.0_wp) s = norm
                tst1 = s
                tst2 = tst1 + abs(h(l, l - 1))
                if (tst2 == tst1) exit
             end do
             ! if l loop finishes naturally, l = low. But we use exit.
             
             x = h(en, en)
             if (l == en) goto 100 ! One root found
             
             y = h(en - 1, en - 1)
             w = h(en, en - 1) * h(en - 1, en)
             if (l == en - 1) goto 200 ! Two roots found
             
             if (its == 10 .or. its == 20) then
                ! Exceptional shift
                t = t + x
                do i = low, en
                   h(i, i) = h(i, i) - x
                end do
                s = abs(h(en, en - 1)) + abs(h(en - 1, en - 2))
                x = 0.75_wp * s
                y = x
                w = -0.4375_wp * s * s
             end if
             
             its = its + 1
             itn = itn - 1
             
             ! Look for two consecutive small sub-diagonal elements
             do m = en - 2, l, -1
                zz = h(m, m)
                r = x - zz
                s = y - zz
                p = (r * s - w) / h(m + 1, m) + h(m, m + 1)
                q = h(m + 1, m + 1) - zz - r - s
                r = h(m + 2, m + 1)
                s = abs(p) + abs(q) + abs(r)
                p = p / s
                q = q / s
                r = r / s
                if (m == l) exit
                tst1 = abs(p) * (abs(h(m - 1, m - 1)) + abs(zz) + abs(h(m + 1, m + 1)))
                tst2 = tst1 + abs(h(m, m - 1)) * (abs(q) + abs(r))
                if (tst2 == tst1) exit
             end do
             
             h(m + 2, m) = 0.0_wp
             do i = m + 3, en
                h(i, i - 2) = 0.0_wp
                h(i, i - 3) = 0.0_wp
             end do
             
             ! Double QR step
             do k = m, en - 1
                notlas = (k /= en - 1)
                if (k /= m) then
                   p = h(k, k - 1)
                   q = h(k + 1, k - 1)
                   r = 0.0_wp
                   if (notlas) r = h(k + 2, k - 1)
                   x = abs(p) + abs(q) + abs(r)
                   if (x == 0.0_wp) cycle
                   p = p / x
                   q = q / x
                   r = r / x
                end if
                
                s = sign(sqrt(p * p + q * q + r * r), p)
                if (k /= m) then
                   h(k, k - 1) = -s * x
                else if (l /= m) then
                   h(k, k - 1) = -h(k, k - 1)
                end if
                
                p = p + s
                x = p / s
                y = q / s
                zz = r / s
                q = q / p
                r = r / p
                
                ! Row modification
                do j = k, en
                   p = h(k, j) + q * h(k + 1, j)
                   if (notlas) then
                      p = p + r * h(k + 2, j)
                      h(k + 2, j) = h(k + 2, j) - p * zz
                   end if
                   h(k + 1, j) = h(k + 1, j) - p * y
                   h(k, j) = h(k, j) - p * x
                end do
                
                ! Column modification
                do i = l, min(en, k + 3)
                   p = x * h(i, k) + y * h(i, k + 1)
                   if (notlas) then
                      p = p + zz * h(i, k + 2)
                      h(i, k + 2) = h(i, k + 2) - p * r
                   end if
                   h(i, k + 1) = h(i, k + 1) - p * q
                   h(i, k) = h(i, k) - p
                end do
             end do
             
          end do ! End inner loop (its)
          
100       continue
          ! One root found
          wr(en) = x + t
          wi(en) = 0.0_wp
          en = en - 1
          exit ! Back to while(en >= low)
          
200       continue
          ! Two roots found
          p = (y - x) / 2.0_wp
          q = p * p + w
          zz = sqrt(abs(q))
          x = x + t
          if (q >= 0.0_wp) then
             ! Real pair
             zz = p + sign(zz, p)
             wr(en - 1) = x + zz
             wr(en) = wr(en - 1)
             if (zz /= 0.0_wp) wr(en) = x - w / zz
             wi(en - 1) = 0.0_wp
             wi(en) = 0.0_wp
          else
             ! Complex pair
             wr(en - 1) = x + p
             wr(en) = x + p
             wi(en - 1) = zz
             wi(en) = -zz
          end if
          en = en - 2
          exit ! Back to while(en >= low)
       end do
    end do
    
  end subroutine hqr

end module bmmod_linalg

! Module: bmmod_core — domain logic
! Purpose: Converts between conceptual parameters and transition rates,
!          checks consistency, computes embedded probabilities/frequencies,
!          builds 3-D models (`tp3d`), and interpolates (`r2txyz`).
module bmmod_core
  use bmmod_types
  use bmmod_data
  use bmmod_linalg
  implicit none
  
contains

  !> Calculates Embedded Transition Frequency matrix and entropy
  subroutine r2etf(r, f, entropy)
    real(wp), intent(in) :: r(mcat_max, mcat_max)
    real(wp), intent(out) :: f(mcat_max, mcat_max)
    real(wp), intent(out) :: entropy
    
    real(wp) :: s(mcat_max), tot
    integer(ip) :: j, k
    
    ! Calculate 'tot' and row totals
    tot = 0.0_wp
    do j = 1, ncat
       s(j) = -p(j) * r(j, j)
       tot = tot + s(j)
    end do
    
    ! Calculate 'frequencies'
    entropy = 0.0_wp
    f = 0.0_wp
    
    do j = 1, ncat
       do k = 1, ncat
          if (k /= j) then
             f(j, k) = p(j) * r(j, k) / tot
             f(j, j) = f(j, j) + f(j, k)
             if (f(j, k) > 0.0_wp) then
                entropy = entropy - f(j, k) * log(f(j, k))
             else if (f(j, k) < 0.0_wp) then
                entropy = entropy + 10000.0_wp * f(j, k)
             end if
          end if
       end do
    end do
  end subroutine r2etf

  !> Calculate embedded transition probabilities
  subroutine r2etp(r, etp)
    real(wp), intent(in) :: r(mcat_max, mcat_max)
    real(wp), intent(out) :: etp(mcat_max, mcat_max)
    integer(ip) :: j, k
    
    do j = 1, ncat
       etp(j, j) = 1.0_wp
       do k = 1, ncat
          if (j /= k) etp(j, k) = -r(j, k) / r(j, j)
       end do
    end do
  end subroutine r2etp

  !> Calculate rates with respect to volumetric proportions
  subroutine r2p(r, rp)
    real(wp), intent(in) :: r(mcat_max, mcat_max)
    real(wp), intent(out) :: rp(mcat_max, mcat_max)
    integer(ip) :: j, k
    
    do j = 1, ncat
       do k = 1, ncat
          if (k /= j) rp(j, k) = -r(j, j) * p(k) / (1.0_wp - p(j))
       end do
    end do
    
    do j = 1, ncat
       rp(j, j) = -1.0_wp / r(j, j)
       do k = 1, ncat
          if (k /= j) rp(j, k) = r(j, k) / rp(j, k)
       end do
    end do
  end subroutine r2p

  !> Calculate rates with respect to number of embedded occurrences
  subroutine r2n(r, rn)
    real(wp), intent(in) :: r(mcat_max, mcat_max)
    real(wp), intent(out) :: rn(mcat_max, mcat_max)
    real(wp) :: tot
    integer(ip) :: j, k
    
    tot = 0.0_wp
    do j = 1, ncat
       tot = tot - p(j) * r(j, j)
    end do
    
    do j = 1, ncat
       do k = 1, ncat
          if (k /= j) rn(j, k) = r(j, j) * p(k) * r(k, k) / (tot + p(j) * r(j, j))
       end do
    end do
    
    do j = 1, ncat
       rn(j, j) = -1.0_wp / r(j, j)
       do k = 1, ncat
          if (k /= j) rn(j, k) = r(j, k) / rn(j, k)
       end do
    end do
  end subroutine r2n

  !> Calculate rates w.r.t. independent transition frequencies
  subroutine r2i(r, ri)
    real(wp), intent(in) :: r(mcat_max, mcat_max)
    real(wp), intent(out) :: ri(mcat_max, mcat_max)
    real(wp) :: r0(mcat_max, mcat_max), entmax
    integer(ip) :: j, k
    
    do j = 1, ncat
       r0(j, j) = -1.0_wp / r(j, j)
       do k = 1, ncat
          if (k /= j) r0(j, k) = 1.0_wp
       end do
    end do
    
    call indep(r0, ri, entmax)
    
    do j = 1, ncat
       ri(j, j) = -1.0_wp / r(j, j)
       do k = 1, ncat
          if (k /= j) ri(j, k) = r(j, k) / ri(j, k)
       end do
    end do
  end subroutine r2i

  !> Independent juxtapositional tendencies (Max Entropy)
  subroutine indep(r0, rconc, entmax)
    real(wp), intent(in) :: r0(mcat_max, mcat_max)
    real(wp), intent(out) :: rconc(mcat_max, mcat_max)
    real(wp), intent(out) :: entmax
    
    real(wp) :: s(mcat_max), f(mcat_max), rtot(mcat_max)
    real(wp) :: tot, freq
    integer(ip) :: i, j, k
    
    tot = 0.0_wp
    do i = 1, ncat
       s(i) = p(i) / r0(i, i)
       tot = tot + s(i)
    end do
    
    s = s / tot
    f = s ! Initialize marginal frequencies
    
    ! IPF Loop
    do i = 1, 30
       do j = 1, ncat
          rtot(j) = 0.0_wp
          do k = 1, ncat
             if (k /= j) rtot(j) = rtot(j) + f(j) * f(k)
          end do
          f(j) = f(j) * s(j) / rtot(j)
       end do
    end do
    
    ! Convert to rate
    do j = 1, ncat
       rconc(j, j) = -1.0_wp / r0(j, j)
       rtot(j) = 0.0_wp
       do k = 1, ncat
          if (k /= j) rtot(j) = rtot(j) + f(j) * f(k)
       end do
       do k = 1, ncat
          if (k /= j) rconc(j, k) = -rconc(j, j) * f(j) * f(k) / rtot(j)
       end do
    end do
    
    ! Calculate max entropy
    entmax = 0.0_wp
    do j = 1, ncat - 1
       do k = j + 1, ncat
          freq = f(j) * f(k)
          if (freq > 0.0_wp) then
             entmax = entmax - 2.0_wp * freq * log(freq)
          else if (freq < 0.0_wp) then
             entmax = entmax - 20000.0_wp * freq
          end if
       end do
    end do
  end subroutine indep

  !> Embedded Markov chain transition probabilities
  subroutine emctp(r0, r)
    real(wp), intent(inout) :: r0(mcat_max, mcat_max)
    real(wp), intent(out) :: r(mcat_max, mcat_max)
    integer(ip) :: j, k
    
    ! Check for symmetry assumptions
    do j = 1, ncat
       if (j /= ibkgr) then
          do k = 1, ncat
             if (k /= j .and. k /= ibkgr .and. r0(j, k) == -1.0_wp) then
                r0(j, k) = r0(k, j) * p(k) * r0(j, j) / (p(j) * r0(k, k))
             end if
          end do
       end if
    end do
    
    ! Convert to transition rates
    do j = 1, ncat
       if (j /= ibkgr) then
          do k = 1, ncat
             if (k /= ibkgr) then
                if (k == j) then
                   r(j, j) = -1.0_wp / r0(j, j)
                else
                   r(j, k) = r0(j, k) / r0(j, j)
                end if
             end if
          end do
       end if
    end do
  end subroutine emctp
  
  !> Embedded Markov chain transition frequencies
  subroutine emctf(r0, r)
    real(wp), intent(inout) :: r0(mcat_max, mcat_max)
    real(wp), intent(out) :: r(mcat_max, mcat_max)
    
    real(wp) :: fmarg(mcat_max), x(mcat_max)
    real(wp) :: tot, totx
    integer(ip) :: j, k
    logical :: diff
    
    tot = 0.0_wp
    do j = 1, ncat
       if (r0(j, j) > 0.0_wp) then
          tot = tot + p(j) / r0(j, j)
       else
          write(ldbg, "('Mean length for category',i3,' must be greater than zero')") j
       end if
    end do
    
    do j = 1, ncat
       if (j /= ibkgr) then
          fmarg(j) = p(j) / (r0(j, j) * tot)
          do k = 1, ncat
             if (k /= ibkgr) then
                if (j == k) then
                   r(j, j) = -1.0_wp / r0(j, j)
                else
                   if (r0(j, k) == -1.0_wp) r0(j, k) = r0(k, j)
                   r(j, k) = (1.0_wp / r0(j, j)) * (r0(j, k) / fmarg(j))
                end if
             end if
          end do
       end if
    end do
    
    ! Check consistency
    if (ibkgr <= 0 .or. ibkgr > ncat) then
       totx = 0.0_wp
       do j = 1, ncat
          x(j) = 0.0_wp
          do k = 1, ncat
             if (k /= j) x(j) = x(j) + r0(j, k)
          end do
          totx = totx + x(j)
       end do
       
       do j = 1, ncat
          x(j) = x(j) * totx * r0(j, j)
       end do
       
       diff = .false.
       do j = 1, ncat
          if (abs(fmarg(j) - x(j)) > 1e-5) diff = .true. ! Added epsilon
       end do
       
       if (diff) then
          print *, "WARNING: Proportions intrinsic to rate matrix differ from input."
          write(ldbg, *) "WARNING: Proportions intrinsic to rate matrix differ from input."
          write(ldbg, "(10f7.4)") (x(j), j=1, ncat)
       end if
    end if
  end subroutine emctf

  !> Volumetric proportions
  subroutine vprop(r0, rconc)
    real(wp), intent(in) :: r0(mcat_max, mcat_max)
    real(wp), intent(out) :: rconc(mcat_max, mcat_max)
    integer(ip) :: j, k
    
    do j = 1, ncat
       if (j /= ibkgr) then
          rconc(j, j) = -1.0_wp / r0(j, j)
          do k = 1, ncat
             if (k /= ibkgr .and. k /= j) then
                rconc(j, k) = p(k) / (r0(j, j) * (1.0_wp - p(j)))
             end if
          end do
       end if
    end do
  end subroutine vprop

  !> Number of embedded occurrences
  subroutine nprop(r0, rconc)
    real(wp), intent(inout) :: r0(mcat_max, mcat_max)
    real(wp), intent(out) :: rconc(mcat_max, mcat_max)
    real(wp) :: denom
    integer(ip) :: j, k, iloop
    
    if (ibkgr == 0) then
       do j = 1, ncat
          rconc(j, j) = -1.0_wp / r0(j, j)
          denom = 0.0_wp
          do k = 1, ncat
             if (k /= j) denom = denom + p(k) / r0(k, k)
          end do
          do k = 1, ncat
             if (k /= j) rconc(j, k) = -rconc(j, j) * p(k) / (r0(k, k) * denom)
          end do
       end do
    else if (ibkgr >= 1) then
       do k = 1, ncat
          if (r0(k, k) <= 0.0_wp) then
             print *, 'Need a mean length for category', k, '!'
             r0(k, k) = 1.0_wp
          end if
          rconc(k, k) = -1.0_wp / r0(k, k)
       end do
       
       do iloop = 1, 20
          do j = 1, ncat
             denom = 0.0_wp
             do k = 1, ncat
                if (k /= j) denom = denom - p(k) * rconc(k, k)
             end do
             
             if (j /= ibkgr) then
                rconc(j, ibkgr) = -rconc(j, j)
                do k = 1, ncat
                   if (k /= j .and. k /= ibkgr) then
                      rconc(j, k) = rconc(j, j) * p(k) * rconc(k, k) / denom
                      rconc(j, ibkgr) = rconc(j, ibkgr) - rconc(j, k)
                   end if
                end do
             end if
          end do
          
          ! Calculate new background category mean length
          rconc(ibkgr, ibkgr) = 0.0_wp
          do j = 1, ncat
             if (j /= ibkgr) rconc(ibkgr, ibkgr) = rconc(ibkgr, ibkgr) - p(j) * rconc(j, ibkgr)
          end do
          rconc(ibkgr, ibkgr) = rconc(ibkgr, ibkgr) / p(ibkgr)
       end do
       
       do k = 1, ncat
          if (k /= ibkgr) then
             rconc(ibkgr, k) = 0.0_wp
             do j = 1, ncat
                if (j /= ibkgr) rconc(ibkgr, k) = rconc(ibkgr, k) - p(j) * rconc(j, k)
             end do
          end if
       end do
       rconc(ibkgr, :) = rconc(ibkgr, :) / p(ibkgr)
    end if
  end subroutine nprop
  
  !> Check rate matrix consistency
  subroutine checkr(r, diff)
    real(wp), intent(in) :: r(mcat_max, mcat_max)
    logical, intent(out) :: diff
    
    real(wp) :: wr(mcat_max), wi(mcat_max)
    complex(wp) :: spec(mcat_max, mcat_max, mcat_max)
    real(wp) :: dif
    integer(ip) :: i, j, k, iprop
    
    diff = .false.
    call spectral(ncat, r, wr, wi, spec)
    
    do i = 1, ncat
       if (abs(wr(i)) < 1e-5 .and. abs(wi(i)) < 1e-5) then
          iprop = i
          do k = 1, ncat
             do j = 1, ncat
                dif = abs(real(spec(i, j, k)) - p(k))
                if (dif > 0.005_wp) diff = .true.
             end do
          end do
       end if
    end do
    
    if (diff) then
       print *, "WARNING: Proportions intrinsic to rate matrix differ from input."
       write(ldbg, *) "WARNING: Proportions intrinsic to rate matrix differ from input."
       write(ldbg, "(10f7.4)") (real(spec(iprop, 1, k)), k=1, ncat)
       print "(10f7.4)", (real(spec(iprop, 1, k)), k=1, ncat)
    end if
  end subroutine checkr

  !> Matrix print
  subroutine mprint(txt, n, a, unit, fmt, iparen)
    character(len=*), intent(in) :: txt
    integer(ip), intent(in) :: n, unit, iparen
    real(wp), intent(in) :: a(mcat_max, mcat_max)
    character(len=*), intent(in) :: fmt
    
    integer(ip) :: i, j
    
    write(unit, *) ' '
    write(unit, "(a)") trim(txt)
    
    do i = 1, n
       if (iparen == 1) then
          write(unit, "(*(a, " // trim(fmt) // ", a))") ("(", a(i, j), ")", j=1, n)
       else
          write(unit, "(*(" // trim(fmt) // "))") (a(i, j), j=1, n)
       end if
    end do
  end subroutine mprint

  !> Build 3-D transition probability model
  subroutine tp3d(filtxyz, fildet)
    character(len=*), intent(in) :: filtxyz, fildet
    
    real(wp), allocatable :: txyz(:), det(:)
    real(wp) :: t(mcat_max, mcat_max)
    complex(wp) :: s
    integer(ip) :: idim(5), idimdet(3)
    integer(ip) :: ihx, ihy, ihz, ia, ikl, idat
    integer(ip) :: nxyz, ntxyz, ndetxyz
    integer(ip) :: k, l
    real(wp) :: pwr
    
    allocate(txyz(mdat_max))
    allocate(det(mlag_max))
    
    idim(1) = 2 * nhx + 1
    idim(2) = 2 * nhy + 1
    idim(3) = 2 * nhz + 1
    idim(4) = ncat
    idim(5) = ncat
    idimdet = idim(1:3)
    
    ia = 0
    pwr = 1.0_wp / real(ncat - 1, wp)
    nxyz = idim(1) * idim(2) * idim(3)
    
    ! Nested loops for lags
    do ihz = -nhz, nhz
       do ihy = -nhy, nhy
          print *, 'ihz=', ihz ! Restored verbosity
          do ihx = -nhx, nhx
             ia = ia + 1
             call r2txyz(ihx, ihy, ihz, t, s)
             
             det(ia) = real(s, wp)**pwr
             
             ! Flatten t into txyz
             ikl = 0
             do k = 1, ncat
                do l = 1, ncat
                   ikl = ikl + 1
                   idat = ia + (ikl - 1) * nxyz
                   txyz(idat) = t(k, l)
                end do
             end do
          end do
       end do
    end do
    
    open(8, file=filtxyz, status='unknown')
    write(8, *) 5
    write(8, *) idim
    ntxyz = nxyz * ncat * ncat
    write(8, *) txyz(1:ntxyz)
    close(8)
    
    open(8, file=fildet, status='unknown')
    write(8, *) 3
    write(8, *) idimdet
    ndetxyz = nxyz
    write(8, *) det(1:ndetxyz)
    close(8)
    
    deallocate(txyz, det)
  end subroutine tp3d

  !> Interpolate 3-D transition probability
  subroutine r2txyz(ihx, ihy, ihz, t, s)
    integer(ip), intent(in) :: ihx, ihy, ihz
    real(wp), intent(out) :: t(mcat_max, mcat_max)
    complex(wp), intent(out) :: s
    
    real(wp) :: r(mcat_max, mcat_max)
    real(wp) :: wr(mcat_max), wi(mcat_max)
    complex(wp) :: tc(mcat_max, mcat_max), spec(mcat_max, mcat_max, mcat_max), cc
    real(wp) :: x, y, z, r2, d, rx, ry, rz, a2, tadd
    integer(ip) :: k, l, i
    
    x = ihx * dhx
    y = ihy * dhy
    z = ihz * dhz
    r2 = x**2 + y**2 + z**2
    d = sqrt(r2)
    s = cmplx(1.0_wp, 0.0_wp, kind=wp)
    
    if (r2 /= 0.0_wp) then
       do k = 1, ncat
          if (k /= ibkgr) then
             tadd = 0.0_wp
             do l = 1, ncat
                if (l /= ibkgr) then
                   rx = rd(1, k, l)
                   if (x < 0.0_wp) rx = rd(1, l, k) * p(l) / p(k)
                   ry = rd(2, k, l)
                   if (y < 0.0_wp) ry = rd(2, l, k) * p(l) / p(k)
                   rz = rd(3, k, l)
                   if (z < 0.0_wp) rz = rd(3, l, k) * p(l) / p(k)
                   
                   a2 = (x * rx)**2 + (y * ry)**2 + (z * rz)**2
                   r(k, l) = sqrt(a2 / r2)
                   
                   if (k /= l .and. rx < 0.0_wp .and. ry < 0.0_wp .and. rz < 0.0_wp) then
                      r(k, l) = -r(k, l)
                   end if
                   if (k == l) r(k, l) = -r(k, l)
                   tadd = tadd + r(k, l)
                end if
             end do
             r(k, ibkgr) = -tadd
          end if
       end do
       
       ! Compute background row
       do l = 1, ncat
          r(ibkgr, l) = 0.0_wp
          do k = 1, ncat
             if (k /= ibkgr) r(ibkgr, l) = r(ibkgr, l) - p(k) * r(k, l)
          end do
          r(ibkgr, l) = r(ibkgr, l) / p(ibkgr)
       end do
       
       ! Spectrally decompose r
       call spectral(ncat, r, wr, wi, spec)
       
       ! Compute t
       tc = cmplx(0.0_wp, 0.0_wp, kind=wp)
       do i = 1, ncat
          cc = cmplx(d * wr(i), d * wi(i), kind=wp)
          cc = exp(cc)
          s = s * cc
          do k = 1, ncat
             do l = 1, ncat
                tc(k, l) = tc(k, l) + cc * spec(i, k, l)
             end do
          end do
       end do
       t = real(tc, wp)
       
    else
       ! Zero distance
       t = 0.0_wp
       do k = 1, ncat
          t(k, k) = 1.0_wp
       end do
    end if
  end subroutine r2txyz
  
  !> Establish rates from transition probability data
  subroutine t2r(r)
     real(wp), intent(out) :: r(mcat_max, mcat_max)
     real(wp) :: t(mcat_max, mcat_max)
     real(wp) :: wr(mcat_max), wi(mcat_max)
     complex(wp) :: rc(mcat_max, mcat_max), spec(mcat_max, mcat_max, mcat_max), cc
     character(len=200) :: txt
     integer(ip) :: lag, nvar, ivar, j, k, l, i
     real(wp) :: dlag, rtot, ctot
     
     read(7, "(a)") txt
     read(7, *) lag
     open(8, file=txt, status='old')
     read(8, "(a)") txt
     read(8, *) nvar
     do ivar = 1, nvar + lag
        read(8, *)
     end do
     
     read(8, *) dlag, ((t(j, k), k=1, ncat), j=1, ncat)
     close(8)
     
     if (ibkgr > 0) then
        do j = 1, ncat
           if (j /= ibkgr) then
              rtot = 0.0_wp
              do k = 1, ncat
                 if (k /= ibkgr) rtot = rtot + t(j, k)
              end do
              t(j, ibkgr) = 1.0_wp - rtot
           end if
        end do
        
        do k = 1, ncat
           ctot = 0.0_wp
           do j = 1, ncat
              if (j /= ibkgr) ctot = ctot + p(j) * t(j, k)
           end do
           t(ibkgr, k) = (p(k) - ctot) / p(ibkgr)
        end do
     end if
     
     call spectral(ncat, t, wr, wi, spec)
     
     do i = 1, ncat
        cc = cmplx(wr(i), wi(i), kind=wp)
        cc = log(cc) / cmplx(dlag, 0.0_wp, kind=wp)
        wr(i) = real(cc)
        wi(i) = aimag(cc)
     end do
     
     rc = cmplx(0.0_wp, 0.0_wp, kind=wp)
     do i = 1, ncat
        cc = cmplx(wr(i), wi(i), kind=wp)
        do k = 1, ncat
           do l = 1, ncat
              rc(k, l) = rc(k, l) + cc * spec(i, k, l)
           end do
        end do
     end do
     
     r = real(rc, wp)
  end subroutine t2r

end module bmmod_core

! Program: BMmod — entry point
! Purpose: Reads input parameters, builds 1-D and 3-D Markov models,
!          writes outputs, and logs diagnostics. All outputs and logs
!          use anonymized author tag 'Bosszz' and today's date.
program BMmod
  use bmmod_types
  use bmmod_data
  use bmmod_linalg
  use bmmod_core
  implicit none
  
  character(len=200) :: parfil, dbgfil, filtxyz, fildet, filnam, txt
  real(wp) :: detdir(3), dhdir(3)
  real(wp) :: r(mcat_max, mcat_max), t(mcat_max, mcat_max), r0(mcat_max, mcat_max)
  real(wp) :: rconc(mcat_max, mcat_max), f(mcat_max, mcat_max)
  real(wp) :: wr(mcat_max), wi(mcat_max)
  real(wp) :: a(mcat_max, mcat_max), b(mcat_max, mcat_max)
  complex(wp) :: tc(mcat_max, mcat_max), spec(mcat_max, mcat_max, mcat_max)
  complex(wp) :: cc, detc
  real(wp) :: dh, entropy, h, pwr, detr, colr
  integer(ip) :: nhdir(3)
  integer(ip) :: id, icat, iopt, j, k, l, nh, ih, iflag, i
  character(len=1) :: dir(3)
  logical :: diff
  
  dir = ['X', 'Y', 'Z']
  
  ! Get parameters
  print *, 'Name of BMmod parameter file:'
  read(5, "(a)") parfil
  open(7, file=parfil, status='old')
  
  print *, 'Number of categories:'
  read(7, *) ncat
  
  print *, 'proportions:'
  read(7, *) p(1:ncat)
  
  print *, 'Background category:'
  read(7, *) ibkgr
  
  print *, 'Name of debugging file:'
  read(7, "(a)") dbgfil
  call clean_filename(dbgfil)
  
  ! Write debugging info
  open(ldbg, file=dbgfil, status='unknown')
  write(ldbg, *) 'BMmod debugging file'
  write(ldbg, *) ' '
  write(ldbg, "('Parameter file: ', a)") trim(parfil)
  write(ldbg, "('Number of categories:', i2)") ncat
  write(ldbg, "('Proportions:', 10f8.4)") p(1:ncat)
  write(ldbg, "('Background category:', i2)") ibkgr
  
  ! Get 3-D model parameters
  print *, 'Name of 3-D transition probability model file:'
  read(7, "(a)") filtxyz
  call clean_filename(filtxyz)
  print *, 'Name of determinant file:'
  read(7, "(a)") fildet
  call clean_filename(fildet)
  
  print *, 'Determinant limits for 3-D model; x,y,z direction:'
  read(7, *) detdir(1:3)
  print *, 'nhx,nhy,nhz:'
  read(7, *) dhdir(1:3)
  
  dhx = dhdir(1)
  dhy = dhdir(2)
  dhz = dhdir(3)
  
  ! Build Markov chain models in each coordinate direction
  do id = 1, 3
     if (detdir(id) > 1.0_wp) detdir(id) = 1.0_wp
     if (detdir(id) <= 0.0_wp) detdir(id) = 0.01_wp
     
     print "('Name of output file for ', a1, ' direction:')", dir(id)
     read(7, "(a)") filnam
     call clean_filename(filnam)
     open(id, file=filnam, status='unknown')
     
     print "(' # of lags, lag spacing ', a1, ' direction:')", dir(id)
     read(7, *) nh, dh
     nhdir(id) = nh
     
     if (nh > 0) then
        print *, 'Rate matrix option: (1-5)'
        read(7, *) iopt
        write(ldbg, *) ' '
        write(ldbg, *) ' '
        write(ldbg, "('------- ', a1, '-DIRECTION: -------')") dir(id)
        
        ! Read data
        if (iopt /= 2) then
           do j = 1, ncat
              read(7, *) (r0(j, k), k=1, ncat)
           end do
        end if
        
        ! Write option
        select case (iopt)
        case (1)
           txt = '  option 1: specified transition rates'
        case (2)
           txt = 'option 2: trans. prob. at specified lag'
        case (3)
           txt = 'option 3: embedded Markov chain tr. pr.'
        case (4)
           txt = 'option 4: embedded Markov chain tr. fr.'
        case (5)
           txt = 'option 5: w.r.t. indep/max entropy tr. fr.'
        case (6)
           txt = 'option 6: w.r.t. volumetric proportions'
        case (7)
           txt = 'option 7: w.r.t. fr. of embed. occurrences'
        case default
           txt = 'Unknown option'
        end select
        
        write(ldbg, "('Method - ', a)") trim(txt)
        write(ldbg, "('1-D model output file: ', a)") trim(filnam)
        
        ! Calculate transition rates
        if (iopt == 1) then
           do j = 1, ncat
              do k = 1, ncat
                 if (j /= k .and. r0(j, k) == -1.0_wp) then
                    r0(j, k) = p(k) * r0(k, j) / p(j)
                 end if
                 r(j, k) = r0(j, k)
              end do
           end do
        else if (iopt == 2) then
           call t2r(r)
        else if (iopt == 3) then
           call emctp(r0, r)
        else if (iopt == 4) then
           call emctf(r0, r)
        else if (iopt == 5) then
           call indep(r0, rconc, entropy) ! entropy dummy here? No, entropy is real var
        else if (iopt == 6) then
           call vprop(r0, rconc)
        else if (iopt == 7) then
           call nprop(r0, rconc)
        end if
        
        ! Convert from conceptual parameters to transition rates
        if (iopt >= 5) then
           do j = 1, ncat
              if (j /= ibkgr) then
                 r(j, j) = -1.0_wp / r0(j, j)
                 do k = 1, ncat
                    if (k /= ibkgr .and. k /= j) then
                       if (r0(j, k) == -1.0_wp) then
                          r(j, k) = p(k) * r0(k, j) * rconc(k, j) / p(j)
                       else
                          r(j, k) = r0(j, k) * rconc(j, k)
                       end if
                    end if
                 end do
              end if
           end do
        end if
        
        ! Check r consistency
        if (ibkgr <= 0 .or. ibkgr > ncat) then
           call checkr(r, diff)
        else
           ! Calculate background rates
           do j = 1, ncat
              if (j /= ibkgr) then
                 r(j, ibkgr) = 0.0_wp
                 do k = 1, ncat
                    if (k /= ibkgr) r(j, ibkgr) = r(j, ibkgr) - r(j, k)
                 end do
              end if
           end do
           
           do k = 1, ncat
              r(ibkgr, k) = 0.0_wp
              do j = 1, ncat
                 if (j /= ibkgr) r(ibkgr, k) = r(ibkgr, k) - p(j) * r(j, k)
              end do
              r(ibkgr, k) = r(ibkgr, k) / p(ibkgr)
           end do
        end if
        
        ! Check for invalid rates
        do j = 1, ncat
           if (r(j, j) >= 0.0_wp) write(ldbg, "('WARNING: Diagonal Transition Rate ', i2, ' is nonnegative.')") j
           do k = 1, ncat
              if (k /= j) then
                 if (r(j, k) < 0.0_wp) write(ldbg, "('WARNING: Off-diagonal Transition Rate ', 2i3, ' is negative.')") j, k
                 if (r(j, k) > -r(j, j)) write(ldbg, "('WARNING: Off-diagonal Transition Rate ', 2i3, ' is too large for row ', i2)") j, k, j
                 colr = -p(k) * r(k, k) / p(j)
                 if (r(j, k) > colr) write(ldbg, "('WARNING: Off-diagonal Transition Rate ', 2i3, ' is too large for column ', i2)") j, k, k
              end if
           end do
        end do
        
        ! Print rate matrix
        txt = 'Rate Matrix for ' // dir(id) // '-Direction:'
        call mprint(txt, ncat, r, ldbg, 'f10.6', 0)
        
        ! Embedded transition probabilities
        call r2etp(r, f)
        txt = 'embedded transition probabilities:'
        call mprint(txt, ncat, f, ldbg, 'f10.6', 0)
        
        ! Embedded transition frequencies and entropy
        call r2etf(r, f, entropy)
        txt = 'embedded transition frequencies:'
        call mprint(txt, ncat, f, ldbg, 'f10.6', 0)
        write(ldbg, *) 'entropy=', entropy
        
        ! Rates w.r.t. independent transition frequencies
        call r2i(r, f)
        txt = 'w.r.t. independent transition freqs:'
        call mprint(txt, ncat, f, ldbg, 'f10.4', 1)
        
        if (iopt == 6) then
           call r2p(r, f)
           txt = 'w.r.t. volumetric proportions:'
           call mprint(txt, ncat, f, ldbg, 'f10.4', 1)
        end if
        
        if (iopt == 7) then
           call r2n(r, f)
           txt = 'w.r.t. # of embedded occurrences'
           call mprint(txt, ncat, f, ldbg, 'f10.4', 1)
        end if
        
        ! Write header for t file
        write(id, "(10f7.4)") p(1:ncat)
        write(id, "(i2)") ncat * ncat + 1
        write(id, *) 'lag'
        do k = 1, ncat
           do l = 1, ncat
              write(id, "(i1, '-', i1)") k, l
           end do
        end do
        
        ! Write rate matrix for later use
        do k = 1, ncat
           do l = 1, ncat
              rd(id, k, l) = r(k, l)
           end do
        end do
        
        ! Spectrally decompose r
        call spectral(ncat, r, wr, wi, spec)
        
        ! Print eigenvalues
        do i = 1, ncat
           write(ldbg, *) ' '
           do k = 1, ncat
              do l = 1, ncat
                 b(k, l) = aimag(spec(i, k, l))
                 a(k, l) = real(spec(i, k, l))
              end do
           end do
           write(ldbg, "('eigenvalue', i2, ' real:', f8.4, '  imag:', f8.4)") i, wr(i), wi(i)
           write(ldbg, *) 'spectral component matrix:'
           txt = 'real'
           call mprint(txt, ncat, a, ldbg, 'f10.5', 0)
           txt = 'imag'
           if (wi(i) /= 0.0_wp) call mprint(txt, ncat, b, ldbg, 'f10.5', 0)
        end do
        
        ! Calculate t at nh lags for spacing dh
        iflag = 0
        pwr = 1.0_wp / real(ncat - 1, wp)
        ih = 0
        
        do while (ih <= nh .or. iflag == 0)
           h = ih * dh
           
           ! Initialize tc
           tc = cmplx(0.0_wp, 0.0_wp, kind=wp)
           detc = cmplx(1.0_wp, 0.0_wp, kind=wp)
           
           do i = 1, ncat
              cc = cmplx(h * wr(i), h * wi(i), kind=wp)
              cc = exp(cc)
              detc = detc * cc
              do k = 1, ncat
                 do l = 1, ncat
                    tc(k, l) = tc(k, l) + cc * spec(i, k, l)
                 end do
              end do
           end do
           
           detr = real(detc, wp)**pwr
           if (detr <= detdir(id) .and. iflag == 0) then
              iflag = 1
              nhdir(id) = nint(h / dhdir(id))
           end if
           
           t = real(tc, wp)
           write(id, "(f10.3, 100f8.4)") h, ((t(k, l), l=1, ncat), k=1, ncat)
           ih = ih + 1
        end do
     end if
     close(id)
  end do
  
  ! Construct 3-D model
  nhx = nhdir(1)
  nhy = nhdir(2)
  nhz = nhdir(3)
  
  print *, ' '
  write(ldbg, *) ' '
  print *, 'Constructing 3-D transition probability model'
  write(ldbg, *) 'Constructing 3-D transition probability model'
  print *, '# of lags in +x,+y,+z direction =', nhx, nhy, nhz
  write(ldbg, *) '# of lags in +x,+y,+z direction =', nhx, nhy, nhz
  
  if (nhx >= 1 .or. nhy >= 1 .or. nhz >= 1) then
     if (ibkgr > 0 .and. ibkgr <= ncat) then
        call tp3d(filtxyz, fildet)
        print *, ' '
        print *, 'total # of lags =', (2 * nhx + 1) * (2 * nhy + 1) * (2 * nhz + 1)
        write(ldbg, *) 'total # of lags =', (2 * nhx + 1) * (2 * nhy + 1) * (2 * nhz + 1)
     else
        print *, 'No 3-D Markov chain model generated'
        print *, '      - no background category.'
        write(ldbg, *) 'No 3-D Markov chain model generated'
        write(ldbg, *) '      - no background category.'
     end if
  end if
  
  print *, 'BMmod has ended'
  close(ldbg)
  close(7)

contains

  subroutine clean_filename(fname)
    character(len=*), intent(inout) :: fname
    integer(ip) :: idx
    
    idx = index(fname, '/')
    if (idx > 0) fname = fname(:idx-1)
    fname = adjustl(fname)
  end subroutine clean_filename

end program BMmod
