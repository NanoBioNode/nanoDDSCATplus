program FEMGeo_Wr

!   Program pip1: Point in Polyhedron
!   Version 0.80125 (2008 January 25) (Subversion:9)
!   Roman Schuh
!   1 Command line parameter: Maximum shape size in dipoles
!
!  The program calls the ivread_wr.f90 subroutine.
!  This routine is based on the ivread.f90 routine by John Burkardt (1999).
!  The program converts various computer graphics files formates
!  (dxf, obj, oogl, smf, vmrl) into the .FEM format needed with TNONAXSYM.
!  (Have a look at the comments in ivread).
!
!  But the focus of the program FEMGeo_Wr is on the wavefront .obj file format.
!  The input .obj file should such be such that it only consists of triangular (!!)
!  surface patches! No free form curves are supported.
!  All dimensions are in microns.
!
!  This .obj file format will also be generated by the SScaTT (superellipsoid
!  scattering tool), which is also included on this CD.
!  FEMGeo has not been tested with the other file formats that can be read by the
!  ivread.f90 routine.
!
!  The Hyperfun program (www.hyperfun.org) is suitable for generation of other
!  particle shapes. For conversion to .obj, visualization and scaling you may use
!  Deep Exploration (www.righthemisphere.com), for grid reduction you may use
!  Rational Reducer Professional (www.rational-reducer.com).
!  To increase the number of faces of a body you can use a divide by four
!  subdivision scheme implemented in the Triangles
!  DOS program (www.geocities.com/Athens/Academy/8764/triangdoc.html).
!  A divide by three or by four scheme is also included in MilkShape-1.5.7.

  !
  !   only particles consisting of one (!) closed surface are considered
  !

  integer, parameter :: face_max = 500000
  integer, parameter :: node_max = 500000
  integer, parameter :: face_order_max = 3
  character(len = 100) :: filein_name, fileout_name

  integer omp_get_num_threads, omp_get_thread_num

  integer node_num, face_num
  real, dimension(:,:), allocatable :: face_point, face_normal, v4
  real, dimension(:), allocatable :: face_area
  integer, dimension(:), allocatable :: face_order
  integer, dimension(:,:), allocatable :: face
  real (kind = 8), dimension(:,:), allocatable :: vv
  real (kind = 8) pp(3)
  logical inside
  integer num_threads, dipole_index, i

  integer NAT_without_offset
  integer, dimension(:,:), allocatable :: IXYZ
  real A1(3), A2(3), DX(3)
  integer JX,JY,JZ,NX2,NY1,NY2,NZ1,NZ2
  real YCM, ZCM

  real xave, xrange, xmax, xmin
  real yave, yrange, ymax, ymin
  real zave, zrange, zmax, zmin
  real maxxyz, xyzscale
  real array(3)
  integer maxpos
  integer NBX, NBY, NBZ
  integer xsh, ysh, zsh, shape_size
  character*80 strafg

  ! Use the default shape size if no size is given
  if (iargc() .ne. 3) then
     print *, 'Usage: pip MAX_DIPOLE INPUT OUTPUT'
     call exit(1)
  else
     call getarg(1, strafg)
     read (unit=strafg, fmt=*) shape_size

     call getarg(2, filein_name)
     call getarg(3, fileout_name)
  end if

  print *, 'Maximum shape size = ', shape_size
  print *, 'Input filename = ', filein_name
  print *, 'Output filename = ', fileout_name

  xsh = shape_size
  ysh = shape_size
  zsh = shape_size

  DX = 1.0

  allocate (face_point(3, face_max))
  allocate (face_normal(3, face_max))
  allocate (face_area(face_max))
  allocate (face_order(face_max))
  allocate (face(face_order_max, face_max))
  allocate (v4(3, node_max))

  face_order = 3

  call ivread_wr(filein_name, face_point, face_normal, face_area, face_num,&
       node_num, v4, face)

  deallocate (face_point)
  deallocate (face_normal)
  deallocate (face_area)

  ! Find the nodal coordinate range
  call cor3_limits(node_max, node_num, v4,&
       xmin, xave, xmax, ymin, yave, ymax, zmin, zave, zmax)

  ! Create new array that contains the list of nodes as doubles
  allocate (vv(3, node_num))
  vv = dble(v4(1:3, 1:node_num))

  deallocate (v4)

  A1 = 0.0
  A2 = 0.0
  A1(1) = 1.0
  A2(2) = 1.0

  xrange = xmax - xmin
  yrange = ymax - ymin
  zrange = zmax - zmin

  maxxyz = max(xrange, yrange, zrange)
  array(1) = xrange
  array(2) = yrange
  array(3) = zrange
  maxpos = maxloc(array, 1)

  if (maxpos == 1) then
     NBX = xsh
     NBY = int(yrange/xrange*xsh)
     NBZ = int(zrange/xrange*xsh)
     xyzscale = xsh/xrange
  endif
  if (maxpos == 2) then
     NBX = int(xrange/yrange*ysh)
     NBY = ysh
     NBZ = int(zrange/yrange*ysh)
     xyzscale = ysh/yrange
  endif
  if(maxpos == 3)then
     NBX = int(xrange/zrange*zsh)
     NBY = int(yrange/zrange*zsh)
     NBZ = zsh
     xyzscale = zsh/zrange
  endif

  if (2 * (NBX/2) .lt. NBX) then
     XCM = 0.0
     NX1 = -NBX/2
     NX2 = NBX/2
  else
     XCM = 0.5
     NX1 = -NBX/2+1
     NX2 = NBX/2
  endif
  if(2*(NBY/2).lt.NBY)then
     YCM=0.
     NY1=-NBY/2
     NY2=NBY/2
  else
     YCM=0.5
     NY1=-NBY/2+1
     NY2=NBY/2
  endif
  if(2*(NBZ/2).lt.NBZ)then
     ZCM=0.
     NZ1=-NBZ/2
     NZ2=NBZ/2
  else
     ZCM=0.5
     NZ1=-NBZ/2+1
     NZ2=NBZ/2
  endif

  allocate (IXYZ(3,NBX*NBY*NBZ))
  NAT_without_offset = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(JX,JY,JZ, pp, inside)

!$OMP BARRIER
  
!$OMP DO SCHEDULE(STATIC)
  do JZ=NZ1,NZ2
     do JY=NY1,NY2
        do JX=NX1,NX2
           pp(1) = 1.0 * JX / xyzscale + (xmax + xmin) / 2.0
           pp(2) = 1.0 * JY / xyzscale + (ymax + ymin) / 2.0
           pp(3) = 1.0 * JZ / xyzscale + (zmax + zmin) / 2.0

           call polyhedron_contains_point_3d (node_num, face_num,&
                & face_order_max, vv, face_order, face, pp, inside)

           if (inside .eqv. .true.) then
!$OMP CRITICAL
              NAT_without_offset = NAT_without_offset + 1
              IXYZ(1, NAT_without_offset) = JX
              IXYZ(2, NAT_without_offset) = JY
              IXYZ(3, NAT_without_offset) = JZ
!$OMP END CRITICAL

           endif
        enddo
     enddo
  enddo
!$OMP END DO

!$OMP END PARALLEL

  deallocate (face_order)
  deallocate (face)

  ! Write the output file
  open (UNIT=12, FILE=fileout_name, STATUS='UNKNOWN')

  ! WARNING::: NBY, NBZ and IXYZ(2,i),IXYZ(3,i) values were flipped
  ! I have not figured out why, but for now I am just flipping them back
  ! during the print statement to the dat file. The data is good, just flipped.

  ! Write the file header
  write (12, FMT=92) NBX, NBZ, NBY, NAT_without_offset, A1, A2, DX

  ! Iterate over the results from each thread, stored in IXYZ
  dipole_index = 1
  do i = 1, NAT_without_offset
     write (12, FMT=93) dipole_index, IXYZ(1,i), IXYZ(3,i), IXYZ(2,i)
     dipole_index = dipole_index + 1
  end do

  close(UNIT=12)

  deallocate(IXYZ)

92 format(' >PIPOBJ: point-in-polyhedron: NBX, NBY, NBZ=',3I4,/,&
       I11,' = NAT',/,&
       3F9.4,' = A_1 vector',/,&
       3F9.4,' = A_2 vector',/,&
       3F9.6,' = lattice spacings (d_x,d_y,d_z)/d',/,&
       ' 0.0 0.0 0.0',/,&
       '         JA    IX    IY    IZ   ICOMP(x,y,z)')
93 format(I11,3I6, '    1 1 1') ! Note that the dipole polarization is fixed

end program FEMGeo_Wr
