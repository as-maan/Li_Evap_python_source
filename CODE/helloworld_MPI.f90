program mpisimple

	implicit none
	integer ierr, rank, nprocs, resultlen, src_val, elem, elem_r, n, nlines, io
	integer, parameter :: dg=99, rd=98
	include 'mpif.h'
! shell elements
	real, dimension(:,:), allocatable :: shells_x
	real, dimension(:,:), allocatable :: shells_y
	real, dimension(:,:), allocatable :: shells_z
	real, dimension(:,:), allocatable :: shells_norm
! obstruction elements - rukawat = sanskrit(obstruction in the way of)
	real, dimension(:,:), allocatable :: ruk_x
	real, dimension(:,:), allocatable :: ruk_y
	real, dimension(:,:), allocatable :: ruk_z
	real, dimension(:,:), allocatable :: ruk_norm
! source 
	real, dimension(3) :: src_NO, src_FG, src
! Get hostname

	integer :: i
	character(len=500) :: fname,dg_fname
	character (len=8) :: hname
	INTEGER DATE_TIME (8)
	CHARACTER (LEN = 12) REAL_CLOCK (3)
!----------------------------------------------------------------------------------
!				Initiate Run
!----------------------------------------------------------------------------------	
	CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                    REAL_CLOCK (3), DATE_TIME)

	call mpi_init(ierr)
	call mpi_get_processor_name(hname, resultlen, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	
	src_val = 1
	src_NO = 25.4*(/-3.3601,-16.8923,0.3908/)
	src_FG = 25.4*(/3.6954,18.5778,0.3908/)  
! set 1 for NO and 2 for FG
	if (rank == 0) then
		write(6,*) 'Hostname - ', hname
    		write(6,*) 'Run Commenced at -', DATE_TIME
		write(6,*) 'Running on ', nprocs, ' processors'
		if (src_val .EQ. 1) then
			write(6,*) 'Source - ', 'NO'
		elseif (src_val .EQ. 2) then
			write(6,*) 'Source - ', 'FG'
		end if 
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
!----------------------------------------------------------------------------------
!				Read Geometry
!----------------------------------------------------------------------------------
! All processes read all geometry
! start with reading shell elements
	! read x coords of the shell
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.x.DAT'
		
	open (rd, file = fname, status = 'old')
	
	! for the first file in shell elements, find the total number of elements
	nlines = 0
	do
		read(rd,*, END=11) 
		nlines = nlines + 1
	end do
	11 close (rd)
	!write(6,*) 'nlines - ', nlines	
	allocate (shells_x(nlines,3))
	open (rd, file = fname, status = 'old')
	do i = 1,nlines
		read(rd,*) shells_x(i,1), shells_x(i,2), shells_x(i,3)
	end do 
	close(rd)
	
	! read y coords of the shell
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.y.DAT'
	open (rd, file = fname, status = 'old')
		allocate (shells_y(nlines,3))
		do i = 1,nlines
			read(rd,*) shells_y(i,1), shells_y(i,2), shells_y(i,3)
		end do 
	close(rd)
	! read z coords of the shell
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.z.DAT'
	open (rd, file = fname, status = 'old')
		allocate (shells_z(nlines,3))
		do i = 1,nlines
			read(rd,*) shells_z(i,1), shells_z(i,2), shells_z(i,3)
		end do 
	close(rd)
	! read normals of the shell elements
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/shells.norm.DAT'
	open (rd, file = fname, status = 'old')
		allocate (shells_norm(nlines,3))
		do i = 1,nlines
			read(rd,*) shells_norm(i,1), shells_norm(i,2), shells_norm(i,3)
		end do
	close(rd)
	 	
! now read obstruction elements
	! read x coords of the obstruction
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.x.DAT'
	open (rd, file = fname, status = 'old')
	nlines = 0		
		do
			READ (rd,*, END=10) 
			nlines = nlines + 1 	
		end do
		10 close (rd)
		allocate (ruk_x(nlines,3))
		open (rd, file = fname, status = 'old')
		do i = 1,nlines
			read(rd,*) ruk_x(i,1), ruk_x(i,2), ruk_x(i,3)
		end do
		close(rd)
	
	! read y coords of the obstruction
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.y.DAT'
	open (rd, file = fname, status = 'old', ACTION='READ')
		allocate (ruk_y(nlines,3))
		do i = 1,nlines
!			read(6,*) ruk_y(i,:)
			read(rd,*) ruk_y(i,1), ruk_y(i,2), ruk_y(i,3)
!			write(dg,*) ruk_y(i,:)
		end do
	close(rd) 

	! read z coords of the obstruction
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.z.DAT'
	open (rd, file = fname, status = 'old')
		allocate (ruk_z(nlines,3))
		do i = 1,nlines
			read(rd,*) ruk_z(i,1), ruk_z(i,2), ruk_z(i,3)
		end do
	close(rd) 
	
	! read normals of the obstruction elements
	fname = '/media/amaan/Storage/Work/PPPL/Analysis/Li_Evap/shells/fortran_implementation/CAD_Source/Geom/ruk.norm.DAT'
	open (rd, file = fname, status = 'old')
		allocate (ruk_norm(nlines,3))
		do i = 1,nlines
			read(rd,*) ruk_norm(i,1), ruk_norm(i,2), ruk_norm(i,3)
		end do
	close(rd) 	

! Announce that all data has been read
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if (rank == 0) then
		write(6,*) 'All read operations completed'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	
!----------------------------------------------------------------------------------
!				setup pre-run parms
!----------------------------------------------------------------------------------

	call get_src(src_val,src,src_NO,src_FG)
	elem = size(shells_norm,1)
	elem_r = size(ruk_norm,1)
	if (rank == 0) then
		write(6,*) 'Number of elements', elem
		write(6,*) 'Number of obstructions', elem_r
!	else if (rank .NE. 0) then
!		write(6,*) 'Rank - ', rank, 'has ', elem, 'elements'
!		write(6,*) 'Rank - ', rank, 'has ', elem_r, 'obstructions'
		
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	n = elem/nprocs
	write(6,*) size(shells_x,1)
	write(6,*) 'Element 0 x - ', shells_x(1,:)
	write(6,*) 'Element 1 x - ', shells_x(2,:)

	call calc_OBS_para(elem,elem_r,n,rank,nprocs,shells_x,shells_y,shells_z,shells_norm,ruk_x,ruk_y,ruk_z,ruk_norm,src)	

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	write(6,*) "Obstruction Calc barrier crossed by process - ", rank
	call mpi_finalize(ierr)
	end


!----------------------------------------------------------------------------------
!				main code ends here
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
!				source subroutine
!----------------------------------------------------------------------------------

subroutine get_src(src_val,src,src_NO,src_FG)
	implicit none
	integer, intent(in) :: src_val
	real, dimension(3), intent(in) :: src_NO, src_FG
	real, dimension(3), intent(inout) :: src
	 
	if (src_val .EQ. 1) then
		src = src_NO
	else if (src_val .EQ. 2) then
		src = src_FG
	end if
	end subroutine get_src

!----------------------------------------------------------------------------------
!			main obstruction calc subroutine
!----------------------------------------------------------------------------------
	
subroutine calc_OBS_para(elem,elem_r,n,rank,nprocs,shells_x,shells_y,shells_z,shells_norm,ruk_x,ruk_y,ruk_z,ruk_norm,src)
	implicit none
	integer, intent(in) :: elem, elem_r, n, rank, nprocs
	real, dimension(elem,3), intent(in) :: shells_x, shells_y, shells_z, shells_norm
	real, dimension(elem_r,3), intent(in) :: ruk_x, ruk_y, ruk_z, ruk_norm
	real, dimension(3), intent(in) :: src
	real, dimension(3) :: p1,p2,p3,c_p,N_p,p_s,d,rwat1,rwat2,rwat3,N_rwat,p_Q
	real :: t, theta, t_Q
	integer :: elem_d, n_d, m , j, o
	logical :: OBS
	integer :: ii

	write(6,*) 'diagnostic run - run only 10 elements'
		elem_d = 10
		n_d = elem_d/nprocs

	if (rank .EQ. nprocs-1) then
		do ii = 1,elem_d
			call get_p(p1,p2,p3,shells_x(ii,:),shells_y(ii,:),shells_z(ii,:))
			call get_cp(c_p,shells_x(ii,:),shells_y(ii,:),shells_z(ii,:))
			N_p = shells_norm(ii,:)
			call get_srcvec(p_s,d,t,src,c_p)
			m = 0
			o = 0
			call get_ray_angle(theta,N_p,d)
			if (theta < 90.0) then
				do j = 0,elem_r
					call get_p(rwat1,rwat2,rwat3,ruk_x(j,:),ruk_y(j,:),ruk_z(j,:))
					N_rwat = ruk_norm(j,:)
					call get_rwatQ(p_Q,t_Q,N_rwat,p_s,d,t,rwat1)
					OBS = get_obs(rwat1,rwat2,rwat3,p_Q,t_Q,t,N_rwat,d)
				end do
			end if
			!write(6,*) N_p
		end do
	end if
			


	end subroutine calc_OBS_para

!----------------------------------------------------------------------------------
!			 	get points of the element
!----------------------------------------------------------------------------------
	
subroutine get_p(p1,p2,p3,x,y,z)
	real, dimension(3), intent(in) :: x,y,z
	real, dimension(3), intent(inout) :: p1,p2,p3
	
	p1 = (/x(1),y(1),z(1)/)
	p2 = (/x(2),y(2),z(2)/)
	p3 = (/x(3),y(3),z(3)/)

	end subroutine get_p
!----------------------------------------------------------------------------------
!				get center of the element
!----------------------------------------------------------------------------------

subroutine get_cp(c,x,y,z)
	real, dimension(3), intent(in) :: x,y,z
	real, dimension(3), intent(inout) :: c

	c = (/(sum(x)/3),(sum(y)/3),(sum(z)/3)/)

	end subroutine get_cp
!----------------------------------------------------------------------------------
!				get ray from source to element
!----------------------------------------------------------------------------------
	
subroutine get_srcvec(p_s,d,t,src,c_p)
	real, dimension(3), intent(in) :: src, c_p
	real, dimension(3), intent(inout) :: p_s,d
	real, intent(inout) :: t

	t = norm2(c_p-src)
	write(6,*) t, size(c_p,1), size(src,1)
	d = c_p-src
	d = (/d(1)/t,d(2)/t,d(3)/t/)
	
	end subroutine get_srcvec

!----------------------------------------------------------------------------------
!				get angle between element normal and source
!----------------------------------------------------------------------------------

subroutine get_ray_angle(theta,N_p,d)
	real, dimension(3), intent(in) :: N_p, d
	real, intent(inout) :: theta
	real(16), parameter :: PI_16 = 4 * atan (1.0_16)

	theta = dot_product(N_p,d)/(norm2(N_p)*(norm2(d)))
	theta = (180/pi_16)*acos(theta)	

	end subroutine get_ray_angle

!----------------------------------------------------------------------------------
!				get point Q
!----------------------------------------------------------------------------------

subroutine get_rwatQ(p_Q,t_Q,N_rwat,p_s,d,t,rwat1)
	real, dimension(3), intent(in) :: N_rwat, p_s, d, rwat1
	real, intent(in) :: t
	real, dimension(3), intent(inout) :: p_Q
	real, intent(inout) :: t_Q
	real :: d_scalar	
	
	d_scalar = dot_product(N_rwat,rwat1)
	t_Q = (d_scalar-dot_product(N_rwat,p_s))/(dot_product(N_rwat,d))
	p_Q = p_s + t_Q*d
	
	
	end subroutine get_rwatQ

!----------------------------------------------------------------------------------
!				check if Q lies inside the obstructing element
!----------------------------------------------------------------------------------

FUNCTION get_obs(rwat1,rwat2,rwat3,p_Q,t_Q,t,N_rwat,d)
	real, dimension(3), intent(in) :: rwat1, rwat2,rwat3,p_Q,N_rwat,d
	real, intent(in) :: t_Q,t
	logical, intent(out) :: get_obs






	
	
