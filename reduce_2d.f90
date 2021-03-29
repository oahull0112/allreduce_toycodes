PROGRAM reduce_2d
include 'mpif.h'
  ! Want to get fancier than the 1d array, and now generate *rows* of data
  ! where the root task needs *columns* of data
  ! In this case, say we have an array:
  ! [ 0  1  2  3  4  5  6  7 
  !   9 10 11 12 13 14 15 16
  !  17 18 19 20 21 22 23 24 
  !  25 26 27 28 29 30 31 32 ]
  ! and each mpi task is assigned contiguous rows and round-robin columns
  ! and each task has its rows, but needs its columns
  ! so we gather the row data and send to the appropriate task for its columns
  ! once this is working, it's the biggest piece of the puzzle
  ! the next step will be to do this as a 3D array for extra staging

  integer :: ntasks, my_id 
  integer :: ierror, tag
  integer :: bufsize
  integer :: my_start
  integer :: my_nvals = 2 ! number of ints each task gets ! THIS WILL LIKELY GO AWAY
  integer :: root = 0     ! send to mpi task 0

  integer :: ii, jj, kk
  integer :: nrows = 8    ! maybe change to 4 later for easier matrix size
  integer :: ncols = 8
  integer :: n_myrows, n_mycols, mycols, myrows
  integer, allocatable :: tot_data(:,:) ! cheating...just to pick my_data from
  integer, allocatable :: my_rowdata(:,:)    !
  integer, allocatable :: send_buffer(:) ! buffer to fill and reduce
  integer, allocatable :: rec_buffer(:)
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id,  ierror)

  ! generate total data:
  allocate(tot_data(nrows, ncols))
  do ii = 1, nrows
    do jj = 1, ncols
      tot_data(ii, jj) = ii*jj
    enddo
  enddo

  ! determine who gets which rows (contiguous):
  ! note that this is a bad method of load balancing
  n_myrows = nrows / ntasks
  diff = mod(nrows, ntasks)
  myrow_start = n_myrows*my_id  ! remember to account for the non-zero indexing somewhere
  if (my_id == ntasks - 1) then
    n_myrows = n_myrows+diff ! give the leftovers to the last task
  endif  

  allocate(my_rowdata(n_myrows, ncols))
  do ii = 1, n_myrows
    do jj = 1, ncols
      my_rowdata(ii, jj) = tot_data(ii+myrow_start, jj)
    enddo
  enddo

  do kk = 1, ntasks
    if (my_id == kk-1) then
      write (*,*) "my id: ", my_id
      do ii = 1, n_myrows
        write(*,*) (my_rowdata(ii,jj), jj=1,ncols)
      enddo
    endif
  enddo

  ! determine who gets which columns (contiguous for now, then round-robin):
  ! note that doing contiguous for now just saves the staging step
  ! We can add the staging step in once everything else is working



  ! Now, give each task their row data:



!  bufsize = my_nvals*ntasks
!  allocate(send_buffer(bufsize))
!  allocate(rec_buffer(bufsize))
!  send_buffer = 0
!  rec_buffer = 0
!  ! hold temp data:
!  allocate(my_data(my_nvals))
!  ! each task gets its own integers of data
!  do ii = 1, my_nvals
!    my_data(ii) = my_id*my_nvals+ii
!  enddo
!
!  ! place data in appropriate spot
!  my_start = my_id*my_nvals
!  do ii = 1, my_nvals
!    send_buffer(my_start+ii) = my_data(ii)
!  enddo
!
!  call MPI_REDUCE(send_buffer, rec_buffer, bufsize, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!
!  if (my_id == root) then
!    write(*,*) "rec_buffer for root: ", rec_buffer
!  endif
!  if (my_id == 1) then
!    write(*,*) "rec_buffer for non-root: ", rec_buffer
!  endif
!
!  deallocate(send_buffer)
!  deallocate(my_data)

call MPI_FINALIZE(ierror)
END PROGRAM
