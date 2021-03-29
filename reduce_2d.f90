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
  ! the next step will be to do this as a 4D array for extra staging

  integer :: ntasks, my_id 
  integer :: ierror, tag
  integer :: bufsize
  integer :: my_start
  integer :: root 

  integer :: ii, jj, kk, irec
  integer :: nrows = 8    ! maybe change to 4 later for easier matrix size
  integer :: ncols = 8
  integer :: myrow_start, n_myrows, n_mycols, mycols, myrows, max_n_mycols, root_colstart
  integer, allocatable :: tot_data(:,:) ! cheating...just to pick my_data from
  integer, allocatable :: my_rowdata(:,:)    !
  integer, allocatable :: send_buffer(:,:) ! buffer to fill and reduce
  integer, allocatable :: rec_buffer(:,:)
  
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

  deallocate(tot_data) ! No cheating!

  ! determine who gets which columns (contiguous for now, then round-robin):
  ! note that doing contiguous for now just saves the staging step
  ! We can add the staging step in once everything else is working

  n_mycols = ncols / ntasks
  diff = mod(ncols, ntasks)
  mycol_start = n_mycols*my_id
  max_n_mycols = n_mycols+diff
  if (my_id == ntasks-1) then
    n_mycols = max_n_mycols
  endif


  ! note: after this works in the conceptual way, then need to think about
  ! the actual row/column efficiency

  ! for now, just get working with even divisibility
  ! then, change to round robin, then implement the indexing scheme
  ! then can do non-even divisibility
  allocate(send_buffer(nrows, max_n_mycols))
  allocate(rec_buffer(nrows, max_n_mycols))
  send_buffer = 0
  bufsize = nrows*max_n_mycols
  do irec = 1, ntasks ! just do first for now
    root = irec-1 ! reduce to the root
    ! now, need to put my row data in the correct spot in the buffer
    root_colstart = root*n_mycols ! columns start in same place for each mpi task, depending on
                                  ! the root task
    do ii = 1, n_myrows ! let i be rows first. You are sending ALL the row data you have,
                        ! just not for all columns.
      do jj = 1, n_mycols ! we are not sending all columns! this is wrong ! ! number of columns to send let j be columns first
                  ! the number of columns to send depends on the RECEIVING task
        send_buffer(myrow_start+ii, jj) = my_rowdata(ii,jj+root_colstart)
      enddo
    enddo 

    call MPI_REDUCE(send_buffer, rec_buffer, bufsize, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierror)
  enddo ! irec
  do kk = 1, ntasks
    if (my_id == kk-1) then
      write (*,*) "my id: ", my_id
      do ii = 1, nrows
        write(*,*) (rec_buffer(ii,jj), jj=1,max_n_mycols)
      enddo
    endif
  enddo
  ! the size of the send buffer is actually going to depend on which task is receiving
  ! (if a task has more columns than the others...)

  deallocate(my_rowdata)
  deallocate(send_buffer)
  deallocate(rec_buffer)

call MPI_FINALIZE(ierror)
END PROGRAM
