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

  integer :: ii, jj, kk, irec, test
  integer :: nrows = 8    
  integer :: ncols = 8
  integer :: n_root_cols, root_col_ind
  integer :: myrow_start, n_myrows, n_mycols, mycols, myrows, max_n_mycols, root_colstart
  integer, allocatable :: ind_col(:)            ! ind_col(ncols) = mpi task that owns the column
  integer, allocatable :: global_ncols(:) 
  integer, allocatable :: tot_data(:,:)         ! cheating...just to pick my_data from
  integer, allocatable :: my_rowdata(:,:)    
  integer, allocatable :: send_buffer(:,:)      ! buffer to fill and reduce
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
  myrow_start = n_myrows*my_id     ! remember to account for the non-zero indexing somewhere
  if (my_id == ntasks - 1) then
    n_myrows = n_myrows+diff       ! give the leftovers to the last task
  endif  

  ! give each MPI task its row data
  allocate(my_rowdata(n_myrows, ncols))
  do ii = 1, n_myrows
    do jj = 1, ncols
      my_rowdata(ii, jj) = tot_data(ii+myrow_start, jj)
    enddo
  enddo

  deallocate(tot_data) ! No cheating!

  ! determine who gets which columns by distributing round-robin
  allocate(global_ncols(ntasks))
  allocate(ind_col(ncols))
  global_ncols = 0
  do ii = 1, ncols
    jj = mod(ii-1, ntasks)                        ! ii-1 to start round robin at task 0
    ind_col(ii) = jj                              ! says which column (ii) belongs to which task (jj)
    global_ncols(jj+1) = global_ncols(jj+1) + 1   ! says how many columns each task owns
  enddo

  ! note: after this works in the conceptual way, then need to think about
  ! the actual row/column efficiency

  max_n_mycols = MAXVAL(global_ncols)
  allocate(send_buffer(nrows, max_n_mycols))
  allocate(rec_buffer(nrows, max_n_mycols))
  send_buffer = 0
  rec_buffer = 0
  bufsize = nrows*max_n_mycols

  do irec = 1, ntasks 
    root = irec-1                             ! reduce to the root
    n_root_cols = global_ncols(irec)          ! how many cols the root needs
    do ii = 1, n_myrows 
      do jj = 1, n_root_cols 
        root_col_ind = irec + (ntasks*(jj-1)) ! this tells you which columns to grab
                                              ! for the root, i.e. grabs the round-robin cols
        send_buffer(myrow_start+ii, jj) = my_rowdata(ii,root_col_ind)
      enddo
    enddo 

   call MPI_REDUCE(send_buffer, rec_buffer, bufsize, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierror)
   send_buffer = 0 ! reset the buffer (this may not be necessary?)
  enddo


  do kk = 1, ntasks
    if (my_id == kk-1) then
      write (*,*) "my id: ", my_id
      do ii = 1, nrows
        write(*,*) (rec_buffer(ii,jj), jj=1,max_n_mycols)
      enddo
    endif
  enddo

  deallocate(my_rowdata)
  deallocate(send_buffer)
  deallocate(rec_buffer)
  deallocate(global_ncols)
  deallocate(ind_col)

call MPI_FINALIZE(ierror)
END PROGRAM
