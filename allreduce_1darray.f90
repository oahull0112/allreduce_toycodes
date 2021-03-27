PROGRAM allreduce_1darray
include 'mpif.h'
  ! Goal of this program is to take an array, fill out task-specific data, and then send to task 1
  ! e.g. say we have 3 mpi tasks, and task 1 has ints 0, 1; 2 has ints 2, 3; 3 has ints 4, 5
  ! we want each task to make an array:
  ! task 1: [0, 1, 0, 0, 0, 0]
  ! task 2: [0, 0, 2, 3, 0, 0]
  ! task 3: [0, 0, 0, 0, 4, 5]
  ! and then send all data to task 1 so task 1 has:
  ! task 1: [0, 1, 2, 3, 4, 5]

  integer :: ntasks, my_id 
  integer :: ierror, tag
  integer :: bufsize
  integer :: my_start
  integer :: my_nvals = 2 ! number of ints each task gets
  integer :: root = 0     ! send to mpi task 0
  integer, allocatable :: my_data(:)    !
  integer, allocatable :: send_buffer(:) ! buffer to fill and reduce
  integer, allocatable :: rec_buffer(:)
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id,  ierror)

  bufsize = my_nvals*ntasks
  allocate(send_buffer(bufsize))
  allocate(rec_buffer(bufsize))
  send_buffer = 0
  rec_buffer = 0
  ! hold temp data:
  allocate(my_data(my_nvals))
  ! each task gets its own integers of data
  do ii = 1, my_nvals
    my_data(ii) = my_id*my_nvals+ii
  enddo

  ! place data in appropriate spot
  my_start = my_id*my_nvals
  do ii = 1, my_nvals
    send_buffer(my_start+ii) = my_data(ii)
  enddo

  call MPI_REDUCE(send_buffer, rec_buffer, bufsize, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierror)

  if (my_id == root) then
    write(*,*) "rec_buffer for root: ", rec_buffer
  endif
  if (my_id == 1) then
    write(*,*) "rec_buffer for non-root: ", rec_buffer
  endif

  deallocate(send_buffer)
  deallocate(my_data)

call MPI_FINALIZE(ierror)
END PROGRAM
