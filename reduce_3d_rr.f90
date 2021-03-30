PROGRAM reduce_3d_rr
include 'mpif.h'
  ! Now, we're working with a 3d array that more closely simulates the structure
  ! of the BGW wavefunction, which is that we have:
  ! tot_data(kpoints, spin, bands)
  ! and we start with each mpi task having all kpoints and all spins for some bands
  ! and want to distribute such that we have all bands and all spins for some kpoints
  ! The only "non-transferable" part will be that kpoints are a "hidden index"
  ! of g-vectors, so will be slightly more complicated when we go to gather
  ! all the data in the staging steps.


  ! variable name changes:
  ! row --> band
  ! column --> kpoint

  integer :: ntasks, my_id 
  integer :: ierror, tag
  integer :: bufsize
  integer :: my_start
  integer :: root 

  integer :: ii, jj, kk, irec, counter
  integer :: nband = 8                          ! data size (can change)
  integer :: nkpt = 8
  integer :: nspin = 2
  integer :: n_root_kpt, root_kpt_ind          ! root task = task the reduce is being sent to
  integer :: myband_start, n_mybands, n_mykpts, mykpts, mybands, max_n_mykpts, root_kptstart

  integer, allocatable :: ind_kpt(:)            ! ind_col(ncols) = mpi task that owns the column
  integer, allocatable :: global_nkpts(:)       ! global_ncols(ntasks) = number of columns each task owns
  integer, allocatable :: tot_data(:,:,:)         ! (nkpt, nspin, nband)
  integer, allocatable :: my_banddata(:,:,:)      ! (nkpt, nspin, n_mybands) 
  integer, allocatable :: send_buffer(:,:,:)      ! buffer to fill and reduce
  integer, allocatable :: rec_buffer(:,:,:)
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id,  ierror)

  ! Generate total data:
  allocate(tot_data(nkpt,  nspin, nband))
  counter = 1
  do ii = 1, nband
    do kk = 1, nspin
      do jj = 1, nkpt
        tot_data(jj,kk,ii) = counter
        counter = counter + 1
      enddo
    enddo
  enddo

  ! Display total data:
  if (my_id == 0) then
    write(*,*) "Matrix of all data for spin = 1:"
    do ii = 1, nkpt
      write(*,*) (tot_data(ii,1,jj), jj=1,nband)
    enddo
    write(*,*) "Matrix of all data for spin = 2:"
    do ii = 1, nkpt
      write(*,*) (tot_data(ii,2,jj), jj=1,nband)
    enddo
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)

  ! determine which task gets which bands (contiguous):
  ! note that this is a bad method of load balancing
  n_mybands = nband / ntasks
  diff = mod(nband, ntasks)
  myband_start = n_mybands*my_id     ! remember to account for the non-zero indexing somewhere
  if (my_id == ntasks - 1) then
    n_mybands = n_mybands+diff       ! give the leftovers to the last task
  endif  

  ! give each MPI task its row data
  allocate(my_banddata(nkpt, nspin, n_mybands))
  do jj = 1, n_mybands
    do kk = 1, nspin
      do ii = 1, nkpt
        my_banddata(ii, kk, jj) = tot_data(ii,kk,jj+myband_start)
      enddo
    enddo
  enddo

  ! Display band data:
  do kk = 1, ntasks
    if (my_id == kk-1) then
      write (*,*) "my id: ", my_id, "my band data for spin 1: "
      do ii = 1, nkpt
        write(*,*) (my_banddata(ii,1,jj), jj=1,n_mybands)
      enddo
      write (*,*) "my id: ", my_id, "my band data for spin 2: "
      do ii = 1, nkpt
        write(*,*) (my_banddata(ii,2,jj), jj=1,n_mybands)
      enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  enddo

  deallocate(tot_data) ! Now we just have the broken up band data available

  ! determine which task gets which kpoints by distributing round-robin
  allocate(global_nkpts(ntasks))
  allocate(ind_kpt(nkpt))
  global_nkpts = 0
  do ii = 1, nkpt
    jj = mod(ii-1, ntasks)                        ! ii-1 to start round robin at task 0
    ind_kpt(ii) = jj                              ! says which kpt (ii) belongs to which task (jj)
    global_nkpts(jj+1) = global_nkpts(jj+1) + 1   ! says how many kpts each task owns
  enddo

  max_n_mykpts = MAXVAL(global_nkpts)
  allocate(send_buffer(max_n_mykpts, nspin, nband))
  allocate(rec_buffer(max_n_mykpts, nspin, nband))
  send_buffer = 0
  rec_buffer = 0
  bufsize = max_n_mykpts*nspin*nband

  ! Loop over MPI tasks and send each task its data
  do irec = 1, ntasks 
    root = irec-1                             ! reduce to the root
    n_root_kpt = global_nkpts(irec)           ! how many kpoints the root needs
    ! Fill in the buffer:
    do ii = 1, n_mybands 
      do kk = 1, nspin
        do jj = 1, n_root_kpt 
          root_kpt_ind = irec + (ntasks*(jj-1)) ! this tells you which kpoint to grab for root
                                                ! grabs the round-robin kpt index belonging to root
          send_buffer(jj, kk, myband_start+ii) = my_banddata(root_kpt_ind, kk, ii)
        enddo
      enddo
    enddo 
    call MPI_REDUCE(send_buffer, rec_buffer, bufsize, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierror)
    send_buffer = 0 ! reset the buffer (this may not be necessary?)
  enddo

  ! Display the results:
  do kk = 1, ntasks
    if (my_id == kk-1) then
      write (*,*) "my id: ", my_id, "my kpoint data for spin 1: "
      do ii = 1, max_n_mykpts
        write(*,*) (rec_buffer(ii,1,jj), jj=1,nband)
      enddo
      write (*,*) "my id: ", my_id, "my kpoint data for spin 2: "
      do ii = 1, max_n_mykpts
        write(*,*) (rec_buffer(ii,2,jj), jj=1,nband)
      enddo
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  enddo

  deallocate(my_banddata)
  deallocate(send_buffer)
  deallocate(rec_buffer)
  deallocate(global_nkpts)
  deallocate(ind_kpt)

call MPI_FINALIZE(ierror)
END PROGRAM
