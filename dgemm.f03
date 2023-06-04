program hip_dgemm

  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipblas

  implicit none

  integer(kind(HIPBLAS_OP_N)), parameter :: transa = HIPBLAS_OP_N, transb = HIPBLAS_OP_N;
  double precision, parameter ::  alpha = 1.1d0, beta = 0.9d0;

  integer(c_int) ::  m = 46340, n = 46340, k = 46340;
  integer(c_int) :: lda, ldb, ldc, size_a, size_b, size_c;

  double precision, allocatable, target, dimension(:) :: ha, hb, hc
  double precision, allocatable, dimension(:) :: hc_exact

  type(c_ptr) :: da = c_null_ptr, db = c_null_ptr, dc = c_null_ptr
  type(c_ptr) :: handle = c_null_ptr

! integer, parameter :: bytes_per_element = 8 !double precision
  integer(c_size_t) :: Nabytes, Nbbytes, Ncbytes
  integer(c_size_t) :: free, total

  integer :: i
  double precision :: error
  double precision, parameter :: error_max = 10*epsilon(error)

  call hipblasCheck(hipblasCreate(handle))

  lda        = m;
  size_a     = k * lda;
  Nabytes = size_a*8_c_size_t

  ldb        = k;
  size_b     = n * ldb;
  Nbbytes = size_b*8_c_size_t

  ldc    = m;
  size_c = n * ldc;
  Ncbytes = size_c*8_c_size_t

  allocate(ha(size_a))
  allocate(hb(size_b))
  allocate(hc(size_c))
  allocate(hc_exact(size_c))

  ! Use these constant matrices so the exact answer is also a
  ! constant matrix and therefore easy to check
  ha(:) = 1.d0
  hb(:) = 2.d0
  hc(:) = 3.d0
  hc_exact = alpha*k*2.d0 + beta*3.d0

  call hipCheck(hipDeviceSynchronize())
  call hipCheck(hipMemGetInfo(free,total))
  write(6,'(A30,":",1X,F8.3,1X,"GB")') 'before hipMalloc',free/1073741824.d0
  call flush(6)

  ! Allocate device memory
  call hipCheck(hipMalloc(da,Nabytes))
  call hipCheck(hipMalloc(db,Nbbytes))
  call hipCheck(hipMalloc(dc,Ncbytes))

  call hipCheck(hipDeviceSynchronize())
  call hipCheck(hipMemGetInfo(free,total))
  write(6,'(A30,":",1X,F8.3,1X,"GB")') 'after hipMalloc',free/1073741824.d0
  call flush(6)

  !Transfer from host to device
  call hipCheck(hipMemcpy(da, c_loc(ha(1)), Nabytes, hipMemcpyHostToDevice))
  call hipCheck(hipMemcpy(db, c_loc(hb(1)), Nbbytes, hipMemcpyHostToDevice))
  call hipCheck(hipMemcpy(dc, c_loc(hc(1)), Ncbytes, hipMemcpyHostToDevice))


  call hipblasCheck(hipblasDgemm(handle,transa,transb,m,n,k,alpha,da,lda,db,ldb,beta,dc,ldc))

  call hipCheck(hipDeviceSynchronize())

  ! Transfer data back to host memory
  call hipCheck(hipMemcpy(c_loc(hc(1)), dc, Ncbytes, hipMemcpyDeviceToHost))


  do i = 1,size_c
     error = abs((hc_exact(i) - hc(i))/hc_exact(i))
     if( error > error_max )then
        write(*,*) "FAILED! Error bigger than max! Error = ", error
        call exit(1)
     end if
  end do

  call hipCheck(hipFree(da))
  call hipCheck(hipFree(db))
  call hipCheck(hipFree(dc))

  call hipCheck(hipDeviceSynchronize())
  call hipCheck(hipMemGetInfo(free,total))
  write(6,'(A30,":",1X,F8.3,1X,"GB")') 'after hipFree',free/1073741824.d0
  call flush(6)

  call hipblasCheck(hipblasDestroy(handle))

  deallocate(ha)
  deallocate(hb)
  deallocate(hc)
  deallocate(hc_exact)

  write(*,*) "PASSED!"

end program hip_dgemm
