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

  type(c_ptr) :: handle = c_null_ptr

  integer(c_size_t) :: free, total

  integer :: i
  double precision :: error
  double precision, parameter :: error_max = 10*epsilon(error)

  call hipblasCheck(hipblasCreate(handle))

  lda        = m;
  size_a     = k * lda;

  ldb        = k;
  size_b     = n * ldb;

  ldc    = m;
  size_c = n * ldc;

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
  write(6,'(A30,":",1X,F8.3,1X,"GB")') 'before data map',free/1073741824.d0
  call flush(6)

  !$omp target data map(ha,hb,hc)

  call hipCheck(hipDeviceSynchronize())
  call hipCheck(hipMemGetInfo(free,total))
  write(6,'(A30,":",1X,F8.3,1X,"GB")') 'after data map',free/1073741824.d0
  call flush(6)

  !$omp target data use_device_ptr(ha,hb,hc)

  call hipblasCheck(hipblasDgemm(handle,transa,transb,m,n,k,alpha,c_loc(ha(1)),lda,c_loc(hb(1)),ldb,beta,c_loc(hc(1)),ldc))

  !$omp end target data

  !$omp end target data

  do i = 1,size_c
     error = abs((hc_exact(i) - hc(i))/hc_exact(i))
     if( error > error_max )then
        write(*,*) i, "FAILED! Error bigger than max! Error = ", error
        call exit(1)
     end if
  end do

  !$omp target exit data map(delete:ha,hb,hc)

  deallocate(ha)
  deallocate(hb)
  deallocate(hc)

  call hipCheck(hipDeviceSynchronize())
  call hipCheck(hipMemGetInfo(free,total))
  write(6,'(A30,":",1X,F8.3,1X,"GB")') 'after map(delete:ha,hb,hc)',free/1073741824.d0
  call flush(6)

  deallocate(hc_exact)

  call hipblasCheck(hipblasDestroy(handle))

  write(6,'(A)') "PASSED!"

end program hip_dgemm

