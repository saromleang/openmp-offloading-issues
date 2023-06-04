# openmp-offloading-issues
Repository for OpenMP Offloading Issues

Issue: Device memory created from an `omp target data` region is not being deallocated after leaving that region. Adding a `target exit data map(delete:...)` statement has no effect either.

Systems: Crusher and Frontier
Compiler: Cray Fortran

Reproducer: 
- dgemm.f03 provides an example where OpenMP Offload is not used and the reported memory management with HIP is shown to work.
- dgemm-offload.f03 provides an example where OpenMP Offload is used and GPU free memory at different parts of the code is printed. This example shows that device memory is still consumed after leaving a target data region.

