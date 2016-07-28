! DO3SE assertion helpers, depends on DO3SE_Util_ml
#define ASSERT(check) call assert((check), "ASSERTION FAILED: check")
#define ASSERT_DEFINED(id) call assert(is_def(id), "id not defined")
#define UNKNOWN_STRING(id) call assert(.false., "unknown id: "//trim((id)))
#define ERROR(msg) call assert(.false., (msg))
! Allocation helpers
#define SAFE_DEALLOC(var) if (allocated(var)) deallocate(var)
#define SAFE_ALLOC_1D(var, n) SAFE_DEALLOC(var); allocate(var(n))
#define SAFE_ALLOC_2D(var, n, m) SAFE_DEALLOC(var); allocate(var(n, m))
