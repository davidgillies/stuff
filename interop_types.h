#ifdef INTEROP_TYPES
use ISO_C_BINDING
#define TYPE type, bind(C)
#define REAL real(c_float)
#define INTEGER integer(c_int)
#define LOGICAL logical(c_bool)
#define CHARACTER(LEN) character(kind=c_char, LEN)
#endif
