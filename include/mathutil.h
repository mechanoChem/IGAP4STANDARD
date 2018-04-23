#if !defined (MATHUTIL)
#define MATHUTIL

inline int factorial(int n) { int val=1; while (n>0) {val*=n--;} return val; }
inline int int_pow(int base, int pow) { int val=1; while (pow>0) { val*=base; pow--; } return val;  }
inline double double_abs(double a) { return ( a>0 ? a:-a); }
inline int int_max(int a, int b)  { return (a>b ? a:b);  }
inline int int_min(int a, int b)  { return (a<b ? a:b);  }

void matinv_dot_vec(double *mat, double *vec, int n);
void matinv(double *mat, int n);

#endif
