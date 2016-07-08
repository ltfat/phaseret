
typedef enum
{
    PHASERET_DGT_FILTERBANK_ALG,
    PHASERET_DGT_FACTLONG_ALG
} phaseret_dgt_alg;

int
gla_long(double* s, double* g, int L, int W, int a, int M,
         int iter, complex double* c);

int
force_magnitude(complex double* cin, double* s, int L, complex double* cout);

