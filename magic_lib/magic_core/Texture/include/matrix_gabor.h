/** \file
 *
 * \brief ??? Gabor filter related data structures and functions
 *
 */

/// Floating point (double) Matrix structure for matrix_gabor
typedef struct MatrixStruct 
{
	double **data;
	int height, width;
} Matrix;

/// Unsigned char Matrix structure for matrix_gabor
typedef struct uc_MatrixStruct 
{
	unsigned char **data;
	int height, width;
} uc_Matrix;

/// Integer Matrix structure for matrix_gabor
typedef struct i_MatrixStruct 
{
	int **data;
	int height, width;
} i_Matrix;

/// Region ?? structure for matrix_gabor
typedef struct RegionStruct
{
	int id;
	int number;
	int *member_i;
	int *member_j;
} Region;

#ifndef PI_CONST
#define PI_CONST
#define PI 3.141592653589793115997963468544185161590576171875
#endif

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

#define ANSI 1

#ifdef ANSI
void nrerror(const char *error_text);
double *dvector(int nl, int nh);
float *fvector(int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_fvector(float *v, int nl, int nh);
float **matrix(int nrl,int nrh,int ncl,int nch);
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
double logBase2(double a);
void nrsort(double *Y, int *I, double *A, int length);
void minimun(double *Y, int *I, double *A, int length);
void Mat_Abs(Matrix *A);
void Mat_Mean(double *mean, Matrix *A);
void Mat_Vector(Matrix *A, float *a);
void Mat_Shift(Matrix *A, Matrix *B, int side);
void Mat_Zeros(Matrix *A);
void Mat_Zeros_uc(uc_Matrix *A);
void Mat_Zeros_i(i_Matrix *A);
void CreateMatrix(Matrix **M, int hei, int wid);
void FreeMatrix(Matrix *M);
void Create_i_Matrix(i_Matrix **M, int hei, int wid);
void Free_i_Matrix(i_Matrix *M);
void Create_uc_Matrix(uc_Matrix **M, int hei, int wid);
void Free_uc_Matrix(uc_Matrix *M);
void Mat_FFT2(Matrix *Output_real, Matrix *Output_imag, Matrix *Input_real, Matrix *Input_imag);
void Mat_IFFT2(Matrix *Output_real, Matrix *Output_imag, Matrix *Input_real, Matrix *Input_imag);
void four2(double **fftr, double **ffti, double **rdata, double **idata, int rs, int cs, int isign);
void four1(double *data, int nn, int isign);
void Mat_Copy(Matrix *A, Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end,
int w_end);
void Mat_uc_Copy(uc_Matrix *A, uc_Matrix *B, int h_target, int w_target, int h_begin, int w_begin,
int h_end, int w_end);
void Mat_i_Copy(i_Matrix *A, i_Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int
h_end, int w_end);
void Mat_Product(Matrix *A, Matrix *B, Matrix *C);
void Mat_Sum(Matrix *A, Matrix *B, Matrix *C);
void Mat_Substract(Matrix *A, Matrix *B, Matrix *C);
void Mat_Fliplr(Matrix *A);
void Mat_Flipud(Matrix *A);
void Mat_uc_Fliplr(uc_Matrix *A);
void Mat_uc_Flipud(uc_Matrix *A);
#endif
#ifdef TRADITIONAL
void nrerror();
double *dvector();
float *vector();
void free_dvector();
void free_vector();
float **matrix();
void free_matrix();
double **dmatrix();
void free_dmatrix();
double logBase2();
void sort();
void minimun();
void Mat_Abs();
void Mat_Mean();
void Mat_Vector();
void Mat_Shift();
void Mat_Zeros();
void Mat_Zeros_uc();
void Mat_Zeros_i();
void CreateMatrix();
void FreeMatrix();
void Create_i_Matrix();
void Free_i_Matrix();
void Create_uc_Matrix();
void Free_uc_Matrix();
void Mat_FFT2();
void Mat_IFFT2();
void four2();
void four1();
void Mat_Copy();
void Mat_uc_Copy();
void Mat_i_Copy();
void Mat_Product();
void Mat_Sum();
void Mat_Substract();
void Mat_Fliplr();
void Mat_Flipud();
void Mat_uc_Fliplr();
void Mat_uc_Flipud();

#endif
