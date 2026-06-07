//#pragma warning( disable : 4234 )
#pragma warning( disable : 4244 )


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SignalProcessing.h"
//#include "../header/matrix_gabor.h"

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

#define PI 3.141592653589793115997963468544185161590576171875

/*******************************************************************/
/*	FUNCTION fft				                                   */
/*	This function performs fast fourier transform                  */
/*  with padding and origin shifted                                */
/*******************************************************************/
void fft(float** &imgR, float** &imgI, int &nRows, int &nCols, int &pRow, int &pCol)
{
	int	ndim[2];	/* size of image array [nCols, nRows]*/

	/* Pad up to nearest power of 2 (if necessary) and +1 if requested */
	padding( imgR, imgI, nRows, nCols, pRow, pCol, false, false);

	ndim[0] = nRows;
	ndim[1] = nCols;

	fftn( *imgR, *imgI, ndim, 2, -1 );

		/* Shift origin form [0][0] to [nRows/2][nCols/2] and NORMALIZED by nRows*nCols*/
	fft_shift_org(imgR, imgI, nRows, nCols );
}


/*******************************************************************/
/*	FUNCTION ifft				                                   */
/*	This function performs inverse fast fourier transform	       */
/*  with origin shifted back                                       */
/*******************************************************************/
void ifft(float** imgR, float** imgI, int nRows, int nCols)
{
	int	ndim[2];	/* size of image array [nCols, nRows]*/

	ndim[0] = nRows;
	ndim[1] = nCols;

	/* Shift origin form [nRows/2][nCols/2] to [0][0] */
	ifft_shift_org(imgR, imgI, nRows, nCols );

	fftn( *imgR, *imgI, ndim, 2, 1 );
}

/*******************************************************************/
/*	FUNCTION padding					                           */
/*	This function pads out the image up to the next power of	   */
/*  2 if necessary and an extra power of 2 if requested.           */
/*  This prevents circular convolution effects                     */
/*******************************************************************/
void padding(float** &imgR, float** &imgI, int &nRows, int &nCols, int &pRow, int &pCol, bool bExtraPadFlag, bool bDcpFftFlag)
{
	int	row, col;
	int	nNewRows, nNewCols;
	float **newbuf;
	int	i;

	/* find power of 2 */
	i = 0;
	while ( (1 << i) < nRows )
		i++;
	nNewRows = (1 << i);
	if ( bExtraPadFlag )
		nNewRows *= 2;

	i = 0;
	while ( (1 << i) < nCols )
		i++;
	nNewCols =	1 << i;
	if ( bExtraPadFlag )
		nNewCols *= 2;

	pRow = ( nNewRows - nRows ) / 2;
	pCol = ( nNewCols - nCols ) / 2;

	if ( nNewCols == nCols && nNewRows	== nRows )
		return;

	newbuf = (float **) malloc( nNewRows*sizeof( float* ) );
    if ( newbuf == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    newbuf[0] = (float *) calloc( nNewRows*nNewCols, sizeof( float ) );
    if ( newbuf[0] == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    for ( i = 1; i < nNewRows; i++ )
		newbuf[i] = newbuf[i-1] + nNewCols;

	for ( row = 0; row < nRows; row++ )
		for ( col = 0; col < nCols; col++ ) 
			newbuf[row + pRow][col + pCol] = imgR[row][col];

	if ( bDcpFftFlag ) 
		dcp_fft( newbuf, pRow, pCol, nNewRows, nNewCols, nRows, nCols );
	
	free(imgR[0]);
	free(imgR);
	imgR = newbuf;

	newbuf = (float **) malloc( nNewRows*sizeof( float* ) );
    if ( newbuf == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    newbuf[0] = (float *) calloc( nNewRows*nNewCols, sizeof( float ) );
    if ( newbuf[0] == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    for ( i = 1; i < nNewRows; i++ )
		newbuf[i] = newbuf[i-1] + nNewCols;

	for ( row = 0; row < nRows; row++ )
		for ( col = 0; col < nCols; col++ ) 
			newbuf[row + pRow][col + pCol] = imgI[row][col];

	if ( bDcpFftFlag ) 
		dcp_fft( newbuf, pRow, pCol, nNewRows, nNewCols, nRows, nCols );

	free(imgI[0]);
	free(imgI);
	imgI = newbuf;

	nRows = nNewRows;
	nCols = nNewCols;
}


/*******************************************************************/
/*	Function dcp_fft						                       */
/*	This function dc pads the image to enable FFT.                 */
/*  It is similar in effect to the vip function ``dc_pad(3Vip)''.  */
/*******************************************************************/
void dcp_fft(float **img, int pRow, int pCol, int nNewRows, int nNewCols, int nRows, int nCols )
{
	int	row, col;
	float val1, val2, val3, val4;

	for ( row = pRow; row < pRow + nRows; row++ ) 
	{
		val1 = img[row][pCol];
		val2 = img[row][pCol + nCols -1];
		for ( col = 0; col < pCol; col++ )
		{
			img[row][col] = val1;
			img[row][nNewCols-col-1] = val2;
		}
	}
	for ( col = pCol; col < pCol + nCols; col++ )
	{
		val1 = img[pRow][col];
		val2 = img[pRow + nRows - 1][col];
		for ( row = 0; row < pRow; row++ )
		{
			img[row][col] = val1;
			img[nNewRows - row - 1][col] = val2;
		}
	}
	val1 = img[pRow][pCol];
	val2 = img[pRow][pCol + nCols -1];
	val3 = img[pRow + nRows - 1][pCol];
	val4 = img[pRow + nRows - 1][pCol + nCols - 1];
	for ( row = 0; row <pRow; row++ )
		for ( col = 0; col < pCol; col++ ) 
		{
			img[row][col] = val1;
			img[row][nNewCols - col - 1] = val2;
			img[nNewRows - row - 1][col] = val3;
			img[nNewRows - row - 1][nNewCols - col - 1] = val4;
		}
}


/*******************************************************************/
/*	FUNCTION	FFTN						                       */
/*	General n-dimensional FFT algorithm.  Input is the multi-	   */
/*	dimensional real array Ar and imaginary array Ai.  The size	   */
/*	of each dimension is given by the array nN.  The number of 	   */
/*	dimensions is ndim.						                       */
/*	The origin of the transformed space is at Ar[0][0]...[0]	   */
/*******************************************************************/
void fftn(float *Ar, float *Ai, int *nN, int ndim, int isign )
// Param in: Ar, Ai 1D version of the multidimensional data
//           nN array of dimensions
//           ndim number of dimensions 
//           isign -1 forward FFT, 1 inverse FFT
{
int idim;
unsigned i, jrev, lrev, ibit;
unsigned ip2, ifp1, ifp2, k, N;
unsigned Nprev=1, Nrem, Ntot=1;
register unsigned j, l;
double theta;
double wr,wi;
double wpr,wpi;
float  tempr,tempi;
float  wtr,wti;
float wtemp;
double t1, t2;

for ( idim=0; idim<ndim; ++idim )
	Ntot*=nN[idim];

for ( idim=ndim-1; idim>=0; --idim ) {
	N=nN[idim];
	Nrem=Ntot/(N*Nprev);
	ip2=Nprev*N;        
	jrev=0;
	for ( j=0; j<ip2; j+=Nprev ) {
		if ( j<jrev )
			for ( i=j; i<j+Nprev; ++i )
				for ( l=i; l<Ntot; l+=ip2 ) {
					lrev=l+jrev-j;
					tempr=Ar[l];
					tempi=Ai[l];
					Ar[l]=Ar[lrev];
					Ai[l]=Ai[lrev];
					Ar[lrev]=tempr;
					Ai[lrev]=tempi;
					}
		ibit=ip2;
		do {
			ibit>>=1;
			jrev^=ibit;
			} while ( ibit>=Nprev && ! ( ibit & jrev )  ) ;
		}
	for ( ifp1=Nprev; ifp1<ip2; ifp1 <<=1 )  {
		ifp2=ifp1 << 1;
		theta=isign*2.0*PI/ ( ifp2/Nprev ) ;
		wpr=sin( 0.5*theta ) ;
		wpr*=(-2.0*wpr);
		wpi=sin( theta ) ;
		wr=1.0;
		wi=0.0;
		for ( l=0; l<ifp1; l+=Nprev )  {
			for ( i=l; i<l+Nprev; ++i ) 
				for ( j=i; j<Ntot; j+=ifp2 )  {
					k=j+ifp1;
					wtr=Ar[k];
					wti=Ai[k];
					Ar[k]=Ar[j]-( tempr =
					     ( t1=wr*wtr )-( t2=wi*wti )  );
					Ai[k]=Ai[j]-( tempi =
					     ( wr+wi )*( wtr+wti )-t1-t2 ) ;
					Ar[j]+=tempr;
					Ai[j]+=tempi;
					}
			wtemp=wr;
			wr+=( t1=wr*wpr )-( t2=wi*wpi ) ;
			wi+=( wtemp+wi )*( wpr+wpi )-t1-t2;
			}
		}
	Nprev*=N;
	}
}

/**************************************************************************************/
/*	Function fft_shift_org					                                          */
/*	This function shifts the origin of the spectrum from [0][0] to [nrows/2][ncols/2] */
/*  and normalized elements by nRows*nCols                                            */
/**************************************************************************************/
void fft_shift_org(float **Ar, float **Ai, int nRows, int nCols )
{
	int	row, col;
	float	temp;
	float	n;

	n = nRows*nCols;
	for ( row = 0; row < nRows; row ++)
		for ( col = 0; col < nCols/2; col++ ) 
		{
			temp = Ar[row][col]/n;
			Ar[row][col] = Ar[row][col + nCols/2]/n;
			Ar[row][col + nCols/2] = temp;
			temp = Ai[row][col]/n;
			Ai[row][col] = Ai[row][col + nCols/2]/n;
			Ai[row][col + nCols/2] = temp;
		}
	for ( col = 0; col < nCols; col++ )
		for ( row = 0; row < nRows/2; row++ )
		{
			temp = Ar[row][col];
			Ar[row][col] = Ar[row + nRows/2][col];
			Ar[row + nRows/2][col] = temp;
			temp = Ai[row][col];
			Ai[row][col] = Ai[row + nRows/2][col];
			Ai[row + nRows/2][col] = temp;
		}
}

/**************************************************************************************/
/*	Function ifft_shift_org					                                          */
/*	This function shifts the origin of the spectrum from [nrows/2][ncols/2] to [0][0] */
/**************************************************************************************/
void ifft_shift_org(float **Ar, float **Ai, int nRows, int nCols )
{
	int	row, col;
	float	temp;

	for ( row = 0; row < nRows; ++row )
		for ( col=0; col<nCols/2; ++col ) 
		{
			temp=Ar[row][col];
			Ar[row][col]=Ar[row][col+nCols/2];
			Ar[row][col+nCols/2]=temp;
			temp=Ai[row][col];
			Ai[row][col]=Ai[row][col+nCols/2];
			Ai[row][col+nCols/2]=temp;
		}

	for ( col=0; col<nCols; ++col )
		for ( row=0; row<nRows/2; ++row ) 
		{
			temp=Ar[row][col];
			Ar[row][col]=Ar[row+nRows/2][col];
			Ar[row+nRows/2][col]=temp;
			temp=Ai[row][col];
			Ai[row][col]=Ai[row+nRows/2][col];
			Ai[row+nRows/2][col]=temp;
		}
}

/**************************************************************************************/
/*	Function rec2circ				                                                  */
/*	This function  converts images consisting of X and Y components                   */
/*  e.g. gradient and complex images to images consisting of radiu                    */
/*  and angle components (principal values).			                     		  */
/**************************************************************************************/
void rec2circ(float** &imgR_rad, float** &imgI_ang, int nRows, int nCols)
{
	int	row, col;	/* image row and column indices */
	int i;
	float **bufRad, **bufAng;

	bufRad = (float **) malloc( nRows*sizeof( float* ) );
    if ( bufRad == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    bufRad[0] = (float *) calloc( nRows*nCols, sizeof( float ) );
    if ( bufRad[0] == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}

	bufAng = (float **) malloc( nRows*sizeof( float* ) );
    if ( bufAng == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    bufAng[0] = (float *) calloc( nRows*nCols, sizeof( float ) );
    if ( bufAng[0] == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}

    for ( i = 1; i < nRows; i++ )
	{
		bufRad[i] = bufRad[i-1] + nCols;
		bufAng[i] = bufAng[i-1] + nCols;
	}

	for ( row = 0; row < nRows; row++ )
		for ( col = 0; col < nCols; col++ )
			bufRad[row][col] = sqrt( imgR_rad[row][col]*imgR_rad[row][col] + imgI_ang[row][col]*imgI_ang[row][col] );

	for ( row = 0; row < nRows; row++ )
		for ( col = 0; col < nCols; col++ )
		{
			//Edited by Michael Stone - May 2016 - changed abs() to fabs() to return float absolutes and prevent compilation errors
			//(abs() is an overloaded function and passing in a float is an arbitrary call as all three instances of the function
			//can accept it, none explicitly take a float as the parameter)
			if ( ( fabs( (imgR_rad[row][col]) ) > ZERO_THRESH ) || ( fabs (imgI_ang[row][col] ) > ZERO_THRESH ) )
				bufAng[row][col] = atan2( imgI_ang[row][col], imgR_rad[row][col] );
			else
				bufAng[row][col] = 0.0;
		}

	free(imgR_rad[0]);
	free(imgR_rad);
	imgR_rad = bufRad;

	free(imgI_ang[0]);
	free(imgI_ang);
	imgI_ang = bufAng;
}

/**************************************************************************************/
/*	Function circ2rec				                                                  */
/*	This function  converts images consisting of radius and angle components to	      */
/*	images containing X and Y or real and imaginary components.						  */
/**************************************************************************************/
void circ2rec(float** &rad_imgR, float** &ang_imgI, int nRows, int nCols)
{
	int	row, col;	/* image row and column indices */
	int i;
	float **bufR, **bufI;

	bufR = (float **) malloc( nRows*sizeof( float* ) );
    if ( bufR == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    bufR[0] = (float *) calloc( nRows*nCols, sizeof( float ) );
    if ( bufR[0] == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}

	bufI = (float **) malloc( nRows*sizeof( float* ) );
    if ( bufI == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}
    bufI[0] = (float *) calloc( nRows*nCols, sizeof( float ) );
    if ( bufI[0] == NULL ) 
	{
		fprintf( stderr, "Fatal Error - out of memory.\n");
        exit(-1);
	}

    for ( i = 1; i < nRows; i++ )
	{
		bufR[i] = bufR[i-1] + nCols;
		bufI[i] = bufI[i-1] + nCols;
	}

	for ( row = 0; row < nRows; row++ )
		for ( col = 0; col < nCols; col++ )
			bufR[row][col] = cos( ang_imgI[row][col] )*rad_imgR[row][col];

	for ( row = 0; row < nRows; row++ )
		for ( col = 0; col < nCols; col++ )
			bufI[row][col] = sin( ang_imgI[row][col] )*rad_imgR[row][col];

	free(rad_imgR[0]);
	free(rad_imgR);
	rad_imgR = bufR;

	free(ang_imgI[0]);
	free(ang_imgI);
	ang_imgI = bufI;
}