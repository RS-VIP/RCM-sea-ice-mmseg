/** \file
 *
 * \brief Signal processing functions.
 *
 */


#ifndef SIGNALPROCESSING_FUNCS
#define SIGNALPROCESSING_FUNCS

///This function performs fast fourier transform with padding and origin shifted
void fft(float** &imgR, float** &imgI, int &nRows, int &nCols, int &pRow, int &pCol);

///This function performs inverse fast fourier transform with origin shifted back
void ifft(float** imgR, float** imgI, int nRows, int nCols);

///This function pads out the image up to the next power of 2 if necessary and an extra if requested
void padding(float** &imgR, float** &imgI, int &nRows, int &nCols, int &pRow, int &pCol, bool bExtraPadFlag, bool bDcpFftFlag);

/// This function dc pads the image to enable FFT.
/**
 * It is similar in effect to the vip function ``dc_pad(3Vip)''.
 */
void dcp_fft(float **img, int pRow, int pCol, int nNewRows, int nNewCols, int nRows, int nCols );

/// General n-dimensional FFT algorithm.
/**
 * Input is the multi-dimensional real array Ar and imaginary array Ai.  The
 * size of each dimension is given by the array nN.  The number of dimensions
 * is ndim.  The origin of the transformed space is at Ar[0][0]...[0]
*/
void fftn(float *Ar, float *Ai, int *nN, int ndim, int isign );

/// This function shifts the origin of the spectrum from [0][0] to [nrows/2][ncols/2] and normalized elements by nRows*nCols
void fft_shift_org(float **Ar, float **Ai, int nRows, int nCols );

/// This function shifts the origin of the spectrum from [nrows/2][ncols/2] to [0][0]
void ifft_shift_org(float **Ar, float **Ai, int nRows, int nCols );

/// This function  converts images consisting of X and Y components e.g. gradient and complex images to images consisting of radius and angle components (principal values).
void rec2circ(float** &imgR_rad, float** &imgI_ang, int nRows, int nCols);

/// This function  converts images consisting of radius and angle components to images containing X and Y or real and imaginary components.
void circ2rec(float** &rad_imgR, float** &ang_imgI, int nRows, int nCols);

#endif
