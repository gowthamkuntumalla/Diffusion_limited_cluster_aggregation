

//////#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>   


/* number of particles */
#define n 10000000

#define pi 4.0*atan(1.0)

#define MASS 100000000

#define dMASS 100000000

// for projected are are monomers spheres (mode 0) or cubes (mode 1)
#define mode 0


/////////////////Rotaion //////////////////////////////////////////////////////////////////////////////////////////////////////////

double Inner_Product(double u[], double c[], int x);

void Copy_Matrix(double *A, double *B, int nrows, int ncols);

void Transform_Matrix(double* U, double s, double c, int col, int y);

int Tridiagonalize_Symmetric_Matrix(double *A, int storage_class, double diagonal[], double off_diagonal[], double *U, int x);

int QL_Tridiagonal_Symmetric_Matrix(double diagonal[], double off_diagonal[], double *U, int x, int max_iteration_count);

void aling_major(int lead, double x[], double y[], double z[], double a_x, double a_y, double a_z, int clustnext[], int clustlead[], int clustnum[]);

int aling_minor(int lead, double x[], double y[], double z[], double a_x, double a_y, double a_z, int clustnext[], int clustlead[], int clustnum[]);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double get_length(int lead, double x[], double y[], double z[], double X_length[], double Y_length[], double Z_length[], int clustnext[], int clustlead[], int clustnum[], double rmsX[], double rmsY[], double rmsZ[]);

double get_area(int i, double x[], double y[], double z[], double X_length[], double Y_length[], double Z_length[], int clustnext[], int clustlead[], int clustnum[], double RG_2d[]);



int main()
{

	double *X, *Y, *Z, *X_length, *Y_length, *Z_length, *Rg, junk, fv, *rmsX, *rmsY, *rmsZ, *R_proj, Rms_ep, x, *RG_2d;
	int *clustnext, *clustlead, *clustnum, *mass;
	int i, Nc, lead;

	///////////////////////////////////////////////////////////////////
	                                                                 //
	X = (double *)malloc(n*sizeof(double));                          //                         
	if (X == NULL)                                                   //                
	{                                                                //
		printf("X");                                                 //                  
	}                                                                // 																
	                                                                 //
	Y = (double *)malloc(n*sizeof(double));                          //                         
	if (Y == NULL)                                                   //                
	{                                                                //
		printf("Y");                                                 //                  
	}                                                                // 																
	                                                                 //
	Z = (double *)malloc(n*sizeof(double));                          //                         
	if (Z == NULL)                                                   //                
	{                                                                //
		printf("Z");                                                 //                  
	}                                                                // 																
	                                                                 //
	clustnum = (int*)malloc(n*sizeof(int));                          //
	if (clustnum == NULL)                                            //
	{                                                                //
		printf("clustnum\n");                                        //
		return 0;                                                    //
	}                                                                //
	                                                                 //
	clustlead = (int*)malloc(n*sizeof(int));                         //
	if (clustlead == NULL)                                           //
	{                                                                //
		printf("clustlead\n");                                       //
		return 0;                                                    //
	}                                                                //
	                                                                 //
	clustnext = (int*)malloc(n*sizeof(int));                         //
	if (clustnext == NULL)                                           //
	{                                                                //
		printf("clustnext\n");                                       //
		return 0;                                                    //
	}															     //
	                                                                 //
	///////////////////////////////////////////////////////////////////



	FILE*unfold;
	FILE*vol_frac;

	unfold = fopen("input", "r");

	for (i = 0; i<n; i++)
	{
		/*Too many printf slow down the performance*/
		//printf("%d\n", i);
		//printf("%d\n", i);
		// fscanf(unfold,"%lf %lf %lf %lf %lf %d %d ",&x, &X[i], &Y[i], &Z[i],&x, &clustlead[i], &clustnext[i]);  
		fscanf(unfold, "%lf %lf %lf %d %d ", &X[i], &Y[i], &Z[i], &clustlead[i], &clustnext[i]);

		// single cluster //////////////
		//fscanf(unfold,"%lf %lf %lf",&X[i], &Y[i], &Z[i]  );   clustlead[i]=0; clustnext[i]=i+1; 

		//if( i==(n-1) )
		//{
			//clustnext[i]=-1;
		//}
		///////////////////////////////
	}
	fclose(unfold);



	//set all of clustnum to -1
	for (i = 0; i<n; i++)
	{
		clustnum[i] = -1;
	}
	Nc = 0;
	//assings each clust a number
	for (i = 0; i<n; i++)
	{
		if (i == clustlead[i])
		{
			clustnum[Nc] = i;
			Nc++;
		}
	}



	////////////////////////////////////////////////////////////////////
	X_length = (double *)malloc(Nc*sizeof(double));                   //                         
	if (X_length == NULL)                                             //                
	{                                                                 //
		printf("X_length");                                           //                  
	}                                                                 //
	                                                                  //
	                                                                  //
	Y_length = (double *)malloc(Nc*sizeof(double));                   //                         
	if (Y_length == NULL)                                             //                
	{                                                                 //
		printf("Y_length");                                           //                  
	}                                                                 //
	                                                                  //
	                                                                  //
	                                                                  //
	Z_length = (double *)malloc(Nc*sizeof(double));                   //                         
	if (Z_length == NULL)                                             //                
	{                                                                 //
		printf("Z_length");                                           //                  
	}                                                                 //
	                                                                  //
	Rg = (double *)malloc(Nc*sizeof(double));                         //                         
	if (Rg == NULL)                                                   //                
	{                                                                 //
		printf("Rg");                                                 //                  
	}                                                                 //
	                                                                  //
	R_proj = (double *)malloc(Nc*sizeof(double));                     //                         
	if (R_proj == NULL)                                               //                
	{                                                                 //
		printf("R_proj");                                             //                  
	}                                                                 //
	                                                                  //
	mass = (int *)malloc(Nc*sizeof(int));                             //                         
	if (mass == NULL)                                                 //                
	{                                                                 //
		printf("mass");                                               //                  
	}                                                                 //
	                                                                  //
	rmsX = (double *)malloc(Nc*sizeof(double));                       //                         
	if (rmsX == NULL)                                                 //                
	{                                                                 //
		printf("rmsX");                                               //                  
	}                                                                 // 																
	                                                                  //
	rmsY = (double *)malloc(Nc*sizeof(double));                       //                         
	if (rmsY == NULL)                                                 //                
	{                                                                 //
		printf("rmsY");                                               //                  
	}                                                                 // 																
	                                                                  //
	rmsZ = (double *)malloc(Nc*sizeof(double));                       //                         
	if (rmsZ == NULL)                                                 //                
	{                                                                 //
		printf("rmsZ");                                               //                  
	}                                                                 //
	RG_2d = (double *)malloc(Nc*sizeof(double));                      //                         
	if (R_proj == NULL)                                               //                
	{                                                                 //
		printf("RG_2d");                                              //                  
	}                                                                 //		
	                                                                  //
	////////////////////////////////////////////////////////////////////


	for (i = 0; i<Nc; i++)
	{
		

		lead = clustnum[i];

		aling_major(lead, X, Y, Z, 1, 0, 0, clustnext, clustlead, clustnum);

		mass[i] = aling_minor(lead, X, Y, Z, 0, 1, 0, clustnext, clustlead, clustnum);

		Rg[i] = get_length(i, X, Y, Z, X_length, Y_length, Z_length, clustnext, clustlead, clustnum, rmsX, rmsY, rmsZ);

		R_proj[i] = get_area(i, X, Y, Z, X_length, Y_length, Z_length, clustnext, clustlead, clustnum, RG_2d);

		printf("%d / %d \n", i, Nc);

	}


	FILE *output;
	output = fopen("Lenghts", "w");

	vol_frac = fopen("vol_frac.txt", "w");

	fprintf(output, "i, mass[i], X_length, Z_length, Y_length, Rg, R_proj, RG_2d \n");


	for (i = 0; i<Nc; i++)
	{

		if ((mass[i] <= (MASS + dMASS)) && (mass[i] >= (MASS - dMASS)))
		{
			if (X_length[i] <0.0005 || Y_length[i] <0.005 || Z_length[i] <0.0005)
			{
			}
			else
			{
				Rms_ep = rmsX[i] * rmsY[i] * rmsZ[i]; Rms_ep = pow(Rms_ep, (1.0 / 3.0));
				fv = (mass[i] * 0.5*0.5*0.5) / (rmsX[i] * rmsY[i] * rmsZ[i]);

				fprintf(vol_frac, "%d %lf %lf %lf %lf \n", i, Rms_ep, rmsX[i], rmsY[i], rmsZ[i]);
				fprintf(output, "%d %d %lf %lf %lf %lf %lf %lf \n", i, mass[i], X_length[i], Z_length[i], Y_length[i], Rg[i], R_proj[i], RG_2d[i]);
			}
		}
	}

	fclose(output);
	fclose(vol_frac);

	return 0;
}

////////////// Find egen vals and vec /////////////////////////////////////
double Inner_Product(double u[], double c[], int x)
{
	int i;
	double ans;
	ans = 0;
	for (i = 0; i<x; i++)
	{
		ans = ans + (u[i] * c[i]);
	}

	return ans;
}
void Copy_Matrix(double *A, double *B, int nrows, int ncols)
{
	int i, j;

	for (i = 0; i<nrows; i++)
	{
		for (j = 0; j<ncols; j++)
		{
			*A = *B;
			*B++; *A++;
		}
	}

}
/*
////////////////////////////////////////////////////////////////////////////////
//  int Tridiagonalize_Symmetric_Matrix(double *A, int storage_class,         //
//               double diagonal[], double off_diagonal[], double *U, int n)  //
//                                                                            //
//  Description:                                                                //
//     This program transforms the symmetric square matrix A to a similar     //
//     tridiagonal matrix by a multiplying A on the right and left by a       //
//     sequence of Householder transformations.                               //
//     Def:  Two matrices A and B are said to be similar if there exists a    //
//           nonsingular matrix S such that A S = S B.                        //
//     Def:  A Householder transformation is an orthogonal transformation of  //
//           the form Q = I - 2 uu'/u'u, where u is a n x 1 column matrix and //
//           ' denotes the transpose.                                         //
//     Thm:  If Q is a Householder transformation then Q' = Q  and  Q' Q = I, //
//           i.e. Q is a symmetric involution.                                //
//     First the input matrix A is copied / expanded to form the matrix U.    //
//     The algorithm proceeds by successivly selecting rows i = n - 1, ..., 1 //
//     and then calculating the Householder transformation Q which annihilates//
//     the components to the left of the subdiagonal for that row and leaves  //
//     the previously selected rows invariant.  The algorithm then updates    //
//     the matrix U, in place, by premultiplication by Q followed by          //
//     postmultiplication by Q.                                               //
//     If the i-th row of U is (a[0],...,a[n-1]), then  choose u' =           //
//     (u[0],...,u[n-1]) by u[0] = 0, ... , u[j] = 0, u[j+2] = a[j+2],...,    //
//     u[n-1] = a[n-1].  The remaining component u[j+1] = a[j+1] - s, where   //
//     s^2 = a[j+1]^2 + ... + a[n-1]^2, and the choice of sign for s,         //
//     sign(s) = -sign(a[j+1]) maximizes the number of significant bits for   //
//     u[j+1].                                                                //
//                                                                            //
//     Note: In general accuracy is improved rows and columns of the matrix   //
//     are ordered so that the elements with the smallest magnitude are       //
//     located at the top left-most corner of the matrix and the elements with//
//     the largest magnitude are located at the bottom right-most corner.     //
//                                                                            //
//  Arguments:                                                                //
//     double *A                                                              //
//            Pointer to the first element of the matrix A[n][n].             //
//     int    storage_class                                                   //
//            The storage class of the matrix A, 0 if full symmetric matrix   //
//            -1 if stored in lower triangular form, +1 if stored in upper    //
//            triangular form.                                                //
//     double diagonal[]                                                      //
//            The diagonal of the tridiagonal symmetric matrix.  The array    //
//            diagonal[] should be dimensioned n in the calling routine.      //
//     double off_diagonal[]                                                  //
//            The subdiagonal and superdiagonal of the tridiagonal symmetric  //
//            matrix.  The subdiagonal and superdiagonal are (off_diagonal[1],//
//            ,,, off_diagonal[n-1]), off_diagonal[0] is set to 0.0.  The     //
//            array off_diagonal[] should be dimensioned n in the calling     //
//            routine.  If QL_Tridiagonal_Symmetric_Matrix is called using    //
//            using diagonal[], off_diagonal[], and U, then off_diagonal      //
//            should be dimensioned n + 1 in the calling routine.             //
//     double *U                                                              //
//            The n x n orthogonal matrix such that AU = UT, where T is the   //
//            tridiagonal matrix with diagonal 'diagonal[]' and off-diagonal  //
//            'off_diagonal[]'.                                               //
//     int    n                                                               //
//            The number of rows or columns of the symmetric matrix A and of  //
//            the matrix U.                                                   //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - Illegal storage class.                                    //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N * (N + 1) / 2], U[N][N], diagonal[N];              //
//     double off_diagonal[N];                                                //
//                                                                            //
//     (your code to create the matrices A and B)                             //
//     Tridiagonalize_Symmetric_Matrix((double*)A, 0, diagonal, off_diagonal, //
//                                                          (double*)U, n);   //
//     (if B is the matrix stored in lower triangular form )                  //
//     Tridiagonalize_Symmetric_Matrix((double*)B, -1, diagonal, off_diagonal,//
//                                                          (double*)U, n);   //
//     (if B is the matrix stored in upper triangular form )                  //
//     Tridiagonalize_Symmetric_Matrix((double*)B, 1, diagonal, off_diagonal, //
//                                                          (double*)U, n);   //
////////////////////////////////////////////////////////////////////////////////
*/
//
void Transform_Matrix(double* U, double s, double c, int col, int y)
{
	double x;
	int i;
	int col1 = col + 1;

	for (i = 0; i < y; i++, U += y) {
		x = *(U + col1);
		*(U + col1) = s * *(U + col) + c * x;
		*(U + col) *= c;
		*(U + col) -= s * x;
	}
}
static double Calculate_Shift(double d2, double d1, double off)
{
	double h;

	h = (d2 - d1) / (off + off);
	off = sqrt(h * h + 1.0);
	d2 = (h < 0.0) ? h - off : h + off;
	return d1 - off / d2;
}
int Tridiagonalize_Symmetric_Matrix(double *A, int storage_class, double diagonal[], double off_diagonal[], double *U, int x)

{
	int i, j, k, col;
	double zero_tolerance = DBL_MIN / DBL_EPSILON;
	double *pU, *ppU, *pUrow, *pUcol;
	double ss;                             // signed sqrt of sum of squares

	// n x n matrices for which n <= 2 are already in Tridiagonal form

	if (x <= 2) {
		for (i = 0, k = 0; i < x; i++) {
			for (j = 0; j < i; j++, k++) *(U + k) = 0.0;
			*(U + k++) = 1.0;
			for (j = i + 1; j < x; j++, k++) *(U + k) = 0.0;
		}
		return 0;
	}

	// Initialize U to A.

	switch (storage_class) {
	case -1:  // A stored in lower triangular form
		pUrow = U;
		pUcol = U;
		for (i = 0; i < x; pUcol++, pUrow += x)
		for (j = 0, pU = pUcol; j <= i; j++, pU += x) {
			*(pUrow + j) = *A;
			*pU = *A++;
		}
		break;
	case 0:   // A stored in full matrix form
		Copy_Matrix(U, A, x, x);
		break;
	case 1:   // A stored in upper triangular form
		for (i = 0, pUcol = U; i < x; pUcol += x + 1) {
			pUrow = pUcol;
			for (j = i, pU = pUcol; j < x; j++, pU += x) {
				*pUrow++ = *A;
				*pU = *A++;
			}
		}
		break;
	default:  // Illegal storage class
		return -1;
	}

	// For each row starting with the last row find the 
	// Householder transformation to zero all entries to
	// the left of the subdiagonal.

	for (i = x - 1, pUrow = U + x * (x - 1); i > 0; i--, pUrow -= x) {

		// Calculate the sum of squares of the elements to the left
		// of the subdiagonal.

		pU = pUrow + i - 1;
		ss = Inner_Product(pUrow, pUrow, i - 1);

		// If the sum of squares is too small set, don't transform.

		if (ss < zero_tolerance) {
			off_diagonal[i] = *pU;
			diagonal[i] = 0.0;
			continue;
		}

		// Calculate the signed square root of the sum of squares to the
		// left of the diagonal.

		ss += *pU * *pU;
		off_diagonal[i] = (*pU >= 0.0) ? -sqrt(ss) : sqrt(ss);
		diagonal[i] = ss - *pU * off_diagonal[i];
		*pU -= off_diagonal[i];
		for (j = 0, pU = U, ss = 0.0; j < i; j++, pU += x) {
			*(pU + i) = *(pUrow + j) / diagonal[i];
			off_diagonal[j] = Inner_Product(pU, pUrow, j + 1);
			for (k = j + 1, ppU = pU + x; k < i; k++, ppU += x)
				off_diagonal[j] += *(ppU + j) * *(pUrow + k);
			ss += off_diagonal[j] * *(pU + i);
			off_diagonal[j] /= diagonal[i];
		}
		ss /= (diagonal[i] + diagonal[i]);
		for (j = 0, pU = U; j < i; j++, pU += x) {
			off_diagonal[j] -= ss * *(pUrow + j);
			for (k = 0; k <= j; k++)
				*(pU + k) -= (*(pUrow + j) * off_diagonal[k]
				+ off_diagonal[j] * *(pUrow + k));
		}
	}
	diagonal[0] = 0.0;
	off_diagonal[0] = 0.0;
	for (i = 0, pUrow = U; i < x; i++, pUrow += x) {
		if (diagonal[i] != 0.0)
		for (j = 0; j < i; j++) {
			for (k = 0, pU = U, ss = 0.0; k < i; k++, pU += x)
				ss += *(pUrow + k) * *(pU + j);
			for (k = 0, pU = U; k < i; k++, pU += x)
				*(pU + j) -= ss * *(pU + i);
		}
		diagonal[i] = *(pUrow + i);
		*(pUrow + i) = 1.0;
		for (j = 0, pU = U; j < i; j++, pU += x) {
			*(pUrow + j) = 0.0;
			*(pU + i) = 0.0;
		}
	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
//  int QL_Tridiagonal_Symmetric_Matrix(double diagonal[],                    //
//          double off_diagonal[], double *U, int n, int max_iteration_count) //
//                                                                            //
//  Description:                                                              //
//     This function implements the QL algorithm with implicit shifts of the  //
//     of the origin for a symmetric tridiagonal matrix.  If the original     //
//     symmetric matrix was tridiagonalized using Householder transformations //
//     then the tridiagonal matrix should be organized so that the largest    //
//     elements are at the bottom right-hand of the matrix.  Also, the process//
//     of tridiagonalizing the matrix will introduce round-off errors, so that//
//     the estimates of the eigenvalues and eigenvectors may be improved by   //
//     calling the inverse_power_method.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double diagonal[]                                                      //
//            On input, the diagonal of the tridiagonal symmetric matrix.     //
//            On output, the eigenvalues.                                     //
//     double off_diagonal[]                                                  //
//            The subdiagonal and superdiagonal of the tridiagonal symmetric  //
//            matrix.  The subdiagonal and superdiagonal are (off_diagonal[1],//
//            ,,, off_diagonal[n-1]).                                         //
//     double *U                                                              //
//            On input, U[n][n] is the identity matrix if the tridiagonal     //
//            matrix is the primary data or is the transformation matrix      //
//            of Householder transformations which is the output of           //
//            tridiagonalize_symmetric_matrix.  On output, U[n][n] is the     //
//            matrix of eigenvectors, the ith column being the eigenvector    //
//            corresponding to the eigenvalue d[i].  These are the eigen-     //
//            vectors of the full matrix, if on input U is not the identity.  //
//     int    n                                                               //
//            The number of rows or columns of the symmetric matrix A and of  //
//            the matrix U.                                                   //
//     int    max_iteration_count                                             //
//            The maximum number of iterations to try to annihilate an off-   //
//            diagonal elements.                                              //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - Failed to converge within max_iteration_count iterations. //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double U[N][N], diagonal[N], double off_diagonal[N];                   //
//     int err;                                                               //
//     int max_iteration_count = 30;                                          //
//                                                                            //
//     (your code to create the vectors diagonal and off_diagonal)            //
//     if (QL_Tridiagonal_Symmetric_Matrix(diagonal, off_diagonal, (double*)U,//
//         N, max_iteration_count) ) printf("Failed to converge\n");          //
//     else printf("Success\n");                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int QL_Tridiagonal_Symmetric_Matrix(double diagonal[], double off_diagonal[], double *U, int x, int max_iteration_count)

{
	int i, j, k;
	int iteration;
	double epsilon;
	double *p_off = off_diagonal + 1;
	double s, c, g, h, q;
	double shift;
	double dum;


	p_off[x - 1] = 0.0;
	for (i = 0; i < x; i++) {
		for (iteration = 0; iteration < max_iteration_count; iteration++) {
			for (k = i; k < (x - 1); k++) {
				epsilon = DBL_EPSILON * (fabs(diagonal[k]) + fabs(diagonal[k + 1]));
				if (fabs(p_off[k]) <= epsilon) break;
			}
			if (k == i) break;
			shift = Calculate_Shift(diagonal[i + 1], diagonal[i], p_off[i]);
			q = diagonal[k] - shift;
			c = 1.0;
			s = 1.0;
			shift = 0.0;
			for (j = k - 1; j >= i; j--) {
				h = c * p_off[j];
				g = s * p_off[j];
				if (fabs(g) >= fabs(q)) {
					c = q / g;
					dum = sqrt(c * c + 1.0);
					p_off[j + 1] = g * dum;
					s = 1.0 / dum;
					c *= s;
				}
				else {
					s = g / q;
					dum = sqrt(s * s + 1.0);
					p_off[j + 1] = q * dum;
					c = 1.0 / dum;
					s *= c;
				}
				q = diagonal[j + 1] - shift;
				dum = s * (diagonal[j] - q) + 2.0 * c * h;
				shift = s * dum;
				diagonal[j + 1] = q + shift;
				q = c * dum - h;
				Transform_Matrix(U, s, c, j, x);
			}
			diagonal[i] -= shift;
			p_off[i] = q;
			p_off[k] = 0.0;
		}
	}
	if (iteration >= max_iteration_count) return -1;
	return 0;
}
void aling_major(int lead, double x[], double y[], double z[], double a_x, double a_y, double a_z, int clustnext[], int clustlead[], int clustnum[])
{
	int i, j, mass;
	double Px, Py, Pz;
	double xc, yc, zc;
	double Ixx, Ixy, Ixz, Iyy, Iyz, Izz;
	double A[3][3], U[3][3], R[3][3];
	double diagonal[3], off_diagonal[4];
	double ex, ey, ez, e1, e2, e3;
	double mag1, mag2, mag3;
	double cx, cy, cz, dot;
	double ang, s, c, u, V;
	double multi;



	// make aling vector unit ///////
	mag1 = a_x*a_x + a_y*a_y + a_z*a_z;
	mag1 = sqrt(mag1);

	a_x = a_x / mag1;  a_y = a_y / mag1; a_z = a_z / mag1;
	/////////////////////////////////////


	i = lead;

	//............ find center of mass ...........///

	xc = 0; yc = 0; zc = 0;
	mass = 0;
	while (i != -1)
	{
		xc = xc + x[i];
		yc = yc + y[i];
		zc = zc + z[i];
		mass++;
		i = clustnext[i];

	}
	xc = xc / mass; yc = yc / mass; zc = zc / mass;

	//...... ............................... ......///




	//////// find egen val and vect of inertia tensor///////////////////////////////////////////////


	Ixx = 0; Ixy = 0; Ixz = 0;
	Iyy = 0; Iyz = 0;
	Izz = 0;

	i = lead;
	while (i != -1)
	{
		Ixx = Ixx + (((y[i] - yc)*(y[i] - yc)) + ((z[i] - zc)*(z[i] - zc)));
		Iyy = Iyy + (((x[i] - xc)*(x[i] - xc)) + ((z[i] - zc)*(z[i] - zc)));
		Izz = Izz + (((x[i] - xc)*(x[i] - xc)) + ((y[i] - yc)*(y[i] - yc)));

		Ixy = Ixy + (x[i] - xc)*(y[i] - yc);
		Ixz = Ixz + (x[i] - xc)*(z[i] - zc);
		Iyz = Iyz + (y[i] - yc)*(z[i] - zc);

		i = clustnext[i];
	}
	Ixy = -1 * Ixy; Ixz = -1 * Ixz; Iyz = -1 * Iyz;

	multi = 1.0 / mass;
	A[0][0] = multi*Ixx;   A[0][1] = multi*Ixy;    A[0][2] = multi*Ixz;

	A[1][0] = multi*Ixy;   A[1][1] = multi*Iyy;    A[1][2] = multi*Iyz;

	A[2][0] = multi*Ixz;   A[2][1] = multi*Iyz;    A[2][2] = multi*Izz;


	U[0][0] = U[0][1] = U[0][2] = 0;
	U[1][0] = U[1][1] = U[1][2] = 0;
	U[2][0] = U[2][1] = U[2][2] = 0;

	j = Tridiagonalize_Symmetric_Matrix(*A, 0, diagonal, off_diagonal, *U, 3);
	i = QL_Tridiagonal_Symmetric_Matrix(diagonal, off_diagonal, *U, 3, 3000);
	if (i == -1 || j == -1)
	{
		printf("FAIL \n");
	}


	e1 = diagonal[0]; e2 = diagonal[1]; e3 = diagonal[2];


	if (e1 <= e2 && e1 <= e3)
	{
		ex = U[0][0]; ey = U[1][0]; ez = U[2][0];
		//printf("e1 \n");
	}
	else
	{
		if (e2 <= e3)
		{
			ex = U[0][1]; ey = U[1][1]; ez = U[2][1];
			//printf("e2 \n");
		}
		else
		{
			ex = U[0][2]; ey = U[1][2]; ez = U[2][2];
			//printf("e3 \n");
		}
	}

	mag2 = ex*ex + ey*ey + ez*ez;
	mag2 = sqrt(mag2);

	ex = ex / mag2; ey = ey / mag2; ez = ez / mag2;


	//////// find egen val and vect of inertia tensor...done////////////////////////////////////




	//cross product a X e. cross product is the axis of rotion////
	cx = (a_y*ez) - (a_z*ey);
	cy = (a_z*ex) - (a_x*ez);
	cz = (a_x*ey) - (a_y*ex);

	mag3 = cx*cx + cy*cy + cz*cz;
	mag3 = sqrt(mag3);

	cx = -cx / mag3; cy = -cy / mag3; cz = -cz / mag3;
	if (mag3 == 0)
	{
		cx = 0; cy = 0; cz = 0;
	}
	/////////////////////////////////////////////////////////////


	// dot product a * e  and angle difference //////////
	dot = a_x*ex + a_y*ey + a_z*ez;

	ang = 1 * acos(dot);

	//printf("e_vals (%lf, %lf ,%lf) \n", e1,e2,e3);
	//printf("a X e (%lf, %lf ,%lf) \n", cx, cy, cz);
	//printf("e_vec (%lf, %lf ,%lf) ang: %lf\n", ex,ey,ez,(180*ang/pi));

	//	if(ang> (pi/2))
	//{
	//ang=(pi)-ang;
	//ang=-1.0*ang;
	//}
	///////////////////////////////////////////////////





	// set up rotion matrix //////////////////////////


	s = sin(ang);
	c = cos(ang);
	u = (1 - cos(ang));
	V = cos(ang / 2);

	R[0][0] = c + cx*cx * u;
	R[0][1] = cx*cy * u - cz * s;
	R[0][2] = cx*cz * u + cy * s;

	R[1][0] = cx*cy * u + cz * s;
	R[1][1] = c + cy*cy * u;
	R[1][2] = cy*cz * u - cx * s;

	R[2][0] = cx*cz * u - cy * s;
	R[2][1] = cy*cz * u + cx * s;
	R[2][2] = c + cz*cz * u;
	///////////////////////////////////////////////		

	//printf("U\n");
	//printf("%0.2lf   %0.2lf   %0.2lf \n",U[0][0],U[0][1],U[0][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",U[1][0],U[1][1],U[1][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",U[2][0],U[2][1],U[2][2]);
	//printf("\n");

	//printf("R\n");
	//printf("%0.2lf   %0.2lf   %0.2lf \n",R[0][0],R[0][1],R[0][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",R[1][0],R[1][1],R[1][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",R[2][0],R[2][1],R[2][2]);
	//printf("\n");


	/// time to overwrite x,y,z //////////////////////////////////////////////////////////////




	i = lead;
	while (i != -1)
	{
		Px = (R[0][0] * (x[i] - xc)) + (R[0][1] * (y[i] - yc)) + (R[0][2] * (z[i] - zc));

		Py = (R[1][0] * (x[i] - xc)) + (R[1][1] * (y[i] - yc)) + (R[1][2] * (z[i] - zc));

		Pz = (R[2][0] * (x[i] - xc)) + (R[2][1] * (y[i] - yc)) + (R[2][2] * (z[i] - zc));

		x[i] = Px;

		y[i] = Py;

		z[i] = Pz;

		i = clustnext[i];

	}












}
int aling_minor(int lead, double x[], double y[], double z[], double a_x, double a_y, double a_z, int clustnext[], int clustlead[], int clustnum[])
{
	int i, j, mass;
	double Px, Py, Pz;
	double xc, yc, zc;
	double Ixx, Ixy, Ixz, Iyy, Iyz, Izz;
	double A[3][3], U[3][3], R[3][3];
	double diagonal[3], off_diagonal[4];
	double ex, ey, ez, e1, e2, e3;
	double mag1, mag2, mag3;
	double cx, cy, cz, dot;
	double ang, s, c, u, V;
	double multi;

	char filename[30];

	FILE* output;


	// make aling vector unit ///////
	mag1 = a_x*a_x + a_y*a_y + a_z*a_z;
	mag1 = sqrt(mag1);

	a_x = a_x / mag1;  a_y = a_y / mag1; a_z = a_z / mag1;
	/////////////////////////////////////


	i = lead;

	//............ find center of mass ...........///

	xc = 0; yc = 0; zc = 0;
	mass = 0;
	while (i != -1)
	{
		xc = xc + x[i];
		yc = yc + y[i];
		zc = zc + z[i];
		mass++;
		i = clustnext[i];

	}
	xc = xc / mass; yc = yc / mass; zc = zc / mass;

	//...... ............................... ......///




	//////// find egen val and vect of inertia tensor///////////////////////////////////////////////


	Ixx = 0; Ixy = 0; Ixz = 0;
	Iyy = 0; Iyz = 0;
	Izz = 0;

	i = lead;
	while (i != -1)
	{
		Ixx = Ixx + (((y[i] - yc)*(y[i] - yc)) + ((z[i] - zc)*(z[i] - zc)));
		Iyy = Iyy + (((x[i] - xc)*(x[i] - xc)) + ((z[i] - zc)*(z[i] - zc)));
		Izz = Izz + (((x[i] - xc)*(x[i] - xc)) + ((y[i] - yc)*(y[i] - yc)));

		Ixy = Ixy + (x[i] - xc)*(y[i] - yc);
		Ixz = Ixz + (x[i] - xc)*(z[i] - zc);
		Iyz = Iyz + (y[i] - yc)*(z[i] - zc);

		i = clustnext[i];
	}
	Ixy = -1 * Ixy; Ixz = -1 * Ixz; Iyz = -1 * Iyz;

	multi = 1.0 / mass;
	A[0][0] = multi*Ixx;   A[0][1] = multi*Ixy;    A[0][2] = multi*Ixz;

	A[1][0] = multi*Ixy;   A[1][1] = multi*Iyy;    A[1][2] = multi*Iyz;

	A[2][0] = multi*Ixz;   A[2][1] = multi*Iyz;    A[2][2] = multi*Izz;


	U[0][0] = U[0][1] = U[0][2] = 0;
	U[1][0] = U[1][1] = U[1][2] = 0;
	U[2][0] = U[2][1] = U[2][2] = 0;

	j = Tridiagonalize_Symmetric_Matrix(*A, 0, diagonal, off_diagonal, *U, 3);
	i = QL_Tridiagonal_Symmetric_Matrix(diagonal, off_diagonal, *U, 3, 3000);
	if (i == -1 || j == -1)
	{
		printf("FAIL \n");
	}


	e1 = diagonal[0]; e2 = diagonal[1]; e3 = diagonal[2];


	if (e1 >= e2 && e1 >= e3)
	{
		ex = U[0][0]; ey = U[1][0]; ez = U[2][0];
		//printf("e1 \n");
	}
	else
	{
		if (e2 >= e3)
		{
			ex = U[0][1]; ey = U[1][1]; ez = U[2][1];
			//printf("e2 \n");
		}
		else
		{
			ex = U[0][2]; ey = U[1][2]; ez = U[2][2];
			//printf("e3 \n");
		}
	}

	mag2 = ex*ex + ey*ey + ez*ez;
	mag2 = sqrt(mag2);

	ex = ex / mag2; ey = ey / mag2; ez = ez / mag2;


	//////// find egen val and vect of inertia tensor...done////////////////////////////////////




	//cross product a X e. cross product is the axis of rotion////
	cx = (a_y*ez) - (a_z*ey);
	cy = (a_z*ex) - (a_x*ez);
	cz = (a_x*ey) - (a_y*ex);

	mag3 = cx*cx + cy*cy + cz*cz;
	mag3 = sqrt(mag3);

	cx = -cx / mag3; cy = -cy / mag3; cz = -cz / mag3;
	if (mag3 == 0)
	{
		cx = 0; cy = 0; cz = 0;
	}
	/////////////////////////////////////////////////////////////


	// dot product a * e  and angle difference //////////
	dot = a_x*ex + a_y*ey + a_z*ez;

	ang = 1 * acos(dot);

	//printf("e_vals (%lf, %lf ,%lf) \n", e1,e2,e3);
	//printf("a X e (%lf, %lf ,%lf) \n", cx, cy, cz);
	//printf("e_vec (%lf, %lf ,%lf) ang: %lf\n", ex,ey,ez,(180*ang/pi));

	//	if(ang> (pi/2))
	//{
	//ang=(pi)-ang;
	//ang=-1.0*ang;
	//}
	///////////////////////////////////////////////////





	// set up rotion matrix //////////////////////////


	s = sin(ang);
	c = cos(ang);
	u = (1 - cos(ang));
	V = cos(ang / 2);

	R[0][0] = c + cx*cx * u;
	R[0][1] = cx*cy * u - cz * s;
	R[0][2] = cx*cz * u + cy * s;

	R[1][0] = cx*cy * u + cz * s;
	R[1][1] = c + cy*cy * u;
	R[1][2] = cy*cz * u - cx * s;

	R[2][0] = cx*cz * u - cy * s;
	R[2][1] = cy*cz * u + cx * s;
	R[2][2] = c + cz*cz * u;
	///////////////////////////////////////////////		

	//	printf("U\n");
	//printf("%0.2lf   %0.2lf   %0.2lf \n",U[0][0],U[0][1],U[0][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",U[1][0],U[1][1],U[1][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",U[2][0],U[2][1],U[2][2]);
	//printf("\n");

	//printf("R\n");
	//printf("%0.2lf   %0.2lf   %0.2lf \n",R[0][0],R[0][1],R[0][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",R[1][0],R[1][1],R[1][2]);
	//printf("%0.2lf   %0.2lf   %0.2lf \n",R[2][0],R[2][1],R[2][2]);
	//printf("\n");


	/// time to overwrite x,y,z //////////////////////////////////////////////////////////////




	i = lead;
	while (i != -1)
	{
		Px = (R[0][0] * (x[i] - xc)) + (R[0][1] * (y[i] - yc)) + (R[0][2] * (z[i] - zc));

		Py = (R[1][0] * (x[i] - xc)) + (R[1][1] * (y[i] - yc)) + (R[1][2] * (z[i] - zc));

		Pz = (R[2][0] * (x[i] - xc)) + (R[2][1] * (y[i] - yc)) + (R[2][2] * (z[i] - zc));

		x[i] = Px;
		y[i] = Py;
		z[i] = Pz;

		i = clustnext[i];

	}



	//sprintf(filename, "clust_%d", lead);
	//output = fopen(filename, "w");

	//i = lead;

	//while (i != -1)
	//{

	//	fprintf(output, "%lf %lf %lf \n", x[i], y[i], z[i]);

		//i = clustnext[i];
	//}

	//fclose(output);



	return mass;


}
///////////////////////////////////////////////////////////////////////////

double get_length(int lead, double x[], double y[], double z[], double X_length[], double Y_length[], double Z_length[], int clustnext[], int clustlead[], int clustnum[], double rmsX[], double rmsY[], double rmsZ[])
{

	int i;
	double x_min, y_min, z_min, x_max, y_max, z_max, Rg, mass, r, rms_x, rms_y, rms_z;

	i = clustnum[lead];

	x_max = 0.0; y_max = 0.0; z_max = 0.0;

	x_min = 10000000;
	y_min = 10000000;
	z_min = 10000000;

	Rg = 0; mass = 0; rms_x = 0; rms_y = 0; rms_z = 0;

	//............  ...........///
	while (i != -1)
	{
		if (x[i] >= x_max)
		{
			x_max = x[i];
		}

		if (y[i] >= y_max)
		{
			y_max = y[i];
		}

		if (z[i] >= z_max)
		{
			z_max = z[i];
		}



		if (x[i] <= x_min)
		{
			x_min = x[i];
		}
		if (y[i] <= y_min)
		{
			y_min = y[i];
		}
		if (z[i] <= z_min)
		{
			z_min = z[i];
		}


		mass++;
		i = clustnext[i];

	}


	i = clustnum[lead];
	while (i != -1)
	{
		rms_x = rms_x + x[i] * x[i];
		rms_y = rms_y + y[i] * y[i];
		rms_z = rms_z + z[i] * z[i];

		r = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
		r = sqrt(r);

		Rg = Rg + (r)*(r);


		i = clustnext[i];
	}

	Rg = Rg / mass;
	Rg = sqrt(Rg);

	rms_x = rms_x / mass;
	rms_y = rms_y / mass;
	rms_z = rms_z / mass;

	rmsX[lead] = sqrt(rms_x);
	rmsY[lead] = sqrt(rms_y);
	rmsZ[lead] = sqrt(rms_z);


	X_length[lead] = x_max - x_min;
	Y_length[lead] = y_max - y_min;
	Z_length[lead] = z_max - z_min;



	i = clustnum[lead];
	while (i != -1)
	{
		x[i] = x[i] - x_min;
		y[i] = y[i] - y_min;
		z[i] = z[i] - z_min;

		i = clustnext[i];
	}

	return Rg;

}
double get_area(int i, double x[], double y[], double z[], double X_length[], double Y_length[], double Z_length[], int clustnext[], int clustlead[], int clustnum[], double RG_2d[])
{
	int lead, **grid, k, j, nx;
	int x_len, y_len, z_len, N;
	double count, Rg_2d, x_cm, z_cm ;

	Rg_2d = 0; 
	x_cm = 0; 
	z_cm = 0; 
	N = 0; 

	lead = clustnum[i];

	x_len = (int)X_length[i] + 1;
	y_len = (int)Y_length[i] + 1;
	z_len = (int)Z_length[i] + 1;


// sets aside memory ///////////////////////////////////////////////////////////////
	grid = (int**)malloc(x_len*sizeof (int*));
	if (grid == NULL)
	{
		printf("grid");
	}
	for (k = 0; k<x_len; k++)
	{
		grid[k] = (int*)malloc(z_len*sizeof(int));
		if (grid[k] == NULL)
		{
			printf("grid[%d]", k);
		}

	}
////////////////////////////////////////////////////////////////////////////
	
// sets grid to zero .................................................	
	for (k = 0; k<x_len; k++)
	{
		for (j = 0; j<z_len; j++)
		{
			grid[k][j] = 0;
		}
	}
// ....................................................................

// sets projected area // ...............................................
	nx = lead;
	while (nx != -1)
	{
		grid[(int)x[nx]][(int)z[nx]] = 1;

		nx = clustnext[nx];
	}
/////////////////////....................................................

	count = 0;

	//FILE *output;
	//output = fopen("2d_grid.txt", "w");

	for (k = 0; k<x_len; k++)
	{
		for (j = 0; j<z_len; j++)
		{
			if (grid[k][j] == 1)
			{
				x_cm = x_cm + k;
				z_cm = z_cm + j;
				N++;

				//fprintf(output, "%d %d \n", k, j);  
				if (mode == 0)
				{
					count = count + 0.7854;
				}

				if (mode == 1)
				{
					count = count++;
				}
			}
		}
	}


	x_cm = x_cm / N;
	z_cm = z_cm / N; 

	for (k = 0; k<x_len; k++)
	{
		for (j = 0; j<z_len; j++)
		{
			if (grid[k][j] == 1)
			{
				Rg_2d = Rg_2d + (k - x_cm)*(k - x_cm) + (j - z_cm)*(j - z_cm); 
			}
		}
	}

	Rg_2d = Rg_2d / (1.0*N); 

	Rg_2d = pow(Rg_2d, 0.5);

	RG_2d[i] = Rg_2d; 

	count = count /(1.0*pi);

	count = sqrt(count);

	//fclose(output);


	free(grid);

	return count;
}