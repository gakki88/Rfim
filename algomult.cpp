/* FILE NAME: algomult.cpp
 * 
 * PURPOSE: This file contains function which implement the multiplicative algorithm. (Yaming Yu 2010)
 *
 * DEVELOPERS: PFIM group, pfim@inserm.fr
 *			       Yuxin Tang, INSERM UMR 1137 BIPID team
 *
 * INPUT VARIABLES FOR THE ENTRY POINT FUNCTION (multiplicativeAlgorithm) FROM R: 
 * Name		      		Type		        	Description
 * ----------       ----------      ----------------------------
 * originMatrices   double*         a list of fisher information matrices(FIM) computed by R functions 
 *                                  and read as a vector of values column by column in the following functions
 * weights          double*         a vector of weight corresponding to input FIMs and will evolve in the end
 *                                  of each iteration of algorithm
 * dim			      	int*			      the dimension of matrices (all matrices have same dimension)
 * n		        		int*	      		total number of input matrices
 * dimA		       		int*	      		the dimension of block A of matrix( all matrices have same dimension of 
 *                                  block A)
 * lambda	      		double*		      a parameter of algorithm between 0 and 1
 * delta		      	double*		    	a precision of stop criterion of algorithm 
 * iteration	     	int*		      	the algorithm maximum number of iteration


 ALGORITHM
 * 	associate each origin FIM input from R with a weight value

 *	sum all weighted matrices as one matrix

 *	calculate determinant of the summed matrix, and also separately the determinant of block A and block B
    of the summed matrix

 *	calculate a coefficient as (D-criterion of the summed matrix) divided by (the total dimension of a matrix)
    
 *  inverse the summed matrix
		
 *  multiply each element of the inverse matrix with the coefficient, the matrix obtained is the derivative 
    of phi
		
 *  multiply every origin FIM input from R with the derivative of phi, then compute the trace of each matrix 
    product
		
 *  evolve weights by multiply corresponding (trace^lambda), and divided by its sum

 *  algorithm stops when criterion about delta is satisfied
 */


#include <utility>
#include <limits.h>
#include <math.h>
#include "algomult.h"
 
 
#ifdef _WINDOWS

#include <Windows.h>

// windows dll entry point
BOOL APIENTRY DllMain( HMODULE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

#elif


#endif  

/* FUNCTION NAME: sum
 *
 * DESCRIPTION: This function adds up a list of matrices as one matrix.
 * 
 * INPUT:
 * Name		      		Type		      	Description
 * ----------       ----------      ----------------------------
 * matrices		    	double*	    		a vector of n matrices column by colunm of dimension dim
 * size			      	int		      		total number of element of a matrix (size = dim * dim)
 * n		        		int		      		total number of matrices
 *
 * RETURN VALUE:
 * Name		       		Type		      	Description
 * ----------       ----------      ----------------------------
 * result		      	double*			    one matrix where each element equals the sum of element of input matrices at
 *									                that position 
 */
void sum(double *matrices, int size, int n, double *result)
{
	int i, j;
	// initialise the result matrix 
	// copy the first matrix in the list to the result matrix
	for (j = 0; j < size; j++) {
		*result = *matrices;
		result++;
		matrices++;
	}
	result -= size;

	// add rest matrices in the list to the result matrix
	for (i = 1; i < n; i++) {
		for (j = 0; j < size; j++){
			*result += *matrices;
			result++;
			matrices++;
		}
		result -= size;
	}
}

/* FUNCTION NAME:  MatrixTimesWeight
 *
 * DESCRIPTION: This function multiplies each element in the matrix by a weight.
 *
 * INPUT:
 * Name		      		Type		      	Description
 * ----------       ----------      ----------------------------
 * matrices		      double*	    		a vector of n matrices column by colunm of dimension dim
 * weights	    		double*		    	n values, one for each matrix
 * size				      int		      		total number of element of a matrix (size = dim * dim)
 * n			        	int		      		total number of matrices
 *
 * RETURN VALUE:
 * Name		      		Type		      	Description
 * ----------       ----------      ----------------------------
 * result		      	double*		      weighted matrices, same number of matrices as input
 */
void MatrixTimesWeight(double *matrices, double *weights, int size, int n, double *result)
{
	int size_n = size * n;

	for (int i = 0; i < size_n; i++)
	{
		*result = (*matrices) * (*weights);
		matrices++;
		result++;
		if ((i + 1) % size == 0)
			weights++;
	}

}


/* FUNCTION NAME:  decompositionLU
 *
 * DESCRIPTION: This function decomposite a matrix as a lower traiangular matrix and an upper triangular matrix.
 *
 * INPUT:
 * Name			      	Type		      	Description
 * ----------       ----------      ----------------------------
 * matrix		      	double*		      a symmetric matrix, can be divided into two blocks independent non-zero 
 *                                  (block A,B)
 * dim			        int		      		dimension of an entire matrix
 * dimA			      	int		      		dimension of block A of the matrix (dimension of block B = dim - dimA)
 *
 * RETURN VALUE:
 * Name		      		Type		      	Description
 * ----------       ----------      ----------------------------
 * L			        	double*		    	LU decomposition result -- lower matrix
 * U			        	double*		    	LU decomposition result -- upper matrix
 */
void decompositionLU(double *matrix, double *L, double *U, int dim, int dimA) {
	int i = 0, j, k, i_dim, i_j_dim, jj_dim;

	// decomposition LU for block A  
	do {
		i_dim = i * dim;
		for (j = 0; j < i; j++) {
			jj_dim = j * dim + j;
			if (U[jj_dim] == 0) {
				i = dimA;
				break;
			}
			i_j_dim = i_dim + j;
			L[i_j_dim] = matrix[i_j_dim];
			for (k = 0; k < j; k++) {
				L[i_j_dim] -= L[i_dim + k] * U[k * dim + j];
			}
			L[i_j_dim] /= U[jj_dim];

		}
		for (j = i; j < dimA; j++) {
			i_j_dim = i_dim + j;
			U[i_j_dim] = matrix[i_j_dim];
			for (k = 0; k < i; k++) {
				U[i_j_dim] -= L[i_dim + k] * U[k * dim + j];
			}
		}
		i++;
	} while (i < dimA);

	// decomposition LU for block B, if dimA < dim
	while (i < dim) {
		i_dim = i * dim;
		for (j = dimA; j < i; j++) {
			jj_dim = j * dim + j;
			if (U[jj_dim] == 0) {
				i = dim;
				break;
			}
			i_j_dim = i_dim + j;
			L[i_j_dim] = matrix[i_j_dim];
			for (k = dimA; k < j; k++) {
				L[i_j_dim] -= L[i_dim + k] * U[k * dim + j];
			}
			L[i_j_dim] /= U[jj_dim];
			
		}
		for (j = i; j < dim; j++) {
			i_j_dim = i_dim + j;
			U[i_j_dim] = matrix[i_j_dim];
			for (k = dimA; k < i; k++) {
				U[i_j_dim] -= L[i_dim + k] * U[k * dim + j];
			}
		}
		i++;
	};


}


/* FUNCTION NAME:  inverse
 *
 * DESCRIPTION: This function inverses a symmetric matrix
 *
 * INPUT:
 * Name		      		Type	      		Description
 * ----------       ----------      ----------------------------
 * m		        		double*	      	a symmetric matrix to be inversed, can be divided into two independent 
 *                                  blocks (A, B)
 * dim		      		int			      	dimension of an entire matrix
 * dimA		      		int		      		dimension of block A of the matrix (dimension of block B = dim - dimA)
 * detA			      	double		      determinant of block A
 * detB		      		double	    		determinant of block B
 *
 * RETURN VALUE:
 * Name	      			Type		      	Description
 * ----------       ----------      ----------------------------
 * invM			      	double*		    	inversed matrix
 */
void inverse(double *m, int dim, double *invM, int dimA, double detA, double detB)
{
	int i, j, ii, jj, i1, j1, dim1, i1_dim, ii_dim, k;
	double det, *c, *d, *l, *u;

	// prepare matrix L and U 
	l = (double *)calloc(dim * dim, sizeof(double));
	u = (double *)calloc(dim * dim, sizeof(double));
	// initialise matrix L with all diagonal
	for (i = 0; i < dim; i++)
		l[i*dim + i] = 1;

	//calculate inverse for block A
	dim1 = dimA - 1;
	c = (double *)calloc((dim1) *(dim1), sizeof(double));

	for (j = 0; j < dimA; j++) {
		for (i = 0; i < dimA; i++) {
			i1 = 0;
			for (ii = 0; ii < dimA; ii++) {
				if (ii == i)
					continue;
				i1_dim = i1 * dim1;
				ii_dim = ii * dim;
				j1 = 0;
				for (jj = 0; jj < dimA; jj++) {
					if (jj == j)
						continue;
					*(c + i1_dim + j1) = *(m + ii_dim + jj);
					j1++;
				}
				i1++;
			}

			decompositionLU(c, l, u, dim1, dim1);
			det = 1;
			for (k = 0; k < dim1; k++)
				det *= u[k * dim1 + k];
			invM[i + j * dim] = pow(-1.0, i + j + 2.0) * (det) / (detA);
		}
	}

	//calculate inverse for block B
	dim1 = dim - dimA - 1;
	d = (double *)calloc((dim1) *(dim1), sizeof(double));

	for (j = dimA; j < dim; j++) {

		for (i = dimA; i < dim; i++) {
			i1 = 0;
			for (ii = dimA; ii < dim; ii++) {
				if (ii == i)
					continue;
				i1_dim = i1 * dim1;
				ii_dim = ii * dim;
				j1 = 0;
				for (jj = dimA; jj < dim; jj++) {
					if (jj == j)
						continue;
					*(d + i1_dim + j1) = *(m + ii_dim + jj);
					j1++;
				}
				i1++;
			}

			decompositionLU(d, l, u, dim1, dim1);
			det = 1;
			for (k = 0; k < dim1; k++)
				det *= u[k * dim1 + k];
			invM[i + j * dim] = pow(-1.0, i + j + 2.0) * (det) / (detB);

		}
	}
	free(c);
	free(d);
	free(l);
	free(u);
}


/* FUNCTION NAME:  computeTrace
 *
 * DESCRIPTION: This function multiplies each matrix in the list by the same  matrix, then compute each 
 *              product's trace
 *
 * INPUT:
 * Name		      		Type		      	Description
 * ----------       ----------      ----------------------------
 * matrices	      	double*	    		a list of matrices to produce with theMatrix
 * theMatrix	    	double*		    	the matrix to be multiplied with each matrices in the list
 * dim		      		int		      		dimension of an entire matrix
 * n			        	int			      	number of matrices in the list
 *
 * RETURN VALUE:
 * Name			      	Type		      	Description
 * ----------       ----------      ----------------------------
 * trace		        double*		    	the trace vector of each matrix product
 */
void computeTrace(double *matrices, double *theMatrix, int dim, int n, double *trace)
{
	int im, il, in, i;
	double *productMatrices = NULL;
	int size = dim * dim, i_size;
	// intermediate matrices
	productMatrices = (double *)calloc(n * size, sizeof(double));

	for (i = 0; i < n; i++) {
		*trace = 0;
		for (in = 0; in < dim; in++)
		{
			for (im = 0; im < dim; im++)
			{
				*productMatrices = 0;
				for (il = 0; il < dim; il++)
				{
					*productMatrices += (*matrices) * (*theMatrix);
					matrices++;
					theMatrix += dim;
				}
				matrices -= dim;
				theMatrix -= size;
				if (im == in)
					*trace += *productMatrices;
				productMatrices++;
				theMatrix++;
			}
			matrices += dim;
			theMatrix -= dim;
		}
		trace++;
	}
	productMatrices -= n * size;
	free(productMatrices);
}


/* FUNCTION NAME:  multiplicativeAlgorithm
 *
 * DESCRIPTION: This is the main function of multiplicative algorithm,which is exported and used in PFIM 
 *              implemented in R-S4
 *
 * INPUT:
 * Name	      			Type		      	Description
 * ----------       ----------      ----------------------------
 * originMatrices	  double*		      a list of matrices computed by R functions
 * weights			    double*	      	a vector of weight corresponding to each matrix in the originMatrices list
 * dim				      int			      	dimension of an entire matrix (all same for matrices in the list)
 * n			        	int*	      		number of matrices in the list
 * dimA			      	int*		      	dimension of block A 
 * lambda		      	double*		      convergence power of the algorithm (0, 1]
 * delta		      	double*		    	algorithm stop criterion parameter
 * iteration		    int*	      		algorithm maximum iteration
 *
 * RETURN VALUE:
 * Name			      	Type		      	Description
 * ----------       ----------      ----------------------------
 * weights		    	double*		    	the vector of weights evolved by algorithm
 * iteration	    	int*	      		real iteration run number
 */
void multiplicativeAlgorithm(double *originMatrices, double *weights, int *dim,
	int *n, int *dimA,	double *lambda, double *delta, int *iteration)
{
	double det, coef, sumw, sumwd, maxTrace, *weightedMatrices, *sumWeightedMatrices,
		*inverseSumWeightMat, *Dphi, *tr, *l, *u, detA, detB;
	int idim = *dim, in = *n, idimA = *dimA, isize = idim * idim, i, it;

	weightedMatrices = (double *)calloc(in * isize, sizeof(double));
	if (!weightedMatrices) {
		it = -2;
		return;
	}
	sumWeightedMatrices = (double *)calloc(isize, sizeof(double));
	inverseSumWeightMat = (double *)calloc(isize, sizeof(double));
	Dphi = (double *)calloc(isize, sizeof(double));
	tr = (double *)calloc(in, sizeof(double));
	l = (double *)calloc(isize, sizeof(double));
	u = (double *)calloc(isize, sizeof(double));

	// initialise diagonal elements of the low matrix for decomposition LU as 1
	for (i = 0; i < idim; i++)
		l[i*idim + i] = 1;

	for (it = 0; it < (*iteration); it++)
	{
		// associate each FIM calculated in R with a weight value
		MatrixTimesWeight(originMatrices, weights, isize, in, weightedMatrices);
		// sum all weighted matrices
		sum(weightedMatrices, isize, in, sumWeightedMatrices);

		//calculate determinant of Mw, and separately of block A and B of Mw
		detA = 1;
		detB = 1;
		decompositionLU(sumWeightedMatrices, l, u, idim, idimA);
		for (i = 0; i < idimA; i++)
			detA *= u[i * idim + i];
		for (i = idimA; i < idim; i++)
			detB *= u[i * idim + i];
		det = detA * detB;
		if (det < 0 || isnan(det)) {
			it = -1;
			break;
		}

		// calculate a coefficient as (D-criterion)/(dimension of matrix)
		coef = pow(det, (1.0 / idim)) / idim;
		// inverse the sum matrix of weighted matrices
		inverse(sumWeightedMatrices, idim, inverseSumWeightMat, idimA, detA, detB);
		// multiply each element of the inverse matrix with the coefficient,
		// the result obtained is the derivative of phi
		MatrixTimesWeight(inverseSumWeightMat, &coef, isize, 1, Dphi);
		// multiply every FIM calculated in R with the derivative of phi, 
		// then obtain the trace of each product
		computeTrace(originMatrices, Dphi, idim, in, tr);
		// evolve weights
		sumw = 0.0; // sum of weights
		sumwd = 0.0; // sum of lambda exponential weights
		for (i = 0; i < in; i++)
		{
			sumwd += weights[i] * tr[i];
			weights[i] *= pow(tr[i], *lambda);
			sumw += weights[i];
		}
		maxTrace = tr[0];
		for (i = 0; i < in; i++)
		{
			weights[i] /= sumw;    // make sure weights are proper fractions
			if (maxTrace < tr[i])
				maxTrace = tr[i];
		}
		
		if (maxTrace < (1 + *delta)*sumwd)    // stop criterion of the algorithm
		{
			printf("loop finished with iteration = %d\n", it);
			break;
		}

	}
	*iteration = it;


	free(weightedMatrices);
	free(sumWeightedMatrices);
	free(inverseSumWeightMat);
	free(Dphi);
	free(tr);
	free(l);
	free(u);
}