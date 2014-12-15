#include "matrix.hpp"
#include "MatrixOutputs.hpp"
#include <stdio.h>


int main()
{
	//Problem 1. LU Factorization
	{
		matrix m1, lu_m1, x(3), b; //Create two matrix objects. One for the input matrix, the other for the output of the LU Optimization

		int permutationvector[ 3 ];

		m1 = matrix( 3, 3 );
		//row 1
		m1( 0, 0 ) = 2.0 / 3.0;
		m1( 0, 1 ) = 1.0;
		m1( 0, 2 ) = 1.0 / 2.0;
		//row 2
		m1( 1, 0 ) = 1.0 / 2.0;
		m1( 1, 1 ) = 2.0 / 3.0;
		m1( 1, 2 ) = 1.0;
		//row 3
		m1( 2, 0 ) = 1.0;
		m1( 2, 1 ) = 1.0 / 2.0;
		m1( 2, 2 ) = 2.0 / 3.0;

		x( 0 ) = 3;
		x( 1 ) = 2;
		x( 2 ) = 1;

		lu_m1 = m1; //Copy the input matrix to the LU Optimization matrix
		LU( lu_m1, permutationvector ); //Perform LU Optimization on the lu_m1 object
		
		printf( "\nInput m1\n" );
		PrintMatrix( m1 );
		printf( "\n LU Factorization of m1\n" );
		PrintMatrix( lu_m1 );

		printf( "\nPermutation Matrix \n\n" );
		printf( "%d, %d, %d\n\n", permutationvector[ 0 ], permutationvector[ 1 ], permutationvector[ 2 ] );
		b = x;

		//Display b vector
		printf( "Vector b = \n" );
		PrintMatrix( b );

		// Solve for "x" and display solution.
		x = usolve( lu_m1, lsolve( lu_m1, permutate( b, permutationvector ) ) );
		printf( "Solution Vector x = \n" );
		PrintMatrix( x, 5 );

		printf( "Press enter to continue with Question 2" );
		getchar();

	}

	//Problem 2. 
	//Use code in the Matrix.hpp file to calculate the solutions to the matrix from question 1 using an iterative approach. 
	{
		matrix A, B, Identity_Mat, X, Iter, bIter, Y, e; //Declare matrix objects needed
		double mu; //µ variable used during calculations
		int	k; //loop counter
		mu = 0.3; //We will use an µ value of 0.3
		A = matrix( 3, 3 ); //Initialize the A matrix and then load it with the following values
		//row 1
		A( 0, 0 ) = 2.0 / 3.0;
		A( 0, 1 ) = 1.0;
		A( 0, 2 ) = 1.0 / 2.0;
		//row2
		A( 1, 0 ) = 1.0 / 2.0;
		A( 1, 1 ) = 2.0 / 3.0;
		A( 1, 2 ) = 1.0;
		//row3
		A( 2, 0 ) = 1.0;
		A( 2, 1 ) = 1.0 / 2.0;
		A( 2, 2 ) = 2.0 / 3.0;

		B = matrix( 3 ); //Initialize the B matrix and load it with the following values
		B( 0 ) = 3;
		B( 1 ) = 2;
		B( 2 ) = 1;

		Identity_Mat = eye( 3, 3 ); //Initialize the Identity_Mat object with the Matrix.hpp's eye function

		X = matrix( 3 );
		X( 0 ) = 0;
		X( 1 ) = 0;
		X( 2 ) = 0;

		//Do the actual computation of ( I - mu * A' * A )
		Iter = ( Identity_Mat - A.transpose() * A * mu);
		bIter = ( A.transpose() * B * mu );

		//Loop that performs the iteration and calculates the iterations
		for ( k = 0; k < 100; k++ )
		{
			//Performs the iteration
			Y = Iter * X + bIter;
			e = A * X - B;
			e = e.transpose() * e; //The error will be in location e(0)
			//Print X and the error
			printf( "X=\n" );
			PrintMatrix( X );
			printf( "\nError =\n" );
			PrintMatrix( e );
			X = Y; //Replace x with y so the iteration can continue
		}
	}

	printf( "Press enter to exit" );
	getchar();
}