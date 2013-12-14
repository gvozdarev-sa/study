#include <iostream>
#include <sstream>
#include <cmath>
#include <omp.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>


#include <boost/numeric/ublas/operations.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Timer.hpp"
#include "Image.hpp"

namespace ublas = boost::numeric::ublas;
typedef double  real_t;


typedef ublas::matrix< real_t> Matrix;
typedef ublas::vector< real_t> Vector;


void Solve_1_thread( const Matrix & A, const Vector & b, Vector & x);
void Solve_reduction( const Matrix & A, const Vector & f, Vector & x);

real_t f1( real_t x);
real_t f( real_t x, Vector & y);
const real_t FROM = -10.0;
const real_t TO   =  10.0;


//const int N = 15;
//const int N = 31;
//const int N = 63;
//const int N = 127;
//const int N = 255;
//const int N = 511;
//const int N = 2047;
//const int N = 4095;
const int N = 8191;


const real_t h = ( TO - FROM) / N;

const real_t a = sqrt( 2);
const real_t b = sqrt( 2);

void SetA( Matrix & A)
{
    A( 0,1) = 1.0;
    for( int i = 1; i < N-1; i++)
    {
        A( i, i + 1) = -2.0;
        A( i, i + 2) =  1.0;
        A( i, i + 0) =  1.0;
    }

    A( N-1, N)     = 1.0;
}

inline real_t f( real_t x){    return 1000000 * ( x - x * x * x);}
//inline real_t f( real_t x){    return  exp( -x);}


real_t F( int x, const Vector & y)
{
    Vector l;
    const real_t c_1_12 = 1.0 / 12.0;
    real_t f_    = f( y( x    ));
    real_t f_p_1 = f( y( x + 1));
    real_t f_m_1 = f( y( x - 1));

    real_t res = h * h * ( f_ + c_1_12 * ( f_m_1 - 2 * f_ + f_p_1 ));
    return res;
}


void SetXApprox( Vector & y)
{
    y( 1)   = a;
    for( int i = 2; i < N; i++)
    {
        y( i) = a + ( b - a) * ( real_t)( i - 1) / N;
    }
    y( N) = b;
}

void SetY( Vector & y)
{
    Vector y_c( y);

    y( 1)   =  a;
    for( int i = 2; i < N; i++)
    {
        y( i) = F( i, y_c);
    }
    y( N) =  b;
}



bool IsSim( const Vector & a, const Vector & b)
{
    real_t eps = 0.001;
    for ( int i = 0; i < N; i++)
    {
        if ( fabs( a( i) - b( i)) > eps)
        {
            return false;
        }
    }
    return true;
}

int main ( int argc, char * argv[ ])
{
    Matrix A( N, N + 2);
    Vector x     ( N + 2),
           x_prev( N + 2),
           y     ( N + 2);

    SetA( A);

    // первое приблежение
    SetXApprox( x);
    y = x;
    SetY( y);

    x_prev = x;
    Solve_reduction( A, y, x);

    int i = 3;
    while( ! IsSim( x, x_prev))
    {
        x_prev = y = x;

        SetXApprox( x);
        y = x;
        SetY( y);
        Solve_reduction( A, y, x);
//        Solve_1_thread( A, y, x);

        i++;
    }

    for ( int i = 1; i < x.size() - 1; i++)
    {
        std::cout << FROM + (real_t)i /N*( TO - FROM)  << "\t" << x( i) << "\n";
    }
}




bool InvertMatrix( const Matrix & input_m, Matrix & inverse_m)
{
    using namespace ublas;
    typedef permutation_matrix< std::size_t> pmatrix;
    // create a working copy of the input
    Matrix A( input_m);
    // create a permutation matrix for the LU-factorization
    pmatrix pm( A.size1());

    // perform LU-factorization
    int res = lu_factorize( A, pm);
    if( res != 0 ) return false;

    // create identity matrix of "inverse"
    inverse_m.assign( ublas::identity_matrix< real_t>( A.size1( )));

    // backsubstitute to get the inverse
    lu_substitute( A, pm, inverse_m);

    return true;
}

void Solve_1_thread( const Matrix & A, const Vector & b, Vector & x)
{
    ublas::slice col( 1, 1, A.size2( )-2);
    ublas::slice row( 0, 1, A.size1( ));

    Matrix __A( A);
    Vector __b( b);
    ublas::matrix_slice< Matrix> _A( __A, row,col);
    ublas::vector_slice< Vector> _f( __b, col);

//    std::cout << _A.size1( ) << " " << _A.size2( ) << "!!!\n";
//    std::cout << _A << _f << std::endl;

    Matrix A_inv( _A.size1( ), _A.size2( ));


    InvertMatrix( _A, A_inv);

    x = ublas::prod( A_inv, _f);
}

void ReductBlockBackward( Matrix & A, Vector & f, int i, Matrix & new_A, Vector & new_f)
{
    ublas::slice col( i, 1, 5);
    ublas::slice row( i, 1, 3);
    ublas::matrix_slice< Matrix> _A( A, row,col);
    ublas::vector_slice< Vector> _f( f, row);

    new_A.resize( 3, 5, false);

    // first row minus second
    new_A( 0,0) = ( _A( 0,0) * _A( 1,1)) - ( _A( 1,0) * _A( 0,1));
    new_A( 0,1) = ( _A( 0,1) * _A( 1,1)) - ( _A( 1,1) * _A( 0,1));
    new_A( 0,2) = ( _A( 0,2) * _A( 1,1)) - ( _A( 1,2) * _A( 0,1));
    new_A( 0,3) = ( _A( 0,3) * _A( 1,1)) - ( _A( 1,3) * _A( 0,1));
    new_A( 0,4) = ( _A( 0,4) * _A( 1,1)) - ( _A( 1,4) * _A( 0,1));
    new_f( 0)   = ( _f( 0)   * _A( 1,1)) - ( _f( 1)   * _A( 0,1));

   // third row minus second
    new_A( 2,0) = ( _A( 2,0) * _A( 1,3)) - ( _A( 1,0) * _A( 2,3));
    new_A( 2,1) = ( _A( 2,1) * _A( 1,3)) - ( _A( 1,1) * _A( 2,3));
    new_A( 2,2) = ( _A( 2,2) * _A( 1,3)) - ( _A( 1,2) * _A( 2,3));
    new_A( 2,3) = ( _A( 2,3) * _A( 1,3)) - ( _A( 1,3) * _A( 2,3));
    new_A( 2,4) = ( _A( 2,4) * _A( 1,3)) - ( _A( 1,4) * _A( 2,3));
    new_f( 2)   = ( _f( 2)   * _A( 1,3)) - ( _f( 1)   * _A( 2,3));
}

/*
 * INPUT:
 *  Matrix A 3 x 5
 *  Vector f 3
 * FROM:
 * a1*x1  + b1*x2 + c1*x3                 = f1;
 *          a2*x2 + b2*x3 + c2*x4         = f2;
 *                  a3*x3 + b3*x4 + c3*x5 = f3;
 * TO:
 * a*x1           + b*x3  +       + c*x5  = f;
 */

void ReductBlock( Matrix & A, Vector & f, int i, Matrix & new_A, Vector & new_f)
{
    ublas::slice col( i, 1, 5);
    ublas::slice row( i, 1, 3);
    ublas::matrix_slice< Matrix> _A( A, row,col);

    ublas::slice vec( i + 1, 1,3);
    ublas::vector_slice< Vector> _f( f, vec);

    //std::cout << i << " " << _A << "\t" << _f <<  std::endl;

    new_A.resize( 3, 5, false);

    // second row minus first
    new_A( 1,0) = ( ( _A( 1,0) * _A( 0,1)) - ( _A( 0,0) * _A( 1,1))) /_A( 1,1) ;
    new_A( 1,1) = ( ( _A( 1,1) * _A( 0,1)) - ( _A( 0,1) * _A( 1,1))) /_A( 1,1) ;
    new_A( 1,2) = ( ( _A( 1,2) * _A( 0,1)) - ( _A( 0,2) * _A( 1,1))) /_A( 1,1) ;
    new_A( 1,3) = ( ( _A( 1,3) * _A( 0,1)) - ( _A( 0,3) * _A( 1,1))) /_A( 1,1) ;
    new_A( 1,4) = ( ( _A( 1,4) * _A( 0,1)) - ( _A( 0,4) * _A( 1,1))) /_A( 1,1) ;
    new_f( 1)   = ( ( _f( 1)   * _A( 0,1)) - ( _f( 0)   * _A( 1,1))) /_A( 1,1) ;

    //std::cout << i << "sec - fir" << new_A << "\t" << new_f <<  std::endl;
    Matrix __A( new_A);
    Vector __f( new_f);
    // second row minus third
    new_A( 1,0) = (( __A( 1,0) * _A( 2,3)) - ( _A( 2,0) * __A( 1,3))) / __A( 1,3);
    new_A( 1,1) = (( __A( 1,1) * _A( 2,3)) - ( _A( 2,1) * __A( 1,3))) / __A( 1,3);
    new_A( 1,2) = (( __A( 1,2) * _A( 2,3)) - ( _A( 2,2) * __A( 1,3))) / __A( 1,3);
    new_A( 1,3) = (( __A( 1,3) * _A( 2,3)) - ( _A( 2,3) * __A( 1,3))) / __A( 1,3);
    new_A( 1,4) = (( __A( 1,4) * _A( 2,3)) - ( _A( 2,4) * __A( 1,3))) / __A( 1,3);
    new_f( 1)   = (( __f( 1)   * _A( 2,3)) - ( _f( 2)   * __A( 1,3))) / __A( 1,3);

    //std::cout << i << new_A << "\t" << new_f <<  std::endl;
}

void Solve_reduction( const Matrix & A, const Vector & f, Vector & x)
{
    using std::cout;
    bool T = true;
    int stride = 1;
    if( N == 3)
    {
        T = false;
    }
    Matrix curr_m = A;
    Vector curr_f = f;

    Matrix next_m;
    Vector next_f;

    std::vector< Matrix> Ai;
    std::vector< Vector> fi;

    Ai.push_back( A);
    fi.push_back( f);

//    cout  << "\n???\n\n" << f << "\n\n???\n";

        while( T)
        {
            int curr_size = curr_m.size1( );
            int next_size = ( curr_size + 1) / 2 - 1;

            next_m.resize( next_size, next_size + 2, false);
            next_f.resize( next_size + 2, false);

            next_m.clear( );
            next_f.clear( );


            #pragma omp parallel for shared(A,f,x,Ai,fi,curr_m,curr_f,next_m,next_f)
            for( int i = 0; i < curr_size - 2; i+=2)
            {
                Matrix block_m( 3,5, false);
                Vector block_f( 3, false);

               ReductBlock( curr_m, curr_f, i, block_m, block_f);
               next_m( i / 2, i / 2+ 0) = block_m( 1,0);
                next_m( i / 2, i / 2+ 1) = block_m( 1,2);
                next_m( i / 2, i / 2+ 2) = block_m( 1,4);

                next_f( i / 2 + 1) = block_f( 1);
            }
           if ( next_size == 3)
            {
                T = false;
            }
            else
            {
                curr_f = next_f;
                curr_m = next_m;

                Ai.push_back( next_m);
                fi.push_back( next_f);

                stride*=2;
            }
//            std::cout << "Next :\n" << next_m << " " << next_f << std::endl;
        }
    Matrix block_m( 3,5, true);
    Vector block_f( 3, true);

//    cout << "???" << next_m << " " << next_f << "\n";
    ReductBlock( next_m, next_f, 0, block_m, block_f);
    int j = ( N + 1)/2;
    x( j + 0) = block_f( 1) / block_m( 1,2);

    x( 1*( N +1) /4) = ( next_f( 3) - x( j + 0) * next_m( 2,2)) / next_m( 2,3);
    x( 3*( N +1) /4) = ( next_f( 1) - x( j + 0) * next_m( 0,2)) / next_m( 0,1);
    x( 1) = a;
    x( N) = b;


//    cout << "BEFORE BACKWARD!!! 1 !!! " << stride << " "  << x << "\n";


   {
        Matrix & curr_m = Ai.back( );
        Vector & curr_f = fi.back( );

        j /= stride;

        x( stride * ( j + 1)) = ( curr_f( 5) - x( stride * ( j + 0)) * curr_m( 4,4) - x( stride * ( j + 2)) * curr_m( 4,6)) / curr_m( 4,5);
        x( stride * ( j - 1)) = ( curr_f( 3) - x( stride * ( j + 0)) * curr_m( 2,4) - x( stride * ( j - 2)) * curr_m( 2,2)) / curr_m( 2,3);

        x( stride * ( j + 3)) = ( curr_f( 6) - x( stride * ( j + 2)) * curr_m( 5,6) - x( stride * ( j + 1)) * curr_m( 5,5)) / curr_m( 5,7);
        x( stride * ( j - 3)) = ( curr_f( 2) - x( stride * ( j - 2)) * curr_m( 1,2) - x( stride * ( j - 1)) * curr_m( 1,1)) / curr_m( 1,1);

        Ai.pop_back( );
        fi.pop_back( );
        stride /= 2;
    }
    T = true;

//    cout << "\n\nBEFORE BACKWARD!!! 2 !!! " << x << " "  << stride<< "\n\n";

        while ( T)
        {
            Matrix & curr_m = Ai.back( );
            Vector & curr_f = fi.back( );

            int curr_size = curr_f.size( );

            //cout << "curr_m " << curr_m <<  "\ncurr_size " << curr_size  <<std::endl;

            #pragma omp parallel for shared(A,f,x,Ai,fi,curr_m,curr_f)
            for ( int i = 0; i <  curr_size - 2; i+=2)
            {
//                cout << "i " << i << "\t";
                ublas::slice col( i, 1, 5);
                ublas::slice row( i, 1, 3);
                ublas::matrix_slice< Matrix> _A( curr_m, row,col);
                ublas::vector_slice< Vector> _f( curr_f, row);

                ///cout << "_A " << " "  <<_A << "\n";
                x( (i + 1) * stride) = ( _f( 1) - x( ( i + 0) * stride ) *  _A( 0,0) - x( ( i + 2) * stride) * _A( 0,2)) / _A( 0,1);
            }

            if ( stride == 1)
            {
                T = false;
            }
            else
            {
                Ai.pop_back( );
                fi.pop_back( );
                stride /= 2;
            }
        }
//    cout << "\n\nAFTER BACKWARD " << x << " "  << stride<< "\n\n";
}


