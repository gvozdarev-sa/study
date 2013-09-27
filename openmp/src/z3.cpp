#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

#include "Image.hpp"


class Timer
{
public:
    void Start( )
    {
        gettimeofday( &_start_time, NULL);
    }
    void Stop( )
    {
        gettimeofday( &_end_time, NULL);
    }
    double GetTime( )
    {
        return _end_time.tv_sec - _start_time.tv_sec + 0.000001 * ( _end_time.tv_usec - _start_time.tv_usec);
    }
private:
    struct timeval _start_time,
                   _end_time;

};

int main(int argc, char *argv[])
{
    if ( argc != 2)
    {
        std::cout << "Please specify image file" << std::endl;
    }

    Image I1( argv[ 1]),
          I2( argv[ 1]);


    int w = I1.GetWidth( ),
        h = I2.GetHeight( );



    Timer timer;
    timer.Start( );
    #pragma omp parallel
    {
        #pragma omp segments
        {
            #pragma omp segment
            #pragma omp for schedule(guided, 100)
            for ( int i = 0; i < h; i++)
            {
                for ( int j = 0; j < w; j++)
                {
                    Color col1 = I1.GetPixel( i, j  ),
                          col2 = I1.GetPixel( i, j+1);
                    Color ncol;

                    ncol.r = std::abs( col1.r - col2.r);
                    ncol.g = std::abs( col1.g - col2.g);
                    ncol.b = std::abs( col1.b - col2.r);

                    I1.SetPixel( i, j, ncol);
                }
            }
            #pragma omp segment
            #pragma omp for schedule(guided, 100)
            for ( int j = 0; j < w; j++)
            {
                for ( int i = 0; i < h; i++)
                {
                    Color col1 = I2.GetPixel( i  , j),
                          col2 = I2.GetPixel( i+1, j);
                    Color ncol;

                    ncol.r = std::abs( col1.r - col2.r);
                    ncol.g = std::abs( col1.g - col2.g);
                    ncol.b = std::abs( col1.b - col2.r);

                    I2.SetPixel( i, j, ncol);
                }
            }
        }
        #pragma omp for schedule(guided, 100)
        for ( int i = 0; i < h; i++)
        {
            for ( int j = 0; j < w; j++)
            {
                Color col1 = I1.GetPixel( i, j),
                      col2 = I2.GetPixel( i, j);
                Color ncol;

                ncol.r = 0.5 * col1.r + 0.5 * col2.r;
                ncol.g = 0.5 * col1.g + 0.5 * col2.g;
                ncol.b = 0.5 * col1.b + 0.5 * col2.b;

                I1.SetPixel( i, j, ncol);
            }
        }
    }

    timer.Stop( );

    I1.SaveImage( "2.jpg");

    std::cout << "Time : " << timer.GetTime( ) << std::endl;

    return 0;
}
