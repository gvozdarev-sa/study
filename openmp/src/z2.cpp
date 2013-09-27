#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <omp.h>

#include "Image.hpp"

// Медианный фильтр


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

int compare_red( Color  a, Color   b)
{
    return a.r < b.r;
}

int compare_green( Color  a, Color   b)
{
    return a.g < b.g;
}

int compare_blue( Color  a, Color   b)
{
    return a.b < b.b;
}


int main(int argc, char *argv[])
{
    if ( argc < 2)
    {
        std::cout << "Please specify image file" << std::endl;
    }
    Image I1( argv[ 1]);
    Image I2( argv[ 1]);


    int R = 3;
    if ( argc == 3)
    {
        std::stringstream convert;
        convert << argv[ 2];
        convert >> R;
    }



    int w = I1.GetWidth( ),
        h = I1.GetHeight( );


    Timer timer;
    timer.Start( );
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 100)
        for ( int i = R; i < h - R; i++)
        {
            for ( int j = R; j < w - R; j++)
            {
                std::vector<Color>  array( 8);
                Color col;

                for ( int ii = 0; ii < R; ii++)
                {
                    for ( int jj = 0; jj < R; jj++)
                    {
                        if ( ii * ii + jj * jj < R * R + 0.1)
                        {
                            array.push_back( I1.GetPixel( i + ii,j + jj));
                        }
                    }
                }

                // sort by red
                std::sort( array.begin( ), array.end( ), compare_red);
                col.r = array[ array.size( ) / 2].r;

                std::sort( array.begin( ), array.end( ), compare_green);
                col.g = array[ array.size( ) / 2].g;

                std::sort( array.begin( ), array.end( ), compare_blue);
                col.b = array[ array.size( ) / 2].b;

                I2.SetPixel( i,j, col);
            }
        }
    }
    timer.Stop( );

    I2.SaveImage( "2.jpg");

    std::cout << "Time : " << timer.GetTime( ) << std::endl;

    return 0;
}
