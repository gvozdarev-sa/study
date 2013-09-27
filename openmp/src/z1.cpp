#include <iostream>
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
    Image cam( argv[ 1]);


    int w = cam.GetWidth( ),
        h = cam.GetHeight( );



    Timer timer;
    timer.Start( );
    #pragma omp paraller for
    {
        for ( int i = 0; i < h; i++)
        {
            for ( int j = 0; j < w; j++)
            {
                Color col = cam.GetPixel( i,j);
                Color ncol;

                ncol.r = ncol.g = ncol.b = 0.299f * col.r + 0.587f * col.g + 0.114f * col.b;
                cam.SetPixel( i, j, ncol);
            }
        }
    }
    timer.Stop( );

    cam.SaveImage( "2.bmp");

    std::cout << "Time" << timer.GetTime( ) << std::endl;

    return 0;
}
