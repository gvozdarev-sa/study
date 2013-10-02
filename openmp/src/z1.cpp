#include <iostream>
#include <omp.h>

#include "Image.hpp"
#include "Timer.hpp"

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
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 100)
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

    cam.SaveImage( "2.jpg");

    std::cout << "Time : " << timer.GetTotalTime( ) << std::endl;

    return 0;
}
