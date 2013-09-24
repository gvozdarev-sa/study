#include <iostream>
#include <omp.h>

#include "Image.hpp"

int main(int argc, char *argv[])
{
    if ( argc != 2)
    {
        std::cout << "Please specify image file" << std::endl;
    }
    Image cam( argv[ 1]);


    int w = cam.GetWidth( ),
        h = cam.GetHeight( );

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
    cam.SaveImage( "2.bmp");


    return 0;
}
