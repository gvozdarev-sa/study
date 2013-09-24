#include <stdio.h>
#include <omp.h>

#define N 100000

#include "Camera.hpp"

int main(int argc, char *argv[])
{
    Camera cam( "1.png");


    int w = cam.GetWidth( ),
        h = cam.GetHeight( );

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

    cam.SaveImage( "2.bmp");


    return 0;
}
