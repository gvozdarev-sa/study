#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "Image.hpp"
#include "Timer.hpp"

// Медианный фильтр


bool compare_red( const Color & a, const Color &  b)
{
    return a.r < b.r;
}

bool compare_green( const Color&  a, const Color &  b)
{
    return a.g < b.g;
}

bool compare_blue( const Color &  a, const Color &  b)
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
        for ( int i = 0; i < h; i++)
        {
            for ( int j = 0; j < w; j++)
            {
                std::vector<Color>  array;
                array.reserve( R * R);
                Color col;

                for ( int ii = -R; ii < R; ii++)
                {
                    for ( int jj = -R; jj < R; jj++)
                    {
                        if ( ii * ii + jj * jj <= R * R)
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

    std::cout << "Time : " << timer.GetTotalTime( ) << std::endl;

    return 0;
}
