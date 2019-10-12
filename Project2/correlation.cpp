#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include "stdlib.h"

#include "image.h"
using namespace std;

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);


int main(int argc, char *argv[]){
    int M, N, Q;
    bool type;
    stringstream factorQ(argv[3]);
    int factor=0;
    factorQ>>factor;

    // read image header
    readImageHeader(argv[1], N, M, Q, type);

    // allocate memory for the image array
    ImageType image(N, M, Q);

    // read image
    readImage(argv[1], image);

    // write image
    ImageType newimage(N, M, Q);

    writeImage(argv[2], newimage);

    return (1);

}
