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
void quantization(ImageType, ImageType&, int);
void upscale(ImageType, ImageType&, int);

int main(int argc, char *argv[]){
    int M, N, Q;
    bool type;
    stringstream factorQ(argv[3]);
    int factor=0;
    factorQ>>factor;
    cout<<factor<<endl;
    // read image header
    readImageHeader(argv[1], N, M, Q, type);

    // allocate memory for the image array
    ImageType image(N, M, Q);

    // read image
    readImage(argv[1], image);

    // write image
    ImageType newimage(N, M, Q);
    quantization(image, newimage, factor);

    writeImage(argv[2], newimage);

    return (1);

}

void quantization(ImageType ogImage, ImageType& newimage, int quantization_factor){
    cout << quantization_factor << endl;
    int rows = 0;
    int cols = 0;
    int levels = 0;
    int new_level = 256 / quantization_factor;

    ogImage.getImageInfo(rows, cols, levels);

    int val = 0;
    for(int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            ogImage.getPixelVal(i, j, val);
            val = floor((double)val/new_level)*new_level;
            newimage.setPixelVal(i, j, val);
        }
    }
    //cout<<"done"<<endl;
}
