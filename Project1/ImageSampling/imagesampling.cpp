#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"

#include "image.h"
using namespace std;

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);
void sampling(ImageType, ImageType&, int);
void upscale(ImageType, ImageType&, int);

int main(int argc, char *argv[]){
    int M, N, Q;
    bool type;
    int scale_factor = argv[3][0] - '0';
    // read image header
    readImageHeader(argv[1], N, M, Q, type);

    // allocate memory for the image array
    ImageType image(N, M, Q);

    // read image
    readImage(argv[1], image);

    // write image
    ImageType newimage(N/scale_factor, M/scale_factor, Q);
    sampling(image, newimage, scale_factor);

    ImageType resizedImage(N, M, Q);
    upscale(newimage, resizedImage, scale_factor);

    //cout<<"here"<<endl;
    writeImage(argv[2], resizedImage);

    return (1);

}

void sampling(ImageType ogImage, ImageType& newimage, int samplefactor){
    int rows = 0;
    int cols = 0;
    int levels = 0;

    ogImage.getImageInfo(rows, cols, levels);

    int x = 0;
    int y = 0;
    int val = 0;
    for(int i = 0; i < rows; i+=samplefactor){
        for (int j = 0; j < cols; j+=samplefactor){
            ogImage.getPixelVal(i, j, val);
            newimage.setPixelVal(x, y, val);
            y++;
        }
        x++;
        y=0;
    }
    //cout<<"done"<<endl;
}

void upscale(ImageType ogImage, ImageType& resizeImage, int samplefactor){
    int x = 0;
    int y = 0;
    int val = 0;

    int rows = 0;
    int cols = 0;
    int levels = 0;

    resizeImage.getImageInfo(rows, cols, levels);

    for(int i = 0; i < rows; i ++){
        for(int j = 0; j < cols; j++){
            ogImage.getPixelVal(x, y, val);
            resizeImage.setPixelVal(i, j, val);
            if((j+1) % samplefactor == 0){
                y++;
            }
        }
        if((i+1) % samplefactor == 0){
                x++;
          }
        y=0;
    }
}
