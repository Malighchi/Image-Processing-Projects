#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"

#include "image.h"
using namespace std;

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);
void equalization(ImageType, ImageType&);

int main(int argc, char *argv[]){
    int M, N, Q;
    bool type;
    // read image header
    readImageHeader(argv[1], N, M, Q, type);

    // allocate memory for the image array
    ImageType image(N, M, Q);
    ImageType eqImage(N, M, Q);

    // read image
    readImage(argv[1], image);




    equalization(image, eqImage);
    writeImage(argv[2], eqImage);

    return (1);

}

void equalization(ImageType ogImage, ImageType& eqImage){
  int rows = 0;
  int cols = 0;
  int levels = 0;
  int val = 0;
  ogImage.getImageInfo(rows, cols, levels);

  float hist_arr[levels];
  for(int i = 0; i < levels; i++){
    hist_arr[i] = 0;
  }
  for(int x = 0; x < rows; x++){
    for (int y =0; y<cols; y++){
      ogImage.getPixelVal(x, y, val);
      hist_arr[val]++;
    }
  }

  for(int i = 0; i < levels+1; i++){
      hist_arr[i] /= (rows*cols);
    }
  int k = 1;
  while(k<levels+1){
    float sum = 0;
    for(int i =0; i < k; i++){
      sum +=(hist_arr[i] * levels);
    }
    //cout<<sum<<endl;
    hist_arr[k-1] = sum;
    k++;
    }

  for(int x = 0; x < rows; x++){
    for (int y =0; y<cols; y++){
      ogImage.getPixelVal(x, y, val);
      eqImage.setPixelVal(x, y, hist_arr[val]);
      }
    }
  }
