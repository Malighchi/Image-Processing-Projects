#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
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

  double hist_arr[levels + 1];
  int test = 0;
  for(int i = 0; i < levels + 1; i++){
    hist_arr[i] = 0;
  }
  for(int x = 0; x < rows; x++){
    for (int y =0; y<cols; y++){

      ogImage.getPixelVal(x, y, val);
      if(val == 255){
        test++;
      }

      (hist_arr[val])++;
    }
  }
    cout << hist_arr[255] << endl;
    cout << test << endl;
  for(int i = 0; i < levels+1; i++){
    //  cout << hist_arr[i] << endl;
      hist_arr[i] /= (rows*cols);

    }
  //  cout << hist_arr[147] << endl;
  int k = 1;
  double new_hist[levels];
  while(k<levels+1){
    double sum = 0;
    for(int i =0; i < k; i++){
      sum = sum + (hist_arr[i] * levels);
    }

    new_hist[k-1] = floor(sum);
    k++;
    sum = 0;
    }

  for(int x = 0; x < rows; x++){
    for (int y =0; y<cols; y++){
      ogImage.getPixelVal(x, y, val);
      eqImage.setPixelVal(x, y, new_hist[val]);
      }
    }
  }
