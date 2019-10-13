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
int correlationSum(int, int, int**);
void correlation(ImageType, ImageType&, ImageType);
void normalizeImage(ImageType, ImageType&);
int** convertToMask(ImageType, int, int);
int** subimage(ImageType, int, int, int, int);

int main(int argc, char *argv[]){
    int M, N, Q;
    int M_mask, N_mask, Q_mask;
    bool type;

    // read image header
    readImageHeader(argv[1], N, M, Q, type);
    readImageHeader(argv[2], N_mask, M_mask, Q_mask, type);
    // allocate memory for the image array
    ImageType image(N, M, Q);
    ImageType mask(N_mask, M_mask, Q_mask);
    // read image
    readImage(argv[1], image);
    readImage(argv[2], mask);

    // write image
    ImageType newimage(N, M, Q);
    ImageType normalizedImage(N, M, Q);
    correlation(image, newimage, mask);
    normalizeImage(newimage, normalizedImage);
    writeImage(argv[3], normalizedImage);

    return (1);

}

int correlationSum(int** pixels, int n, int m, int** mask){
  int sum = 0, temp = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      temp = mask[i][j] * pixels[i][j];
      sum += temp;
    }
  }
  return sum;
}

void correlation(ImageType oldImg, ImageType& newImg, ImageType mask){
  int M, N, Q;
  int M_mask, N_mask, Q_mask;
  oldImg.getImageInfo(N,M,Q);
  mask.getImageInfo(N_mask,M_mask,Q_mask);


  int** filter = convertToMask(mask, N_mask, M_mask);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      int** sub = subimage(oldImg, i, j, N_mask, M_mask);
      int sum = correlationSum(sub, N_mask, M_mask, filter);
      newImg.setPixelVal(i,j,sum);
    }
  }
}

int** subimage(ImageType img, int x, int y, int n_mask, int m_mask){
  int** sub = 0;
  sub = new int*[n_mask];
  int M = 0;
  int N = 0;
  int Q = 255;
  int x_count = 0;
  int y_count = 0;
  int pixel = 0;
  img.getImageInfo(N,M,Q);
  //cout << n_mask << endl;
  //cout << m_mask << endl;
  for(int k = 0; k < n_mask; k++){
    sub[k] = new int[m_mask];
  }
  for(int i = -n_mask/2; i<n_mask/2; i++){
    for(int j = -m_mask/2; j< m_mask/2; j++){
  //    cout << y_count << " " << x_count << endl;
      if((x+i>=0 && x+i<N) && (y+j>=0 && y+j<M)){
        img.getPixelVal(x+i, y+j, pixel);
        sub[x_count][y_count] = pixel;
      }
      else{
        sub[x_count][y_count] = 0;
      }
      y_count++;
    }
    y_count = 0;
    x_count++;
}

  return sub;
}

void normalizeImage(ImageType ogImage, ImageType& newImage){
  int M = 0;
  int N = 0;
  int Q = 0;
  int x = 0;
  int counter = 0;
  double newPix = 0;
  int pix = 0;
  ogImage.getImageInfo(N, M, Q);
  int max = 0;
  int min = 1000;
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      ogImage.getPixelVal(i,j,x);
      max = std::max(max, x);
      min = std::min(min, x);
    }
  }
  cout << max << " " << min << endl;
  for(int k = 0; k < N; k++){
    for(int l = 0; l < M; l++){
      ogImage.getPixelVal(k,l,x);
      newPix = (double)(x*255)/(double)max;
      pix = (int)newPix;
      newImage.setPixelVal(k,l,pix);
    }
  }
//  cout << counter << endl;
}

int** convertToMask(ImageType mask, int N, int M){
  int pixel;
  int** new_mask = 0;
  new_mask = new int*[N];

  for(int i = 0; i < N; i++){
    new_mask[i] = new int[M];
    for(int j = 0; j < M; j++){
      mask.getPixelVal(i,j,pixel);
      new_mask[i][j] = pixel;
    }
  }
  return new_mask;
}
