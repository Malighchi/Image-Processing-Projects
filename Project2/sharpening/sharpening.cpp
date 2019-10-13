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
int applyMask(int, int, int**);
void filter(ImageType, ImageType&, int**, int);
void normalizeImage(ImageType, ImageType&);
int** subimage(ImageType, int, int, int);

int main(int argc, char *argv[]){
    int M, N, Q;
    int mask_size = 3;
    bool type;
    int** h_prewitt = 0;
    int** v_prewitt = 0;
    int** h_sobel = 0;
    int** v_sobel = 0;
    int** laplacian = 0;
    h_prewitt = new int*[mask_size];
    v_prewitt = new int*[mask_size];
    h_sobel = new int*[mask_size];
    v_sobel = new int*[mask_size];
    laplacian = new int*[mask_size];
    for(int x = 0; x < mask_size; x++){
      h_prewitt[x] = new int[mask_size];
      v_prewitt[x] = new int[mask_size];
      h_sobel[x] = new int[mask_size];
      v_sobel[x] = new int[mask_size];
      laplacian[x] = new int[mask_size];
    }

    h_prewitt[0][0] = -1, h_prewitt[0][1] = 0, h_prewitt[0][2] = 1;
    h_prewitt[1][0] = -1, h_prewitt[1][1] = 0, h_prewitt[1][2] = 1;
    h_prewitt[2][0] = -1, h_prewitt[2][1] = 0, h_prewitt[2][2] = 1;

    v_prewitt[0][0] = -1, v_prewitt[0][1] = -1, v_prewitt[0][2] = -1;
    v_prewitt[1][0] = 0, v_prewitt[1][1] = 0, v_prewitt[1][2] = 0;
    v_prewitt[2][0] = 1, v_prewitt[2][1] = 1, v_prewitt[2][2] = 1;

    h_sobel[0][0] = -1, h_sobel[0][1] = 0, h_sobel[0][2] = 1;
    h_sobel[1][0] = -2, h_sobel[1][1] = 0, h_sobel[1][2] = 2;
    h_sobel[2][0] = -1, h_sobel[2][1] = 0, h_sobel[2][2] = 1;

    v_sobel[0][0] = -1, v_sobel[0][1] = -2, v_sobel[0][2] = -1;
    v_sobel[1][0] = 0, v_sobel[1][1] = 0, v_sobel[1][2] = 0;
    v_sobel[2][0] = 1, v_sobel[2][1] = 2, v_sobel[2][2] = 1;

    laplacian[0][0] = 0, laplacian[0][1] = 1, laplacian[0][2] = 0;
    laplacian[1][0] = 1, laplacian[1][1] = -4, laplacian[1][2] = 1;
    laplacian[2][0] = 0, laplacian[2][1] = 1, laplacian[2][2] = 0;
    // read image header
    readImageHeader(argv[1], N, M, Q, type);
    // allocate memory for the image array
    ImageType image(N, M, Q);
    // read image
    readImage(argv[1], image);

    // write image
    ImageType hsimg(N, M, Q);
    ImageType vsimg(N,M,Q);
    ImageType hpimg(N,M,Q);
    ImageType vpimg(N,M,Q);
    ImageType lapimg(N,M,Q);

    filter(image, hsimg, h_sobel, mask_size);
    filter(image, vsimg, v_sobel, mask_size);
    filter(image, hpimg, h_prewitt, mask_size);
    filter(image, vpimg, v_prewitt, mask_size);
    filter(image, lapimg, laplacian, mask_size);
    writeImage(argv[2], hsimg);
    writeImage(argv[3], vsimg);
    writeImage(argv[4], hpimg);
    writeImage(argv[5], vpimg);
    writeImage(argv[6], lapimg);

    return (1);

}

int applyMask(int** pixels, int mask_size, int** mask){
  int sum = 0, temp = 0;
  for(int i = 0; i < mask_size; i++){
    for(int j = 0; j < mask_size; j++){
      temp = mask[i][j] * pixels[i][j];
      sum += temp;
    }
  }
  return sum;
}

void filter(ImageType oldImg, ImageType& newImg, int** mask, int mask_size){
  int M, N, Q;
  oldImg.getImageInfo(N,M,Q);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      int** sub = subimage(oldImg, i, j, mask_size);
      int sum = applyMask(sub, mask_size, mask);
      newImg.setPixelVal(i,j,sum);
    }
  }
}

int** subimage(ImageType img, int x, int y, int mask_size){
  int** sub = 0;
  sub = new int*[mask_size];
  int M = 0;
  int N = 0;
  int Q = 255;
  int x_count = 0;
  int y_count = 0;
  int pixel = 0;
  img.getImageInfo(N,M,Q);
  //cout << n_mask << endl;
  //cout << m_mask << endl;
  for(int k = 0; k < mask_size; k++){
    sub[k] = new int[mask_size];
  }
  for(int i = -mask_size/2; i<mask_size/2; i++){
    for(int j = -mask_size/2; j< mask_size/2; j++){
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
