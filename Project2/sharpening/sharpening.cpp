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

void gradientMagn(ImageType, ImageType, ImageType&);
void filter(ImageType img, ImageType& filterimg, int filtersize, int N, int M, int mask[]);

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
      newPix = (double)((x-min)*255)/(double)(max-min);
      pix = (int)newPix;
      newImage.setPixelVal(k,l,pix);
    }
  }

}

int main(int argc, char *argv[]){
    int M, N, Q;
    int mask_size = 3;
    bool type;
    int v_prewitt[9] = {-1,0,1,-1,0,1,-1,0,1};
    int h_prewitt[9] = {-1,-1,-1,0,0,0,1,1,1};
    int v_sobel[9] = {-1,-2,-1,0,0,0,1,2,1};
    int h_sobel[9] = {-1,0,1,-2,0,2,-1,0,2};
    int laplacian[9] = {0,1,0,1,-4,1,0,1,0};



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
    ImageType gradS(N,M,Q);
    ImageType gradP(N,M,Q);
    ImageType normalLap(N,M,Q);
    ImageType normalS(N,M,Q);
    ImageType normalHS(N,M,Q);
    ImageType normalVS(N,M,Q);
    ImageType normalVP(N,M,Q);
    ImageType normalHP(N,M,Q);
    ImageType normalP(N,M,Q);

    filter(image, vsimg,mask_size, N, M, v_sobel);
    filter(image, hsimg,mask_size,N,M,h_sobel);
    filter(image, vpimg,mask_size,N,M,v_prewitt);
    filter(image, hpimg,mask_size,N,M, h_prewitt);
    filter(image, lapimg,mask_size,N,M, laplacian);

    gradientMagn(hsimg,vsimg,gradS);
    gradientMagn(hpimg,vpimg,gradP);

    normalizeImage(vsimg, normalVS);
    normalizeImage(hsimg, normalHS);
    normalizeImage(lapimg, normalLap);
    normalizeImage(gradS, normalS);
    normalizeImage(vpimg, normalVP);
    normalizeImage(hpimg, normalHP);
    normalizeImage(gradP, normalP);

    //gradientMagn(hsimg,vsimg,gradS);
    //gradientMagn(hpimg,vpimg,gradP);
    writeImage(argv[2], normalVS);
    writeImage(argv[3], normalHS);
    writeImage(argv[4], normalVP);
    writeImage(argv[5], normalHP);
    writeImage(argv[6], normalLap);
    writeImage(argv[7], normalS);
    writeImage(argv[8], normalP);

    return (1);

}


void gradientMagn(ImageType x, ImageType y, ImageType& newi)
{
  int xpix;
  int ypix;
  int square;
  int N, M, Q;
  x.getImageInfo(N,M,Q);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      x.getPixelVal(i,j,xpix);
      y.getPixelVal(i,j,ypix);
      square = sqrt((xpix * xpix) + (ypix * ypix));
      newi.setPixelVal(i,j,square);
    }
  }
}
void filter(ImageType img, ImageType& filterimg, int filtersize, int N, int M, int mask[]){
    float applied = 0;
    int curr_mask = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            for(int s = (-filtersize/2); s <= (filtersize/2); s++){
                for(int t = (-filtersize/2); t <= (filtersize/2); t++){
                    int value = 0;
                    if(((i+s)<0 || (i+s)>=N) || ((j+t)<0 || (j+t)>=M)){
                        applied+=0;
                    }
                    else{
                        img.getPixelVal(i+s, j+t, value);
                        //cout<<value<<" "<<mask[curr_mask]<<endl;
                        applied += (value*mask[curr_mask]);

                    }
                    curr_mask++;
                }
            }
            //cout<<applied<<endl;
            filterimg.setPixelVal(i, j, applied);

            applied = 0;
            curr_mask = 0;
        }

    }
}
