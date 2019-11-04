
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include "stdlib.h"

#include "image.h"
#include <iostream>
#include <math.h>
using namespace std;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);

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

void fft(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}




void fft2D(int N, int M, ImageType& real_Fuv, ImageType& imag_Fuv, int isign){
	std::cout <<N << M <<std::endl;
	int x;
	for(int i = 0; i < N; i++){
		float arr[M*2+1];
		for(int count = 0; count < (M*2+1); count++){
			arr[count] = 0;
		}
		int arr_cntr = 1;
		for(int j = 0; j < M; j++){
			real_Fuv.getPixelVal(i,j,x);
			arr[arr_cntr] = x;
			arr_cntr++;
			imag_Fuv.getPixelVal(i,j,x);
			arr[arr_cntr] = x;
			arr_cntr++;
		}
		fft(arr, M, isign);
		arr_cntr=1;
		int forward = 1;
		if(isign == -1){
			forward = M;
		}
		for(int j = 0; j < M; j++){
			real_Fuv.setPixelVal(i,j,(arr[arr_cntr]/forward));
			//real_Fuv[i][j] = arr[arr_cntr]/forward;
			arr_cntr++;
			imag_Fuv.setPixelVal(i,j,(arr[arr_cntr]/forward));
			//imag_Fuv[i][j] = arr[arr_cntr]/forward;
			arr_cntr++;
		}
	}
	for(int i = 0; i < M; i++){
		float arr[N*2+1];
		for(int count = 0; count < (N*2+1); count++){
			arr[count] = 0;
		}
		int arr_cntr = 1;
		for(int j = 0; j < N; j++){
			real_Fuv.getPixelVal(j,i,x);
			arr[arr_cntr] = x;
			arr_cntr++;
			imag_Fuv.getPixelVal(j,i,x);
			arr[arr_cntr] = x;
			arr_cntr++;
		}
		fft(arr, N, isign);
		arr_cntr=1;
		int forward = 1;
		if(isign == -1){
			forward = N;
		}
		for(int j = 0; j < N; j++){
			real_Fuv.setPixelVal(i,j,(arr[arr_cntr]/forward));
			//real_Fuv[i][j] = arr[arr_cntr]/forward;
			arr_cntr++;
			imag_Fuv.setPixelVal(i,j,(arr[arr_cntr]/forward));
			//imag_Fuv[i][j] = arr[arr_cntr]/forward;
			arr_cntr++;
		}
	}
}

#undef SWAP

int main (){
	int N,M,Q;
	bool type;
	/*
  float data[9] = {0,2,0,3,0,4,0,4,0};
  fft(data,4,-1);
  for(int i = 1; i < 9; i++)
    data[i] /= 4;
  for(int i = 1; i < 9; i+=2){

     std::cout << data[i] <<" + "<<data[i+1]<<"j"<<std::endl;
	 std::cout << "Magnitude: "<<sqrt(((data[i]*data[i])+(data[i+1]*data[i+1])))<<std::endl;
  }
  std::cout << std::endl;

  fft(data,4,1);
  for(int i = 1; i < 9; i+=2){
     std::cout << data[i] << std::endl;
  }*/

	readImageHeader("lenna.pgm", N, M, Q, type);
	// allocate memory for the image array
	ImageType image(N, M, Q);
	ImageType imag(N,M,Q);
	ImageType magn(N,M,Q);
	// read image
	readImage("lenna.pgm", image);
	fft2D(N,M,image,imag,-1);
	gradientMagn(image,imag,magn);
	writeImage("lenna1.pgm", magn);
	fft2D(N,M,image,imag,1);
	writeImage("lenna2.pgm", image);
  return 1;


}
