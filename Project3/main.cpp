
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


void gradientMagn(float** x, float** y, ImageType& newi)
{
  float xpix;
  float ypix;
  float square;
  int N, M, Q;
  newi.getImageInfo(N,M,Q);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      xpix = x[i][j];
      //.getPixelVal(i,j,xpix);
      ypix = y[i][j];
      //.getPixelVal(i,j,ypix);
      square = sqrt((xpix * xpix) + (ypix * ypix));
      square = square * pow(-1,i + j);
      square = (log2(1 + square));
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




void fft2D(int N, int M, float** real_Fuv, float** imag_Fuv, int isign){
//  if(isign == -1){
	int x;
  //float real[N][M];
  //float imag[N][M];
  //for(int o = 0; o < N; o++){
    //for(int u = 0; u < M; u++){
    //  real_Fuv.getPixelVal(o,u,x);
    //  cout << x << " ";
    //  real[o][u] = x;
    //  imag_Fuv.getPixelVal(o,u,x);
    //  real[o][u] = x;
  //    cout << x << " ";
    //}
  //cout << real[1][1] << " " << imag[1][1];
	for(int i = 0; i < N; i++){
		float arr[M*2+1];
		for(int count = 0; count < (M*2+1); count++){
			arr[count] = 0;
		}
		int arr_cntr = 1;
		for(int j = 0; j < M; j++){
		//	real_Fuv.getPixelVal(i,j,x);
			arr[arr_cntr] = real_Fuv[i][j];
      //cout << x << " ";
			arr_cntr++;
		//	imag_Fuv.getPixelVal(i,j,x);
			arr[arr_cntr] = imag_Fuv[i][j];
    //  cout << x << " ";
			arr_cntr++;
		}
    //cout << endl;
		fft(arr, M, isign);

		arr_cntr=1;
		float forward = 1;
		if(isign == -1){
			forward = M;
		}
		for(int j = 0; j < M; j++){
			//real_Fuv.setPixelVal(i,j,(arr[arr_cntr]/forward));
			real_Fuv[i][j] = arr[arr_cntr]/forward;
			arr_cntr++;
			//imag_Fuv.setPixelVal(i,j,(arr[arr_cntr]/forward));
			imag_Fuv[i][j] = arr[arr_cntr]/forward;
			arr_cntr++;
		}
//cout << endl << endl << endl;
	}
	for(int i = 0; i < M; i++){
		float arr[N*2+1];
		for(int count = 0; count < (N*2+1); count++){
			arr[count] = 0;
		}
		int arr_cntr = 1;
		for(int j = 0; j < N; j++){
		//	real_Fuv.getPixelVal(j,i,x);
			arr[arr_cntr] = real_Fuv[j][i];
			arr_cntr++;
		//	imag_Fuv.getPixelVal(j,i,x);
			arr[arr_cntr] = imag_Fuv[j][i];
			arr_cntr++;
		}
		fft(arr, N, isign);
		arr_cntr=1;
		float forward = 1;
		if(isign == -1){
			forward = N;
		}
		for(int j = 0; j < N; j++){
			//real_Fuv.setPixelVal(j,i,(arr[arr_cntr]/forward));
			real_Fuv[j][i] = arr[arr_cntr]/forward;
			arr_cntr++;
			//imag_Fuv.setPixelVal(j,i,(arr[arr_cntr]/forward));
			imag_Fuv[j][i] = arr[arr_cntr]/forward;
			arr_cntr++;
		}
  }
  //for(int o = 0; o < N; o++){
    //for(int u = 0; u < M; u++){
    //  real_Fuv.setPixelVal(o,u,real[o][u]);
      //imag_Fuv.setPixelVal(o,u,imag[o][u]);
    //}
  //}
	//}
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
	ImageType image(N,M,Q);
	ImageType imag(N,M,Q);
  for(int x = 0; x < N; x++){
    for(int y = 0; y < M; y++){
      imag.setPixelVal(x,y,0);
    }
  }
	ImageType magn(N,M,Q);
	// read image
	readImage("lenna.pgm", image);
  int i;
  float** img1;
  img1 = new float*[N];
  float** img2;
  img2 = new float*[N];
  for(int x = 0; x < N; x++){
    img1[x] = new float[M];
    img2[x] = new float[M];
    for(int y = 0; y < M; y++){
      image.getPixelVal(x,y,i);
      img1[x][y] = i;
      imag.getPixelVal(x,y,i);
      img2[x][y] = i;
    }
  }

	fft2D(N,M,img1,img2,-1);
	gradientMagn(img1,img2,magn);
	writeImage("lenna1.pgm", magn);
	fft2D(N,M,img1,img2,1);
  for(int x = 0; x < N; x++){
    for(int y = 0; y < M; y++){
      image.setPixelVal(x,y,img1[x][y]);
    }
  }
	writeImage("lenna5.pgm", image);
  return 1;


}
