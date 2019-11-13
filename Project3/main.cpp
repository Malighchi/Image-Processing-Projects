
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include "stdlib.h"

#include "image.h"
#include <iostream>
#include <math.h>
using namespace std;

double pi = 3.14159265358979323846;
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);
void centerIt(ImageType&);

float* generateCos(float u, float N){
  int n = N;
  float* wave = new float[2*n];
  float start = -(N/2);
  int count = 0;
  for(float i = 0; i < N; i++, count++){
    //cout << i << endl;
    wave[count] = 0;
    count++;
    //cout << cos(2*pi*u*(i/N)) << endl;
    wave[count] = (cos(2*pi*u*(i/N))) * pow(-1, i);
  //  cout << wave[count] << endl;
  }
  return wave;
}

float* phase(float* wave, int nn){
  float pix;
  float* ph = new float[nn/2];
  for(int i = 1; i < nn; i+=2){
    pix = atan(wave[i+1] / wave[i]);
    ph[i] = pix;
  }
  return ph;
}

void centerIt(ImageType& img){
  int N, M, Q;
  int pix;
  img.getImageInfo(N,M,Q);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){

      img.getPixelVal(i,j,pix);

      pix = pix * pow(-1, i + j);

      img.setPixelVal(i,j,pix);
    }
    }
  }



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
      //xpix = xpix * pow(-1, i+j);
      //.getPixelVal(i,j,xpix);
      ypix = y[i][j];
      //ypix = ypix * pow(-1, i+j);
      //.getPixelVal(i,j,ypix);
      square = sqrt((xpix * xpix) + (ypix * ypix));
    //  square = square * pow(-1,i + j);
      square = 100 * (log2(1 + square));
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
	//Problem 1.a
  float data[9] = {0,2,0,3,0,4,0,4,0};
  fft(data,4,-1);
  for(int i = 1; i < 9; i++)
    data[i] /= 4;
  for(int i = 1; i < 9; i+=2){

     cout << data[i] <<" + "<<data[i+1]<<"j"<<endl;
	 cout << "Magnitude: "<<sqrt(((data[i]*data[i])+(data[i+1]*data[i+1])))<<endl;
  }

  fft(data,4,1);
  for(int i = 1; i < 9; i+=2){
     cout << data[i] << endl;
  }

  //Problem 1.b

  ofstream cos,cosr,cosi,cosp,cosm;
  cos.open("cos.txt",ios::out);
  cosr.open("cosr.txt",ios::out);
  cosi.open("cosi.txt",ios::out);
  cosp.open("cosp.txt",ios::out);
  cosm.open("cosm.txt",ios::out);
  float* cosine = new float[257];
  cosine = generateCos(8,128);

  int out_count = -64;
  for(int i = 1; i < 257; i = i + 2){
    if(abs(cosine[i]) < .00001){
      cosine[i] = 0;
    }
    cos <<"("<<out_count<<", "<< cosine[i] <<")"<< endl;
    out_count++;
  //  cosine[i] = cosine[i] * pow(-1,i);
  }

  fft(cosine,128, -1);

  out_count = -64;
  for(int i = 0; i < 257; i++){
    cosine[i] /= 128;
  }

  for(int i = 2 ; i < 257; i = i + 2 ){
    if(abs(cosine[i]) < .00001){
      cosine[i] = 0;
    }
    cosi <<"("<<out_count<<", "<< cosine[i] <<")"<< endl;
    out_count++;
  }
  out_count = -64;
  for(int i = 1; i < 257; i = i + 2){
    if(abs(cosine[i]) < .00001){
      cosine[i] = 0;
    }
    cosr <<"("<<out_count<<", "<< cosine[i] <<")"<< endl;
    out_count++;
  }

  float* ph = new float[128];
  ph = phase(cosine, 256);
  out_count = -64;
  for(int j = 0; j < 128; j++){
    if(abs(ph[j]) < .00001){
      ph[j] = 0;
    }
    cosp <<"("<<out_count<<", "<< ph[j] <<")"<< endl;
    out_count++;
  }
  out_count = -64;
  for(int j = 1; j < 256; j += 2){
    cosm <<"("<<out_count<<", "<< sqrt(((cosine[j]*cosine[j])+(cosine[j+1]*cosine[j+1]))) <<")"<< endl;
    out_count++;
  }

  //problem 1.c
  float rectangle[257];
  for(int i = 0; i < 257; i++){
    rectangle[i] = 0;
  }
  for(int i = 1; i < 257; i+=2){
    if(i > 64 && i < 193){
      rectangle[i] = 1.0 * pow(-1,(i/2));
    }
    else{
      rectangle[i] = 0.0;
    }
  }

  ofstream rect,rectr,recti,rectp,rectm;
  rect.open("rect.txt",ios::out);
  rectr.open("rectr.txt",ios::out);
  recti.open("recti.txt",ios::out);
  rectp.open("rectp.txt",ios::out);
  rectm.open("rectm.txt",ios::out);
  out_count = -64;
  for(int i = 1; i < 257; i = i + 2){
    if(abs(rectangle[i]) < .00001){
      rectangle[i] = 0;
    }
    rect <<"("<<out_count<<", "<<rectangle[i] <<")"<< endl;
    out_count++;
  //  cosine[i] = cosine[i] * pow(-1,i);

  }
  fft(rectangle,128, -1);

  for(int i = 0; i < 257; i++){
    rectangle[i] /= 128;
  }
  out_count = -64;
  for(int i = 2 ; i < 257; i = i + 2 ){
    if(abs(rectangle[i]) < .00001){
      rectangle[i] = 0;
    }
    recti <<"("<<out_count<<", "<<rectangle[i] <<")"<< endl;
    out_count++;
  }
  out_count = -64;
  for(int i = 1; i < 257; i = i + 2){
    if(abs(rectangle[i]) < .00001){
      rectangle[i] = 0;
    }
    rectr <<"("<<out_count<<", "<<rectangle[i] <<")"<< endl;
    out_count++;
  }

  float* phr = new float[128];
  phr = phase(rectangle, 256);
  out_count = -64;
  for(int j = 0; j < 128; j++){
    if(abs(phr[j]) < .00001){
      phr[j] = 0;
    }
    rectp <<"("<<out_count<<", "<< phr[j] <<")"<< endl;
    out_count++;
  }

  out_count = -64;
  for(int j = 1; j < 256; j += 2){
    rectm <<"("<<out_count<<", "<< sqrt(((rectangle[j]*rectangle[j])+(rectangle[j+1]*rectangle[j+1]))) <<")"<< endl;
    out_count++;
  }



  //problem 2
	readImageHeader("square128.pgm", N, M, Q, type);
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
	readImage("square128.pgm", image);
  int i;
  float** img1;
  img1 = new float*[N];
  float** img2;
  img2 = new float*[N];
//  centerIt(image);
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
  writeImage("square128_magn_unshift.pgm", magn);

	//writeImage("square32_n.pgm", magn);
	fft2D(N,M,img1,img2,1);
  for(int x = 0; x < N; x++){
    for(int y = 0; y < M; y++){
    //  cout << img1[x][y] << endl;
      image.setPixelVal(x,y,img1[x][y]);
    }
  }
//	writeImage("square64_m.pgm", image);
  return 1;


}
