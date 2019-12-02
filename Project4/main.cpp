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
void fft2D(int N, int M, float** real_Fuv, float** imag_Fuv, int isign);



void gausfilter(ImageType img, ImageType& filterimg, int filtersize, int N, int M, float mask[]){
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
                        applied += (value*mask[curr_mask]);
                    }
                    curr_mask++;
                }
            }

            filterimg.setPixelVal(i, j, applied);

            applied = 0;
            curr_mask = 0;
        }

    }
}

void addMotionBlur(ImageType& img, int N, int M){

}

void addGaussian(ImageType& img, int N, int M){

}

/*extern float ranf();         //ranf() is uniform in 0..1
float box_muller(float m, float s)	//normal random variate generator
{				       //mean m, standard deviation s
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        / use value from previous call
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ranf() - 1.0;
			x2 = 2.0 * ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}
*/
void convolve(ImageType img, int** filter, int N, int M, int n, int m, ImageType& out){
  int pix = 0;
  int in;
  for(int i = 0; i < N; i++){
    for(int j = 0; j< M; j++){
      for(int x = 0; x<n; x++){
        int filt_row = n - 1 - x;
        for(int y = 0;y<m; y++){
          int filt_col = m - 1 - y;

          int input_x = i  + (n/2 - filt_row);
          int input_y = j + (m/2 - filt_col);

          if(input_x >= 0 && input_x < N && input_y >=0 && input_y < M){
            img.getPixelVal(input_x,input_y,in);
            pix += in * filter[filt_row][filt_col];
            out.setPixelVal(i,j,pix);
          }
        }
      }
			pix = 0;
    }
  }

}

int* detectCosineImpulse(float** r_img, int N, int M){
  int *arr = new int[N];
  for(int x = 0; x < N; x++){
    arr[x] = 0;
  }
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M/4; j++){
      if(abs(r_img[i][j]) < .01){
        arr[i] = 0;
        break;
      }
      else{
        arr[i] = 1;
      }
    }
  }
  return arr;
}

void denoiseCosine(float** r_img, int* impulses, int N, int M){
  for(int x = 0; x < N; x++){
    cout << x << " " << impulses[x] << endl;
    for(int y = 0; y < M; y++){
      if(impulses[x] == 1){
        r_img[x][y] = 0;
      }
    }
  }
}

void extractNoise(float** r_img, float** i_img, int* impulses, int N, int M, int Q){
  float** noise_r;
  float** noise_i;
  noise_r = new float*[N];
  noise_i = new float*[N];
  for(int i = 0; i < N; i++){
    noise_r[i] = new float[M];
    noise_i[i] = new float[M];
  }

  ImageType noise(N,M,Q);
  for(int x = 0; x < N; x++){
    for(int y = 0; y < M; y++){
      if(impulses[x] == 1){
        noise_r[x][y] = r_img[x][y];
        noise_i[x][y] = i_img[x][y];
      }
      else{
        noise_r[x][y] = 0;
        noise_i[x][y] = 0;
      }
    }
  }
  fft2D(N,M,noise_r,noise_i,1);
  for(int x = 0; x < N; x++){
    for(int y = 0; y < M; y++){
    //  cout << img1[x][y] << endl;
      noise.setPixelVal(x,y,noise_r[x][y]);
    }
  }
  centerIt(noise);
  writeImage("noise.pgm", noise);
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



void magn(float** x, float** y, ImageType& newi)
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
      square = 1000 * (log2(1 + square));
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

}

#undef SWAP

int main (){
	int N,M,Q;
  bool type;


	readImageHeader("boy_noisy.pgm", N, M, Q, type);
	// allocate memory for the image array
  ImageType output(N,2*M,Q);

	ImageType image(N,M,Q);
	ImageType imag(N,M,Q);
  for(int x = 0; x < N; x++){
    for(int y = 0; y < M; y++){
      imag.setPixelVal(x,y,0);
    }
  }
	ImageType magni(N,M,Q);
	// read image
	readImage("boy_noisy.pgm", image);
  int i;
  float** img1;
  img1 = new float*[N];
  float** img2;
  img2 = new float*[N];
	float** gaus_img1;
	float** gaus_img2;
	gaus_img1 = new float*[N];
	gaus_img2 = new float*[N];

	float gausmask15 [225] = {2, 2, 3, 4, 5, 5, 6, 6, 6, 5, 5, 4, 3, 2, 2,
		                       2, 3, 4, 5, 7, 7, 8, 8, 8, 7, 7, 5, 4, 3, 2,
		                       3, 4, 6, 7, 9, 10, 10, 11, 10, 10, 9, 7, 6, 4, 3,
		                       4, 5, 7, 9, 10, 12, 13, 13, 13, 12, 10, 9, 7, 5, 4,
		                       5, 7, 9, 11, 13, 14, 15, 16, 15, 14, 13, 11, 9, 7, 5,
		                       5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5,
		                       6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6,
		                       6, 8, 11, 13, 16, 18, 19, 20, 19, 18, 16, 13, 11, 8, 6,
		                       6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6,
		                       5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5,
		                       5, 7, 9, 11, 13, 14, 15, 16, 15, 14, 13, 11, 9, 7, 5,
		                       4, 5, 7, 9, 10, 12, 13, 13, 13, 12, 10, 9, 7, 5, 4,
		                       3, 4, 6, 7, 9, 10, 10, 11, 10, 10, 9, 7, 6, 4, 3,
		                       2, 3, 4, 5, 7, 7, 8, 8, 8, 7, 7, 5, 4, 3, 2,
		                       2, 2, 3, 4, 5, 5, 6, 6, 6, 5, 5, 4, 3, 2, 2};

	int sum15 = 0;
    for(int i = 0; i<225; i++){
        sum15+=gausmask15[i];
    }
    for(int i = 0; i<225; i++){
        gausmask15[i]=gausmask15[i]/sum15;
    }
	ImageType boygaus(N,M,Q);
	ImageType gausimg(N,M,Q);
	gausfilter(image, gausimg, 15, N, M, gausmask15);
	for(int x = 0; x < N; x++){
		gaus_img1[x] = new float[M];
		gaus_img2[x] = new float[M];
		for(int y = 0; y < M; y++){
			gausimg.getPixelVal(x,y,i);
			gaus_img1[x][y] = i;
			gausimg.getPixelVal(x,y,i);
			gaus_img2[x][y] = i;
		}
	}
	writeImage("boy_denoised_gaus.pgm", gausimg);
  centerIt(image);
	centerIt(gausimg);
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

	for(int x = 0; x < N; x++){
		for(int y = 0; y < M; y++){
			gausimg.getPixelVal(x,y,i);
			gaus_img1[x][y] = i;
			gausimg.getPixelVal(x,y,i);
			gaus_img2[x][y] = i;
		}
	}


	fft2D(N,M,gaus_img1,gaus_img2, -1);
	fft2D(N,M,img1,img2,-1);

	magn(gaus_img1, gaus_img2, boygaus);
  magn(img1,img2,magni);

  writeImage("boy_noisy_magn_unshift.pgm", magni);
	writeImage("boy_noisy_gaus.pgm", boygaus);
  //Exercise 1

  int* impulse_arr = detectCosineImpulse(img1, N, M);
  extractNoise(img1, img2, impulse_arr, N, M, Q);
  denoiseCosine(img1, impulse_arr, N, M);

	fft2D(N,M,img1,img2,1);
  for(int x = 0; x < N; x++){
    for(int y = 0; y < M; y++){
      image.setPixelVal(x,y,img1[x][y]);
    }
  }
  centerIt(image);
	writeImage("boy_denoised_fft.pgm", image);


  //Exercise 2
	float** sobel_kernel;
	float** imag_kernel;
	sobel_kernel = new float*[6];
	imag_kernel = new float*[6];
	for(int x = 0; x < 6; x++){
		sobel_kernel[x] = new float[6];
		imag_kernel[x] = new float[6];
	}

	sobel_kernel[0][0] = 0, sobel_kernel[0][1] = 0, sobel_kernel[0][2] = 0, sobel_kernel[0][3] = 0, sobel_kernel[0][4] = 0, sobel_kernel[0][5] = 0,
	sobel_kernel[1][0] = 0, sobel_kernel[1][1] = 0, sobel_kernel[1][2] = 0, sobel_kernel[1][3] = 0, sobel_kernel[1][4] = 0, sobel_kernel[1][5] = 0,
	sobel_kernel[2][0] = 0, sobel_kernel[2][1] = 0, sobel_kernel[2][2] = -1, sobel_kernel[2][3] = 0, sobel_kernel[2][4] = 1, sobel_kernel[2][5] = 0,
	sobel_kernel[3][0] = 0, sobel_kernel[3][1] = 0, sobel_kernel[3][2] = -2, sobel_kernel[3][3] = 0, sobel_kernel[3][4] = 2, sobel_kernel[3][5] = 0,
	sobel_kernel[4][0] = 0, sobel_kernel[4][1] = 0, sobel_kernel[4][2] = -1, sobel_kernel[4][3] = 0, sobel_kernel[4][4] = 1, sobel_kernel[4][5] = 0,
	sobel_kernel[5][0] = 0, sobel_kernel[5][1] = 0, sobel_kernel[5][2] = 0, sobel_kernel[5][3] = 0, sobel_kernel[5][4] = 0, sobel_kernel[5][5] = 0;
	ImageType sobel(3,3,Q);
	ImageType freq_filtered(N+2, M+2, Q);
	ImageType spat_filtered(N+2, M+2, Q);
	readImageHeader("lenna.pgm", N, M, Q);
	ImageType lenna(N,M,Q);
	readImage



}
