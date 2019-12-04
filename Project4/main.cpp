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

void normalizeImage(float** image, int N, int M){

  float max = -10000000000;
  float min = 10000000000;
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      max = std::max(max, image[i][j]);
      min = std::min(min, image[i][j]);
    }
  }

  for(int k = 0; k < N; k++){
    for(int l = 0; l < M; l++){
      image[k][l] = ((image[k][l]-min)*255)/(max-min);
    }
  }
}

float** createSobel(int N, int M, float** freq){

	for(int i = 0; i < N; i++){

		for(int j = 0; j < M; j++){
			freq[i][j] = 0;
		}
	}

	for(int i = 0; i < N-2; i++){
		for(int j = 0; j < M-2; j++){
			freq[i][j] = -1;
			freq[i+1][j] = -2;
			freq[i+2][j] = -1;
			freq[i][j+1] = 0;
			freq[i][j+2] = 1;
			freq[i+1][j+1] = 0;
			freq[i+1][j+2] = 2;
			freq[i+2][j+1] = 0;
			freq[i+2][j+2] = 1;
			if(i!=0 && j!= 0 && (freq[i][j] == -freq[N-i][M-j])){
					return freq;
				}
				else{
							freq[i][j] = 0;
							freq[i+1][j] = 0;
							freq[i+2][j] = 0;
							freq[i][j+1] = 0;
							freq[i][j+2] = 0;
							freq[i+1][j+1] = 0;
							freq[i+1][j+2] = 0;
							freq[i+2][j+1] = 0;
							freq[i+2][j+2] = 0;
						}
					}

				}



	return freq;
}

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

void addMotionBlur(float** h, int N, int M, float a, float b, float T){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			h[i][j] = (T/(pi*((i-floor(N/2))*a + (j-floor(M/2))*b)))*sin(pi*((i-floor(N/2))*a+(j-floor(M/2))*b))*exp((-1)*pi*((i-floor(N/2))*a+(j-floor(M/2))*b));
			/*if(abs(h[i][j]) < .0001){
				h[i][j] = 1;
			}*/
			h[i][j] = 100000 * (log2(1 + h[i][j]));
		}

	}
}

float ranf(){        //ranf() is uniform in 0..1
	float r = static_cast <float> (rand()) / static_cast <float> (1);
	return r;
}

float box_muller(float m, float s)	//normal random variate generator
{				       //mean m, standard deviation s
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)
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

void addGaussian(float** f, int N, int M, float m, float s){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			f[i][j] += box_muller(m,s);
		}
	}
}

void inverseFiltering(float** f, float** h, float** g, int N, int M){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			f[i][j] = (1/h[i][j])*g[i][j];
		}
	}
}

void wienerFiltering(float** f, float** h, float** g, int N, int M, float K){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			f[i][j] = ((abs(h[i][j]) * abs(h[i][j]))/(abs(h[i][j])*abs(h[i][j])+K)) * g[i][j]/h[i][j];
		}
	}
}




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

void centerFloat(float** img, int N, int M){
	  for(int i = 0; i < N; i++){
	    for(int j = 0; j < M; j++){


	      img[i][j] = img[i][j] * pow(-1, i + j);

	    }
	 }
}

void magn(float** x, float** y, ImageType& newi, int N, int M)
{
  float xpix;
  float ypix;
  float square;

  for(int i = 0; i < N; i++){

    for(int j = 0; j < M; j++){
      xpix = x[i][j];

      ypix = y[i][j];

      square = sqrt((xpix * xpix) + (ypix * ypix));

      square = 1000 * (log2(1 + square));
      newi.setPixelVal(i,j,square);
    }
  }
}



void float_magn(float** x, float** y, float** newi, int N, int M)
{
  float xpix;
  float ypix;
  float square;

  for(int i = 0; i < N; i++){

    for(int j = 0; j < M; j++){
      xpix = x[i][j];

      ypix = y[i][j];

      square = sqrt((xpix * xpix) + (ypix * ypix));

      square = 1000 * (log2(1 + square));
      newi[i][j] = square;
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

			arr[arr_cntr] = real_Fuv[i][j];

			arr_cntr++;

			arr[arr_cntr] = imag_Fuv[i][j];

			arr_cntr++;
		}

		fft(arr, M, isign);

		arr_cntr=1;
		float forward = 1;
		if(isign == -1){
			forward = M;
		}
		for(int j = 0; j < M; j++){

			real_Fuv[i][j] = arr[arr_cntr]/forward;
			arr_cntr++;

			imag_Fuv[i][j] = arr[arr_cntr]/forward;
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

			arr[arr_cntr] = real_Fuv[j][i];
			arr_cntr++;

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

			real_Fuv[j][i] = arr[arr_cntr]/forward;
			arr_cntr++;

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
	ImageType imag(2*N,2*M,Q);
	ImageType image(N,M,Q);

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

	magn(gaus_img1, gaus_img2, boygaus, N, M);
  magn(img1,img2,magni, N, M);

  writeImage("boy_noisy_magn.pgm", magni);
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


	readImageHeader("lenna.pgm", N, M, Q, type);
	//ImageType image1(2*N,2*M,Q);
	ImageType lenna(N,M,Q);
	readImage("lenna.pgm", lenna);
	ImageType freq_filtered(2*N, 2*M, Q);
	ImageType spat_filtered(N, M, Q);
	ImageType sobel(2*N,2*M,Q);

	float** sobel_kernel;
	float** imag_kernel;
	sobel_kernel = new float*[2*N];
	imag_kernel = new float*[2*N];
	for(int x = 0; x < 2*N; x++){
		sobel_kernel[x] = new float[2*M];
		imag_kernel[x] = new float[2*M];
	}


  float** real_lenna;
	float** imag_lenna;
	real_lenna = new float*[2*N];
	imag_lenna = new float*[2*N];
	for(int x = 0; x < 2*N; x++){
		real_lenna[x] = new float[2*M];
		imag_lenna[x] = new float[2*M];
	}


  for(int i = 0; i < 2*N; i++){
    for(int j = 0; j < 2*M; j++){
      real_lenna[i][j] = 0;
      imag_lenna[i][j] = 0;
			imag_kernel[i][j] = 0;
			sobel_kernel[i][j] = 0;
    }
  }

  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++){
      int value = 0;
      lenna.getPixelVal(i, j, value);
      real_lenna[i][j] = value;
    }
  }


	sobel_kernel = createSobel(2*N,2*M, sobel_kernel);

	centerFloat(real_lenna,2*N,2*M);
	centerFloat(sobel_kernel, 2*N, 2*M);
	//normalizeImage(sobel_kernel,2*N,2*M);
	fft2D(2*N,2*M,sobel_kernel, imag_kernel, -1);
	centerFloat(imag_kernel, 2*N, 2*M);
	for(int i = 0; i < 2*M; i++){
		for(int j = 0; j < 2*N; j++){
			sobel_kernel[i][j] = 0;
			imag_kernel[i][j] = 1000 * (log2(1 + imag_kernel[i][j]));
		}
	}
	fft2D(2*N,2*M,real_lenna,imag_lenna,-1);
	//normalizeImage(imag_kernel,2*N,2*M);
	magn(sobel_kernel, imag_kernel, sobel, 2*N,2*M);

	for(int i = 0; i < 2*N; i++){
		for(int j = 0; j < 2*M; j++){
			float temp = real_lenna[i][j];
			real_lenna[i][j]=sobel_kernel[i][j]*real_lenna[i][j]-imag_kernel[i][j]*imag_lenna[i][j];
			imag_lenna[i][j]=temp*imag_kernel[i][j]-imag_lenna[i][j]*sobel_kernel[i][j];
		}
	}


	fft2D(2*N,2*M, real_lenna, imag_lenna, 1);
	//normalizeImage(real_lenna,2*N,2*M);
	for(int i = 0; i < 2*N; i++){
		for(int j = 0; j < 2*M; j++){
			freq_filtered.setPixelVal(i,j, real_lenna[i][j]);
		}
	}

	centerIt(freq_filtered);

	writeImage("sobel_freq.pgm", sobel);
	writeImage("sobel_freq_lenna.pgm", freq_filtered);


	//Exercise 3
	ImageType filtered(2*N,2*M,Q);
	float** h_real;
	float** h;
	float** real_ex3;
	float** imag_ex3;
	h = new float*[2*N];
	h_real = new float*[2*N];
	real_ex3 = new float*[2*N];
	imag_ex3 = new float*[2*N];

	for(int x = 0; x < 2*N; x++){
		h_real[x] = new float[2*N];
		h[x] = new float[2*M];
		real_ex3[x] = new float[2*M];
		imag_ex3[x] = new float[2*M];
	}
	for(int i = 0; i < 2*N; i++){
		for(int j = 0; j < 2*M; j++){
			h_real[i][j] = 0;
			h[i][j] = 0;
		}
	}
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			int value = 0;
			lenna.getPixelVal(i, j, value);
			real_ex3[i][j] = value;
		}
	}
	centerFloat(real_ex3, 2*N, 2*M);
	fft2D(2*N,2*M, real_ex3, imag_ex3, -1);
	addMotionBlur(h,2*N, 2*M, .1, .1, 1);

	//	real_lenna[i][j] = (real_lenna[i][j] * sobel_kernel[i][j]) - (imag_lenna[i][j] * imag_kernel[i][j]) ;
	//	imag_lenna[i][j] =(real_lenna[i][j] * imag_kernel[i][j]) + (imag_lenna[i][j] * sobel_kernel[i][j]) ;
	for(int i = 0; i < 2*N; i++){
		for(int j = 0; j < 2*M; j++){
			float temp = real_ex3[i][j];
			//real_ex3[i][j]=h_real[i][j]*real_ex3[i][j]-h[i][j]*imag_ex3[i][j];
			//imag_ex3[i][j]=temp*h[i][j]-imag_ex3[i][j]*h_real[i][j];
			real_ex3[i][j] = (h[i][j]*real_ex3[i][j]);
			imag_ex3[i][j] = imag_ex3[i][j]* h[i][j];
		}
	}

fft2D(2*N,2*M, real_ex3, imag_ex3, 1);
//	normalizeImage(real_ex3,2*N,2*M);
//	normalizeImage(imag_ex3,2*N,2*M);
	for(int i = 0; i < 2*N; i++){
		for(int j = 0; j < 2*M; j++){
			filtered.setPixelVal(i,j, real_ex3[i][j]);
		}
	}
	//magn(real_ex3, imag_ex3, filtered, 2*N,2*M);
	centerIt(filtered);
	writeImage("motion_blurred.pgm", filtered);
}
