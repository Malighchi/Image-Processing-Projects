#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <ctime>
#include "image.h"
using namespace std;

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);
void corruptimg(ImageType&, float, int, int);
void medianfilter(ImageType, ImageType&, int, int, int);

int main(int argc, char *argv[]){
    int N, M, Q;
    bool type;

    stringstream strfiltersize(argv[2]);
    int filtersize=0;
    strfiltersize>>filtersize;

    // read image header
    readImageHeader(argv[1], N, M, Q, type);

    // allocate memory for the image array
    ImageType image(N, M, Q);
    
    // read image
    readImage(argv[1], image);

    ImageType corrimage30(image);
    ImageType corrimage50(image);

    corruptimg(corrimage30, 0.30, N, M);
    corruptimg(corrimage50, 0.50, N, M);

    ImageType filtimg30(corrimage30);
    ImageType filtimg50(corrimage50);

    // write image
    char *lenna = (char*)"lenna.pgm";
    if(strcmp(argv[1], lenna) == 0){
        char *corr30name = (char*)"lenna30sp.pgm";
        char *corr50name = (char*)"lenna50sp.pgm";
        writeImage(corr30name, corrimage30);
        writeImage(corr50name, corrimage50);

        medianfilter(corrimage30, filtimg30, filtersize, N, M);
        medianfilter(corrimage50, filtimg50, filtersize, N, M);
        if(filtersize == 7){
            char *median30name = (char*)"lenna30filtered7.pgm";
            char *median50name = (char*)"lenna50filtered7.pgm";
            writeImage(median30name, filtimg30);
            writeImage(median50name, filtimg50);
        }
        else{
            char *median30name = (char*)"lenna30filtered15.pgm";
            char *median50name = (char*)"lenna50filtered15.pgm";
            writeImage(median30name, filtimg30);
            writeImage(median50name, filtimg50);
        }

    }
    else{
        char *corr30name = (char*)"boat30sp.pgm";
        char *corr50name = (char*)"boat50sp.pgm";
        writeImage(corr30name, corrimage30);
        writeImage(corr50name, corrimage50);

        medianfilter(corrimage30, filtimg30, filtersize, N, M);
        medianfilter(corrimage50, filtimg50, filtersize, N, M);
        if(filtersize == 7){
            char *median30name = (char*)"boat30filtered7.pgm";
            char *median50name = (char*)"boat50filtered7.pgm";
            writeImage(median30name, filtimg30);
            writeImage(median50name, filtimg50);
        }
        else{
            char *median30name = (char*)"boat30filtered15.pgm";
            char *median50name = (char*)"boat50filtered15.pgm";
            writeImage(median30name, filtimg30);
            writeImage(median50name, filtimg50);
        }
    }
    

    return (1);

}

void corruptimg(ImageType& corrimg, float c_factor, int N, int M){
    srand(time(0));
    int corr = (M*N)*c_factor;

    for(int i = 0; i < corr; i++){
        int white_black = (rand() % 2);
        int randx = (rand() % N);
        int randy = (rand() % M);
        if(white_black == 1){
            corrimg.setPixelVal(randx, randy, 255);
        } 
        else{
            corrimg.setPixelVal(randx, randy, 0);
        }
    }
}

void medianfilter(ImageType corrimg, ImageType& filtimg, int masksize, int N, int M){
    int mask_arr[masksize*masksize];
    int curr_mask = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            for(int s = (-masksize/2); s <= (masksize/2); s++){
                for(int t = (-masksize/2); t <= (masksize/2); t++){
                    int value = 0;
                    if(((i+s)<0 || (i+s)>=N) || ((j+t)<0 || (j+t)>=M)){
                        mask_arr[curr_mask] = 0;
                    }
                    else{
                        corrimg.getPixelVal(i+s, j+t, value);
                        mask_arr[curr_mask] = value;
                    }
                    curr_mask++;
                }
            }
            sort(&mask_arr[0], &mask_arr[masksize*masksize]);
            
            if((masksize*masksize) % 2 == 0){
                int medianvalue = (mask_arr[(masksize*masksize)/2] + mask_arr[((masksize*masksize)/2)-1])/2;
                filtimg.setPixelVal(i, j, medianvalue);
            }
            else{
                filtimg.setPixelVal(i, j, mask_arr[((masksize*masksize)/2)]);
            }

            
            curr_mask = 0;
        }

    }
}