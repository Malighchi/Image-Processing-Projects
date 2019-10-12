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
void avgfilter(ImageType, ImageType&, int, int, int, float[]);
void gausfilter(ImageType, ImageType&, int, int, int, float[]);

int main(int argc, char *argv[]){
    int N, M, Q;
    bool type;

    stringstream strfiltersize(argv[2]);
    int filtersize=0;
    strfiltersize>>filtersize;
    
    float avgmask7 [49];
    float avgmask15 [225];

    for(int i = 0; i<49; i++){
        avgmask7[i] = ((float)1/(float)49);
    }
    for(int i = 0; i<225; i++){
        avgmask15[i] = ((float)1/(float)225);
    }

    float gausmask7 [49] = {1, 1, 2, 2, 2, 1, 1,
                          1, 2, 2, 4, 2, 2, 1,
                          2, 2, 4, 8, 4, 2, 2,
                          2, 4, 8, 16, 8, 4, 2,
                          2, 2, 4, 8, 4, 2, 2,
                          1, 2, 2, 4, 2, 2, 1,
                          1, 1, 2, 2, 2, 1, 1};
                          
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

    int sum7 = 0;
    for(int i = 0; i<49; i++){
        sum7+=gausmask7[i];
    }
    for(int i = 0; i<49; i++){
        gausmask7[i] = gausmask7[i]/sum7;
    }

    int sum15 = 0;
    for(int i = 0; i<225; i++){
        sum15+=gausmask15[i];
    }
    for(int i = 0; i<225; i++){
        gausmask15[i]=gausmask15[i]/sum15;
    }

    // read image header
    readImageHeader(argv[1], N, M, Q, type);

    // allocate memory for the image array
    ImageType image(N, M, Q);
    
    // read image
    readImage(argv[1], image);

    ImageType filtimgavg(image);
    ImageType filtimggaus(image);

    // write image
    char *lenna = (char*)"lenna.pgm";
    if(strcmp(argv[1], lenna) == 0){
        if(filtersize == 7){
            avgfilter(image, filtimgavg, filtersize, N, M, avgmask7);
            gausfilter(image, filtimggaus, filtersize, N, M, gausmask7);
            char *avgname = (char*)"lennaavgfiltered7.pgm";
            char *gausname = (char*)"lennagausfiltered7.pgm";
            writeImage(avgname, filtimgavg);
            writeImage(gausname, filtimggaus);
        }
        else{
            avgfilter(image, filtimgavg, filtersize, N, M, avgmask15);
            gausfilter(image, filtimggaus, filtersize, N, M, gausmask15);
            char *avgname = (char*)"lennaavgfiltered15.pgm";
            char *gausname = (char*)"lennagausfiltered15.pgm";
            writeImage(avgname, filtimgavg);
            writeImage(gausname, filtimggaus);
        }

    }
    else{
        if(filtersize == 7){
            avgfilter(image, filtimgavg, filtersize, N, M, avgmask7);
            gausfilter(image, filtimggaus, filtersize, N, M, gausmask7);
            char *avgname = (char*)"sfavgfiltered7.pgm";
            char *gausname = (char*)"sfgausfiltered7.pgm";
            writeImage(avgname, filtimgavg);
            writeImage(gausname, filtimggaus);
        }
        else{
            avgfilter(image, filtimgavg, filtersize, N, M, avgmask15);
            gausfilter(image, filtimggaus, filtersize, N, M, gausmask15);
            char *avgname = (char*)"sfavgfiltered15.pgm";
            char *gausname = (char*)"sfgausfiltered15.pgm";
            writeImage(avgname, filtimgavg);
            writeImage(gausname, filtimggaus);
        }
    }
    

    return (1);

}

void avgfilter(ImageType img, ImageType& filterimg, int filtersize, int N, int M, float mask[]){
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