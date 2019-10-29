#include "fft.cpp"
#include <stdio.h>

int main (){
  float data[4] = {2,3,4,4};
  fft(data,4,-1);
  for(int i = 0; i < 4; i++){
     printf("%f\n", data[i]);
  }
}
