from PIL import Image
import numpy as np
import math
import matplotlib.pyplot as plt
from FFT import *


def Experiment2():
  # generateImg(512, 512)
  DFT2D()


def applyFiltering(array1, array2, Height, Width):
    newarray = array1
    for rows in range(Height):
        for cols in range(Width):
            newarray[rows][cols] = (array1[rows][cols]) * (array2[rows][cols])
            # if(cols % 2 and rows % 2):
            #   newarray[rows][cols] = 0.0
    return newarray

# Does the FFT of a 2D image
def DFT2D():

    image = Image.open("newLenna.png")
    pixels = list(image.getdata())
    width, height = image.size
    newImage = Image.new("L", (width, height))
    pixels = [pixels[i * width:(i + 1) * width] for i in range(height)]

    image1 = Image.open("SobelPadded512.png")
    pixels1 = list(image1.getdata())
    width1, height1 = image1.size
    newImage1 = Image.new("L", (width1, height1))
    pixels1 = [pixels1[i * width1:(i + 1) * width1] for i in range(height1)]

    # # N = height x width
    N = (newImage.size[0] + newImage.size[1]) // 4

    print(N, width, height)
    isign = -1

    pixels1[1][1] = 1.0
    pixels1[1][2] = 0.0
    pixels1[1][3] = -1.0

    pixels1[2][1] = 2.0
    pixels1[2][2] = 0.0
    pixels1[2][3] = -2.0

    pixels1[3][1] = 1.0
    pixels1[3][2] = 0.0
    pixels1[3][3] = -1.0

    pixels = normalizeImage(pixels, height, width, N)
    pixels1 = normalizeImage(pixels1, height, width, N)

    # Iterate over all the rows and store it into test2D
    pixels = ApplyFFTRow(pixels, width, N, isign)
    pixels1 = ApplyFFTRow(pixels1, width, N, isign)

    # Iterate over all the columns and store it into pixels
    pixels = ApplyFFTCol(pixels, height, N, isign)
    pixels1 = ApplyFFTCol(pixels1, height, N, isign)

    pixels = applyFiltering(pixels, pixels1, height, width)

    minVal = 1000000000000000000000000000.0
    maxVal = -101010100000000000000000000.0

    # newImage = LinearScaleValues(newImage, maxVal, minVal, pixels, height, width)
    # newImage1 = LinearScaleValues(newImage1, maxVal, minVal, pixels1, height, width)

    # newImage = LogScaleValues(newImage, maxVal, minVal, pixels, height, width)
    # newImage1 = LogScaleValues(newImage1, maxVal, minVal, pixels1, height, width)

    ################################## Do the reverse ##################################
    # Iterate over all the columns and store it into pixels
    pixels = ApplyFFTCol(pixels, height, N, 1)

    # Iterate over all the rows and store it into test2D
    pixels = ApplyFFTRow(pixels, width, N, 1)


    # # pixels = changeAmplitudeSpatial(pixels, height, width, N, 1)
    newImage = LinearScaleValues(newImage, maxVal, minVal, pixels, height, width)


    # # Only if we don't linearScale to see what it looks like.
    # for rows in range(height):
    #     for cols in range(width):
    #         val = int(pixels[rows][cols])
    #         val1 = int(pixels1[rows][cols])
    #         val2 = int(newPixels[rows][cols])

    #         newImage.putpixel((cols,rows), val)
    #         newImage1.putpixel((cols,rows), val1)
    #         newImage2.putpixel((cols,rows), val2)

    newImage.save("lenna_Sobel.png")
    # newImage1.save("./data_output/SobelFFT.png")



Experiment2()
