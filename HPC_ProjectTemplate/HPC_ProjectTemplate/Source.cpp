#include <iostream>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <msclr\marshal_cppstd.h>
#include <ctime>
#pragma once 
#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int* padImage(int* input, int width, int height, int kernelSize) {
    int paddedWidth = width + 2 * (kernelSize / 2);
    int paddedHeight = height + 2 * (kernelSize / 2);

    // Allocate memory for padded image (ARGB for each pixel)
    int* paddedImage = new int[paddedWidth * paddedHeight];

    // Fill the padded image with zeros
    memset(paddedImage, 0, sizeof(int) * paddedWidth * paddedHeight);

    // Copy the original image into the padded image
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            // Get the ARGB value from the input array
            int argb = input[i * width + j];
            // Store the ARGB value in the padded image
            paddedImage[(i + kernelSize / 2) * paddedWidth + (j + kernelSize / 2)] = argb;
        }
    }

    return paddedImage;
}

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
    int OriginalImageWidth, OriginalImageHeight;

    //**Read Image and save it to local arrays**	
    //Read Image and save it to local arrays
    System::Drawing::Bitmap BM(imagePath);

    OriginalImageWidth = BM.Width;
    OriginalImageHeight = BM.Height;
    *w = OriginalImageWidth;
    *h = OriginalImageHeight;

    // Allocate memory for input image (grayscale)
    int* input = new int[OriginalImageWidth * OriginalImageHeight];

    // Read the pixel values and convert to grayscale
    for (int i = 0; i < OriginalImageHeight; i++) {
        for (int j = 0; j < OriginalImageWidth; j++) {
            System::Drawing::Color c = BM.GetPixel(j, i);
            // Calculate grayscale value (average of RGB components)
            int grayscale = (c.R + c.G + c.B) / 3;
            // Store the grayscale value in the input array
            input[i * OriginalImageWidth + j] = grayscale;
        }
    }
    return input;
}

void createImage(int* image, int width, int height, int index)
{
    System::Drawing::Bitmap MyNewImage(width, height);

    // Iterate over each pixel in the image
    for (int i = 0; i < MyNewImage.Height; i++) {
        for (int j = 0; j < MyNewImage.Width; j++) {
            // Set the pixel color in the new image
            int grayscale = image[i * width + j];
            // Create a Color object with grayscale value for RGB components
            System::Drawing::Color c = System::Drawing::Color::FromArgb(grayscale, grayscale, grayscale);
            MyNewImage.SetPixel(j, i, c);
        }
    }

    // Save the image
    MyNewImage.Save("..//Data//Output//outputRes" + ".jpg");
    cout << "result Image Saved " << index << endl;
}
int* applyBlurFilter(int* input, int width, int height,int kernelSize, int sigma) {
    int* paddedImage = padImage(input, width, height, kernelSize);

    // Allocate memory for the filtered image
    int* filteredImage = new int[width * height];
    // Calculate Gaussian kernel
    double* kernel = new double[kernelSize];
    double sum = 0.0;
    int halfSize = kernelSize / 2;
    for (int i = -halfSize; i <= halfSize; ++i) {
        kernel[i + halfSize] = exp(-(i * i) / (2 * sigma * sigma));
        sum += kernel[i + halfSize];
    }
    for (int i = 0; i < kernelSize; ++i) {
        kernel[i] /= sum;
    }
    //double kernel[3][3] = {
      //  {1.0 / 16, 2.0 / 16, 1.0 / 16},
       // {2.0 / 16, 4.0 / 16, 2.0 / 16},
       // {1.0 / 16, 2.0 / 16, 1.0 / 16}
    //};

#pragma omp parallel for collapse(2) schedule(static) shared(input, width, height,kernelSize,sigma )

   // Apply the Gaussian-like blur filter
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double sumR = 0, sumG = 0, sumB = 0;

            // Convolve with the predefined kernel
            for (int m = -1; m <= 1; ++m) {
                for (int n = -1; n <= 1; ++n) {
                    int indexX = j + n;
                    int indexY = i + m;
                    if (indexX >= 0 && indexX < width && indexY >= 0 && indexY < height) {
                        int argb = paddedImage[(i + m) * (width + kernelSize - 1) + (j + n)];
                        double weight = kernel[m + halfSize] * kernel[n + halfSize];
                        sumR += weight * ((argb >> 16) & 0xFF);
                        sumG += weight * ((argb >> 8) & 0xFF);
                        sumB += weight * (argb & 0xFF);
                    }
                }
            }

            // Combine the weighted sums to get the filtered pixel value
            int filteredPixel = ((int)sumR << 16) | ((int)sumG << 8) | (int)sumB;
            filteredImage[i * width + j] = filteredPixel;
        }
    }

    return filteredImage;
}


int main()
{
    int ImageWidth = 4, ImageHeight = 4;

    int start_s, stop_s, TotalTime = 0;

    System::String^ imagePath;
    std::string img;
    img = "..//Data//Input//test1.jpg";

    imagePath = marshal_as<System::String^>(img);
    int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);

    start_s = clock();

    // Apply the blur filter
    int kernelSize = 3; // Adjust kernel size 
    double sigma = 1.0; // Adjust sigma (standard deviation) as needed

    int* blurredImageData = applyBlurFilter(imageData, ImageWidth, ImageHeight, kernelSize,sigma);

    stop_s = clock();
    TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

    // Create the filtered (blurred) image
    createImage(blurredImageData, ImageWidth, ImageHeight, 1);
    cout << "time: " << TotalTime << " milliseconds" << endl;

    // Free memory
    delete[] imageData;
    delete[] blurredImageData;
    system("pause");
    return 0;
}
