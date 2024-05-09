#include <iostream>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>
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

//example 
   //double kernel[3][3] = {
     // {1.0 / 16, 2.0 / 16, 1.0 / 16},
      // {2.0 / 16, 4.0 / 16, 2.0 / 16},
      // {1.0 / 16, 2.0 / 16, 1.0 / 16}
   //};
double** generate2DGaussianKernel(int size, double sigma) { //sigma : standard deviation
    double** kernel = new double* [size];
    double sum = 0.0;
    int halfSize = size / 2;
   
    // Allocate memory for rows
    for (int i = 0; i < size; ++i) {
        kernel[i] = new double[size];
    }
    for (int i = -halfSize; i <= halfSize; ++i) {
        for (int j = -halfSize; j <= halfSize; ++j) {
        // Calculate the Gaussian value for each position in the kernel,represent the relative positions from the center of the kernel
            kernel[i + halfSize][j + halfSize] = exp(-(i * i + j * j) / (2 * sigma * sigma));
            sum += kernel[i + halfSize][j + halfSize];
        }
    }
    // Normalize the kernel
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            kernel[i][j] /= sum;
        }
    }

    return kernel;
}

/*It extracts the width and height of the image and stores them in the variables pointed to by w and h.
Memory is dynamically allocated for the input image array based on the image dimensions.
Pixel values are read from the image using the GetPixel method, and
their grayscale values are calculated by averaging the RGB components and stored in the allocated input array.
The function returns a pointer to the array containing grayscale pixel values.
*/

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
/*
This function creates a new bitmap image using the provided pixel values, width, and height.
It initializes a System::Drawing::Bitmap object with the specified width and height.
Iterating over each pixel in the image, it retrieves the grayscale value from the input pixel array.
Using the grayscale value, it creates a System::Drawing::Color object and sets the corresponding pixel in the new image.
Finally, the new image is saved to our output directory .

*/
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
    MyNewImage.Save("..//Data//Output//outputRes" + ".png");
    cout << "result Image Saved " << index << endl;
}


/*The function fills the padded image array with zeros.
It copies the original image into the center of the padded image array, leaving the borders zero-padded.
The ARGB value of the pixel from the original image is stored in the padded image at an adjusted position.
Since we're padding the image, we need to shift the pixel positions accordingly.
(i + kernelSize / 2) and (j + kernelSize / 2) shift the pixel position by kernelSize / 2 in both the row and column directions, respectively.
This ensures that the original image is placed in the center of the padded image.
*/
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


int* applyBlurFilterSequential(int* input, int width, int height, int kernelSize, int sigma) {
    // Generate Gaussian kernel
    double** kernel = generate2DGaussianKernel(kernelSize, sigma);

    // Allocate memory for filtered image
    int* filteredImage = new int[width * height];

    // Pad the input image
    int paddedWidth = width + 2 * (kernelSize / 2);
    int paddedHeight = height + 2 * (kernelSize / 2);
    int* paddedImage = padImage(input, width, height, kernelSize);

    // Apply Gaussian blur using padded image
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double sumR = 0, sumG = 0, sumB = 0;
            double weightSum = 0;

            for (int m = -kernelSize / 2; m <= kernelSize / 2; ++m) {
                for (int n = -kernelSize / 2; n <= kernelSize / 2; ++n) {
                    int indexX = j + n + (kernelSize / 2);
                    int indexY = i + m + (kernelSize / 2);
                    int argb = paddedImage[indexY * paddedWidth + indexX];
                    double weight = kernel[m + kernelSize / 2][n + kernelSize / 2];
                    sumR += weight * ((argb >> 16) & 0xFF);
                    sumG += weight * ((argb >> 8) & 0xFF);
                    sumB += weight * (argb & 0xFF);
                    weightSum += weight;
                }
            }

            if (weightSum > 0) {
                sumR /= weightSum;
                sumG /= weightSum;
                sumB /= weightSum;
            }

            int filteredPixel = (((int)sumR) << 16) | (((int)sumG) << 8) | ((int)sumB);
            filteredImage[i * width + j] = filteredPixel;
        }
    }

    // Cleanup memory
    delete[] paddedImage;
    for (int i = 0; i < kernelSize; ++i) {
        delete[] kernel[i];
    }
    delete[] kernel;

    // Return filtered image
    return filteredImage;
}

/*
OpenMP parallelization is utilized to iterate over each pixel in the input image.
Within each pixel iteration, the function performs convolution with the Gaussian kernel to calculate the filtered pixel value.
The filtered pixel values are stored in a new array allocated for the filtered image.
After processing all pixels, memory allocated for the padded image and Gaussian kernel is deallocated.
The function returns a pointer to the array containing filtered pixel values.
*/
int* applyBlurFilterOpenMP(int* input, int width, int height, int kernelSize, int sigma) {
    // Generate Gaussian kernel
    double** kernel = generate2DGaussianKernel(kernelSize, sigma);

    // Allocate memory for filtered image
    int* filteredImage = new int[width * height];

    // Pad the input image
    int paddedWidth = width + 2 * (kernelSize / 2);
    int paddedHeight = height + 2 * (kernelSize / 2);
    int* paddedImage = padImage(input, width, height, kernelSize);

    omp_set_num_threads(8); // set to 8 threads

    #pragma omp parallel for collapse(2) schedule(static) shared(input,width,height,kernelSize,sigma)
    // Apply Gaussian blur using padded image
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double sumR = 0, sumG = 0, sumB = 0;
            double weightSum = 0;

            for (int m = -kernelSize / 2; m <= kernelSize / 2; ++m) {
                for (int n = -kernelSize / 2; n <= kernelSize / 2; ++n) {
                    int indexX = j + n + (kernelSize / 2);
                    int indexY = i + m + (kernelSize / 2);
                    int argb = paddedImage[indexY * paddedWidth + indexX];
                    double weight = kernel[m + kernelSize / 2][n + kernelSize / 2];
                    sumR += weight * ((argb >> 16) & 0xFF);
                    sumG += weight * ((argb >> 8) & 0xFF);
                    sumB += weight * (argb & 0xFF);
                    weightSum += weight;
                }
            }

            if (weightSum > 0) {
                sumR /= weightSum;
                sumG /= weightSum;
                sumB /= weightSum;
            }

            int filteredPixel = (((int)sumR) << 16) | (((int)sumG) << 8) | ((int)sumB);
            filteredImage[i * width + j] = filteredPixel;
        }
    }

    // Cleanup memory
    delete[] paddedImage;
    for (int i = 0; i < kernelSize; ++i) {
        delete[] kernel[i];
    }
    delete[] kernel;

    // Return filtered image
    return filteredImage;
}
int main()
{
    int ImageWidth = 4, ImageHeight = 4;

    int start_s, stop_s, TotalTime = 0;

    System::String^ imagePath;
    std::string img;
    img = "..//Data//Input//lena.png";

    imagePath = marshal_as<System::String^>(img);
    int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
    start_s = clock();
    int kernelSize = 5; // Adjust kernel size 
    double sigma =2; // Adjust sigma (standard deviation) as needed
    //int* blurredImageData1 = applyBlurFilterSequential(imageData, ImageWidth, ImageHeight, kernelSize, sigma);
    int* blurredImageData = applyBlurFilterOpenMP(imageData, ImageWidth, ImageHeight, kernelSize, sigma);
    stop_s = clock();    

    TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

    // Create the filtered (blurred) image
    createImage(blurredImageData, ImageWidth, ImageHeight, 0);
     
    cout << "time: " << TotalTime << " milliseconds" << endl;

    // Free memory
    delete[] imageData;
    delete[] blurredImageData;
    system("pause");
    return 0;
}
