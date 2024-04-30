#include <iostream>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <msclr\marshal_cppstd.h>
#include <ctime>
#include <mpi.h>
#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>

using namespace std;
using namespace msclr::interop;

// Function to generate a 1D Gaussian kernel
double* generateGaussianKernel(int size, double sigma) {
    double* kernel = new double[size];
    double sum = 0.0;
    int halfSize = size / 2;

    for (int i = -halfSize; i <= halfSize; ++i) {
        kernel[i + halfSize] = exp(-(i * i) / (2 * sigma * sigma));
        sum += kernel[i + halfSize];
    }

    // Normalize the kernel
    for (int i = 0; i < size; ++i) {
        kernel[i] /= sum;
    }

    return kernel;
}

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
    int OriginalImageWidth, OriginalImageHeight;

    //*********************************************************Read Image and save it to local arrayss*************************	
    //Read Image and save it to local arrayss
    System::Drawing::Bitmap BM(imagePath);

    OriginalImageWidth = BM.Width;
    OriginalImageHeight = BM.Height;
    *w = OriginalImageWidth;
    *h = OriginalImageHeight;

    // Allocate memory for input image
    int* input = new int[OriginalImageWidth * OriginalImageHeight];

    // Read the pixel values and store them in the input array
    for (int i = 0; i < OriginalImageHeight; i++) {
        for (int j = 0; j < OriginalImageWidth; j++) {
            System::Drawing::Color c = BM.GetPixel(j, i);
            // Store the green component of the pixel
            input[(i * OriginalImageWidth + j)] = c.ToArgb();
        }
    }
    return input;
}

void createImage(int* image, int width, int height, int index)
{
    System::Drawing::Bitmap MyNewImage(width, height);

    for (int i = 0; i < MyNewImage.Height; i++) {
        for (int j = 0; j < MyNewImage.Width; j++) {
            // Set the pixel color in the new image
            int argb = image[i * width + j];
            // Create a Color object from the ARGB value
            System::Drawing::Color c = System::Drawing::Color::FromArgb(argb);

            MyNewImage.SetPixel(j, i, c);
        }
    }

    // Save the image
    MyNewImage.Save("./output_image" + ".jpg");
    cout << "result Image Saved " << index << endl;
}

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

int* LowPassOmp(int* input, int width, int height, int kernelSize, double sigma) {
    // Allocate memory for the filtered image
    int* paddedImage = padImage(input, width, height, kernelSize);

    double* kernel = generateGaussianKernel(kernelSize, sigma);
    int* filteredImage = new int[width * height];
    // Apply the blur filter
    // Iterate over each pixel in the image
#pragma omp parallel for collapse(2) schedule(static) shared(input, width, height, kernelSize,sigma)

    for (int i = kernelSize / 2; i < height + kernelSize / 2; ++i) {
        for (int j = kernelSize / 2; j < width + kernelSize / 2; ++j) {
            double sumR = 0, sumG = 0, sumB = 0;

            // Convolve with the Gaussian kernel
            for (int m = -kernelSize / 2; m <= kernelSize / 2; ++m) {
                for (int n = -kernelSize / 2; n <= kernelSize / 2; ++n) {
                    int indexX = j + n;
                    int indexY = i + m;
                    if (indexX >= 0 && indexX < width + kernelSize && indexY >= 0 && indexY < height + kernelSize) {
                        int pixel = paddedImage[indexY * (width + kernelSize) + indexX];
                        double weight = kernel[m + kernelSize / 2] * kernel[n + kernelSize / 2];
                        sumR += weight * ((pixel >> 16) & 0xFF);
                        sumG += weight * ((pixel >> 8) & 0xFF);
                        sumB += weight * (pixel & 0xFF);
                    }
                }
            }

            // Combine the weighted sums to get the filtered pixel value
            int filteredPixel = ((int)sumR << 16) | ((int)sumG << 8) | (int)sumB;
            filteredImage[(i - kernelSize / 2) * width + (j - kernelSize / 2)] = filteredPixel;
        }
    }

    // Free memory allocated for padded image
    delete[] paddedImage;
    delete[] kernel;

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
    //sequentialLowPassFilter(imageData, ImageWidth, ImageHeight, 3, 0, rank, size);
    start_s = clock();
    int kernelsize = 3;
    double sigma = 1.0;
    //parallelLowPassFilter(imageData, ImageWidth, ImageHeight, 3, 0, rank, size);
    int* blurredImageData = LowPassOmp(imageData, ImageWidth, ImageHeight, kernelsize, sigma);
    stop_s = clock();
    TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

    // Create the filtered (blurred) image
    createImage(blurredImageData, ImageWidth, ImageHeight, 1);
    cout << "time: " << TotalTime << " milliseconds" << endl;

    // Free memory
    delete[] blurredImageData;
    system("pause");

    return 0;
}
