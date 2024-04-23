#include <iostream>
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

int* parallelLowPassFilter(int* input, int width, int height, int kernelSize) {
    // Pad the input image
    int* paddedImage = padImage(input, width, height, kernelSize);

    // Allocate memory for the filtered image (ARGB for each pixel)
    int* filteredImage = new int[width * height];

    // Apply the low-pass filter
    for (int i = kernelSize / 2; i < height + kernelSize / 2; ++i) {
        for (int j = kernelSize / 2; j < width + kernelSize / 2; ++j) {
            int sumA = 0, sumR = 0, sumG = 0, sumB = 0;
            for (int m = -kernelSize / 2; m <= kernelSize / 2; ++m) {
                for (int n = -kernelSize / 2; n <= kernelSize / 2; ++n) {
                    // Get the ARGB value from the padded image
                    int argb = paddedImage[(i + m) * (width + kernelSize - 1) + (j + n)];
                    // Extract individual components
                    int a = (argb >> 24) & 0xFF;
                    int r = (argb >> 16) & 0xFF;
                    int g = (argb >> 8) & 0xFF;
                    int b = argb & 0xFF;
                    // Accumulate sum of each component
                    sumA += a;
                    sumR += r;
                    sumG += g;
                    sumB += b;
                }
            }
            // Compute average values for each component
            int avgA = sumA / (kernelSize * kernelSize);
            int avgR = sumR / (kernelSize * kernelSize);
            int avgG = sumG / (kernelSize * kernelSize);
            int avgB = sumB / (kernelSize * kernelSize);
            // Combine average values into ARGB format
            int avgARGB = (avgA << 24) | (avgR << 16) | (avgG << 8) | avgB;
            // Store the average ARGB value in the filtered image
            filteredImage[(i - kernelSize / 2) * width + (j - kernelSize / 2)] = avgARGB;
        }
    }

    // Free memory allocated for padded image
    delete[] paddedImage;

    return filteredImage;
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

int main()
{
    int ImageWidth = 4, ImageHeight = 4;

    int start_s, stop_s, TotalTime = 0;

    System::String^ imagePath;
    std::string img;
    img = "./input_image.jpg";

    imagePath = marshal_as<System::String^>(img);
    int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);

    start_s = clock();

    // Apply the parallel low-pass filter
    int kernelSize = 3;
    int* filteredImageData = parallelLowPassFilter(imageData, ImageWidth, ImageHeight, kernelSize);

    stop_s = clock();
    TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

    // Create the filtered image
    createImage(filteredImageData, ImageWidth, ImageHeight, 1);
    cout << "time: " << TotalTime << " milliseconds" << endl;

    // Free memory
    delete[] imageData;
    delete[] filteredImageData;

    return 0;
}
