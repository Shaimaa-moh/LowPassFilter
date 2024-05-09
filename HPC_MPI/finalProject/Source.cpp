#include <iostream>
#include <cmath>
#include <cstring>
#include <msclr\marshal_cppstd.h>
#include <omp.h>
#include <mpi.h>
#pragma once
#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>

using namespace std;
using namespace msclr::interop;

double** generate2DGaussianKernel(int size, double sigma) {
    double** kernel = new double* [size];
    double sum = 0.0;
    int halfSize = size / 2;

    for (int i = 0; i < size; ++i) {
        kernel[i] = new double[size];
    }
    for (int i = -halfSize; i <= halfSize; ++i) {
        for (int j = -halfSize; j <= halfSize; ++j) {
            kernel[i + halfSize][j + halfSize] = exp(-(i * i + j * j) / (2 * sigma * sigma));
            sum += kernel[i + halfSize][j + halfSize];
        }
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            kernel[i][j] /= sum;
        }
    }
    return kernel;
}

int* inputImage(int* w, int* h, System::String^ imagePath) {
    int* input;
    int OriginalImageWidth, OriginalImageHeight;

    System::Drawing::Bitmap BM(imagePath);

    OriginalImageWidth = BM.Width;
    OriginalImageHeight = BM.Height;
    *w = BM.Width;
    *h = BM.Height;
    input = new int[BM.Height * BM.Width];
    for (int i = 0; i < BM.Height; i++) {
        for (int j = 0; j < BM.Width; j++) {
            System::Drawing::Color c = BM.GetPixel(j, i);
            input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3);
        }
    }
    return input;
}

void createImage(int* image, int width, int height, int index) {
    System::Drawing::Bitmap MyNewImage(width, height);
    for (int i = 0; i < MyNewImage.Height; i++) {
        for (int j = 0; j < MyNewImage.Width; j++) {
            int grayscale = image[i * width + j];
            grayscale = max(0, min(255, grayscale));
            System::Drawing::Color c = System::Drawing::Color::FromArgb(grayscale, grayscale, grayscale);
            MyNewImage.SetPixel(j, i, c);
        }
    }
    MyNewImage.Save("..//Data//Output//outputRes" + ".jpg");
    cout << "result Image Saved " << index << endl;
}

int* padImage(int* input, int width, int height, int kernelSize) {
    int paddedWidth = width + 2 * (kernelSize / 2);
    int paddedHeight = height + 2 * (kernelSize / 2);
    int* paddedImage = new int[paddedWidth * paddedHeight];
    memset(paddedImage, 0, sizeof(int) * paddedWidth * paddedHeight);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int argb = input[i * width + j];
            paddedImage[(i + kernelSize / 2) * paddedWidth + (j + kernelSize / 2)] = argb;
        }
    }
    return paddedImage;
}

int* applyBlurFilter(int* input, int width, int height, int kernelSize, int sigma) {
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




void parallelLowPassFilter(int* imageData, int ImageWidth, int ImageHeight, int kernelSize, int index, int rank, int size) {
   // MPI_Bcast(&ImageWidth, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&ImageHeight, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double start_s, stop_s, TotalTime = 0;
    int rowsPerProcess = ImageHeight / size;
    int startRow = rank * rowsPerProcess;
    int endRow = startRow + rowsPerProcess;

    int* localImageData = new int[ImageWidth * rowsPerProcess];
    start_s = clock();
    MPI_Scatter(imageData + startRow * ImageWidth, ImageWidth * rowsPerProcess, MPI_INT,
        localImageData, ImageWidth * rowsPerProcess, MPI_INT, 0, MPI_COMM_WORLD);

    int sigma = 2;

    int* filteredData = applyBlurFilter(localImageData, ImageWidth, rowsPerProcess, kernelSize, sigma);

  
    
    //if (rank == 0) {
        int* gatheredData = new int[ImageWidth * ImageHeight];
    //}
    MPI_Gather(filteredData, ImageWidth * rowsPerProcess, MPI_INT,
        gatheredData + startRow * ImageWidth, ImageWidth * rowsPerProcess, MPI_INT, 0, MPI_COMM_WORLD);
    stop_s = clock();
    
    if (rank == 0) {
        createImage(gatheredData, ImageWidth, ImageHeight, 1);
        TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
        cout << "time: " << TotalTime << endl;
        delete[] gatheredData;
    }

    delete[] localImageData;
    delete[] filteredData;
    if (rank == 0) {
        delete[] imageData;
    }
}

int main(int argc, char* argv[]) {
    int ImageWidth, ImageHeight;
    System::String^ imagePath;
    std::string img = "..//Data//Input//lena.png";
    imagePath = marshal_as<System::String^>(img);
    int *imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
   

    

    //------------------------------------- Low Pass Filter -------------------------------------

    parallelLowPassFilter(imageData, ImageWidth, ImageHeight, 5, 0, rank, size);

    //-------------------------------------------------------------------------------------------

    
   
  
    MPI_Finalize();
}