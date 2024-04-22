#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input; //point to the memory block allocated to store the pixel values of the image.


	int OriginalImageWidth, OriginalImageHeight;

	//********************Read Image and save it to local arrayss********	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;

	//integer arrays Red, Green, and Blue, each of size BM.Height * BM.Width, to store the individual color channels of the image.

	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width]; //Allocate memory for the input array to store the grayscale values of the image.
#pragma omp parallel for

	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			// calculate the grayscale value as the average of the RGB value
			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3);

		}

	}
	// points to the memory block containing the grayscale values of the image.

	return input;
}


void createImage(int* image, int width, int height, int index)
{
	// *image : A pointer to an array containing pixel values of the image.

	System::Drawing::Bitmap MyNewImage(width, height); // Creates a new Bitmap object named MyNewImage with the specified widthand height.This will be used to create the new image.

	//it checks each pixel's intensity (stored in the image array) to ensure it falls within the valid range of 0 to 255.
#pragma omp parallel for

	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//If a pixel's intensity is less than 0, it sets it to 0. 
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			//If it's greater than 255, it sets it to 255. 
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			//It then creates a new Color object (c) using the intensity value for all three RGB components. This effectively creates a grayscale image.

			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}
int* lowPassFilter(int* input, int height, int width)
{
	int* output = new int[height * width];

	// Gaussian kernel weights
	double kernel[3][3] = {
		{1.0 / 16, 2.0 / 16, 1.0 / 16},
		{2.0 / 16, 4.0 / 16, 2.0 / 16},
		{1.0 / 16, 2.0 / 16, 1.0 / 16}
	};

	int numThreads = 4; // Adjust as needed
	omp_set_num_threads(numThreads);

	// Apply the filter in parallel
      #pragma omp parallel for shared(input, output) num_threads(numThreads)

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			double total = 0;

			// Convolve with the kernel ,For each pixel loop iterates over the neighboring pixels in a 3x3 region (the kernel).

			//loop over kernel rows 
			for (int m = -1; m <= 1; m++)
			{
				// loop over kernel cols
				for (int n = -1; n <= 1; n++)
				{

					int row = i + m; // adds the relative row coordinate m to the current row coordinate i.
					//It calculates the row index of the pixel in the input image that corresponds to the current position in the kernel.

					int col = j + n;

					//Condition to check  that only valid pixel positions within the image bounds are accessed.

					if (row >= 0 && row < height && col >= 0 && col < width)
					{
						//The input pixel values are multiplied by the corresponding kernel elements.

						//Adding 1 to m and n  map correctly to the indices 0, 1, and 2 in the kernel array.
						// input[ row*width + col] The resulting value represents the index in the input image array where the pixel value corresponding to the current position in the kernel is located.

						total += input[row * width + col] * kernel[m + 1][n + 1];
					}
				}
			}
			// This line assigns the integer value of total to the corresponding pixel in the output image after the convolution is complete.
			output[i * width + j] = (int)total;
		}
	}

	return output;
}

int main()
{
	int ImageWidth = 4, ImageHeight = 4;

	int start_s, stop_s, TotalTime = 0;

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";
	start_s = clock();

	imagePath = marshal_as<System::String^>(img);
	int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
	int* FilterOutput = lowPassFilter(imageData, ImageWidth, ImageHeight);


	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	createImage(FilterOutput, ImageWidth, ImageHeight, 2);
	cout << "time: " << TotalTime << endl;

	free(FilterOutput);
	return 0;

}