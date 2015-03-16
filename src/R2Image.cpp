// Source file for image class



// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"

//#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2,3.5);
	R2Point p2(2.1,2.2);
	R2Point p3(0.2,1.6);
	R2Point p4(0.0,0.5);
	R2Point p5(-0.2,4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1,5,1,6);

	linEquations[1][1] = p1[0]*p1[0];
	linEquations[1][2] = p1[0]*p1[1];
	linEquations[1][3] = p1[1]*p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0]*p2[0];
	linEquations[2][2] = p2[0]*p2[1];
	linEquations[2][3] = p2[1]*p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0]*p3[0];
	linEquations[3][2] = p3[0]*p3[1];
	linEquations[3][3] = p3[1]*p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;
	
	linEquations[4][1] = p4[0]*p4[0];
	linEquations[4][2] = p4[0]*p4[1];
	linEquations[4][3] = p4[1]*p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0]*p5[0];
	linEquations[5][2] = p5[0]*p5[1];
	linEquations[5][3] = p5[1]*p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n",p1[0],p1[1]);
	printf("Point #2: %f,%f\n",p2[0],p2[1]);
	printf("Point #3: %f,%f\n",p3[0],p3[1]);
	printf("Point #4: %f,%f\n",p4[0],p4[1]);
	printf("Point #5: %f,%f\n",p5[0],p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1,6,1,6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] + 
										p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] + 
										p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] + 
										p1[0]*nullspaceMatrix[4][smallestIndex] + 
										p1[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] + 
										p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] + 
										p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] + 
										p2[0]*nullspaceMatrix[4][smallestIndex] + 
										p2[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] + 
										p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] + 
										p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] + 
										p3[0]*nullspaceMatrix[4][smallestIndex] + 
										p3[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] + 
										p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] + 
										p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] + 
										p4[0]*nullspaceMatrix[4][smallestIndex] + 
										p4[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] + 
										p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] + 
										p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] + 
										p5[0]*nullspaceMatrix[4][smallestIndex] + 
										p5[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34,-2.8);

	printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] + 
											test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] + 
											test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] + 
											test_point[0]*nullspaceMatrix[4][smallestIndex] + 
											test_point[1]*nullspaceMatrix[5][smallestIndex] + 
											nullspaceMatrix[6][smallestIndex]);

	return;	
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{
	R2Image tempImgX(width,height);

	for(int i = 1; i < width-1;i++)
	{
		for(int j = 1; j < height-1;j++)
		{
			
			tempImgX.Pixel(i,j) = Pixel(i-1,j-1) + (Pixel(i-1,j)*2) + Pixel(i-1,j+1)- Pixel(i+1,j-1) - (Pixel(i+1,j)*2) - Pixel(i+1,j+1);
			//R2Pixel R2gray_pixel(0.5, 0.5, 0.5, 1.0);
			//tempImgX.Pixel(i,j) += R2gray_pixel;
		}
	}
	for(int i =1; i < width-1;i++)
	{
		for(int j =1; j < height-1;j++)
		{
			Pixel(i,j) = tempImgX.Pixel(i,j);
		//	Pixel(i,j).Clamp();
		}
	}
	// Apply the Sobel oprator to the image in X direction
  

}

void R2Image::
SobelY(void)
{

	// Apply the Sobel oprator to the image in Y direction
	
	R2Image tempImgY(width,height);
	for(int i = 1; i < width-1;i++)
	{
		for(int j = 1; j < height-1;j++)
		{
			tempImgY.Pixel(i,j) = Pixel(i-1,j-1)- Pixel(i-1,j+1) + (Pixel(i,j-1)*2)-(Pixel(i,j+1)*2) + (Pixel(i+1,j-1)) - (Pixel(i+1,j+1)) ;
			//R2Pixel R2gray_pixel(0.5, 0.5, 0.5, 1.0);
			//tempImgY.Pixel(i,j) += R2gray_pixel;
			
		}
	}
	for(int i = 1; i < width-1;i++)
	{
		for(int j = 1; j < height-1;j++)
		{
			Pixel(i,j) = tempImgY.Pixel(i,j);
			//Pixel(i,j).Clamp();
		}
	}
	

	

}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image
  
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred
    R2Image tempIMG(width, height);
    R2Image tempIMG2(width, height);
    int k = (int) (3 * sigma);
    double* sepK = new double[k];
    double weightSum = 0.0;
    for (int i = 0; i < k; i++){
        double distance = (i - k / 2)*(i - k / 2) + (k / 2)*(k / 2);
        double exponet = -(distance) / (2 * pow(sigma, 2));
        sepK[i] = pow(2.718, exponet);
        weightSum += sepK[i];
    }
 
    for (int y = k/2; y < height - k/2; y++) {
        for (int x = k/2; x < width - k/2; x++) {
                double valRed = 0;
                double valGreen = 0;
                double valBlue = 0;
                for (int j = 0; j < k; j++){
                    int index = y + (j - k / 2);
                    valRed = valRed + Pixel(x, index).Red() * sepK[j];
                    valGreen = valGreen + Pixel(x, index).Green() * sepK[j];
                    valBlue = valBlue + Pixel(x, index).Blue() * sepK[j];
                }
                valRed = valRed / weightSum;
                tempIMG.Pixel(x, y).SetRed(valRed);
                valGreen = valGreen / weightSum;
                tempIMG.Pixel(x, y).SetGreen(valGreen);
                valBlue = valBlue / weightSum;
                tempIMG.Pixel(x, y).SetBlue(valBlue);
        }
    }
               
    for (int y = k / 2; y < height - k / 2; y++) {
                    for (int x = k / 2; x < width - k / 2; x++) {
                                    double valRed = 0;
                                    double valGreen = 0;
                                    double valBlue = 0;
                                    for (int i = 0; i < k; i++){
                                        int index = x + (i - k / 2);
                                        valRed = valRed + tempIMG.Pixel(index, y).Red() * sepK[i];
                                        valGreen = valGreen + tempIMG.Pixel(index, y).Green() * sepK[i];
                                        valBlue = valBlue + tempIMG.Pixel(index, y).Blue() * sepK[i];
                                    }
                                    valRed = valRed / weightSum;
                                    tempIMG2.Pixel(x, y).SetRed(valRed);
                                    valGreen = valGreen / weightSum;
                                    tempIMG2.Pixel(x, y).SetGreen(valGreen);
                                    valBlue = valBlue / weightSum;
                                    tempIMG2.Pixel(x, y).SetBlue(valBlue);
                    }
    }
    for (int i = 1; i < width - 1; i++) {
                    for (int j = 1; j < height - 1; j++) {
                                    Pixel(i, j).SetRed(tempIMG2.Pixel(i, j).Red());
                                    Pixel(i, j).SetBlue(tempIMG2.Pixel(i, j).Blue());
                                    Pixel(i, j).SetGreen(tempIMG2.Pixel(i, j).Green());
                            //   Pixel(i, j).Clamp();
                    }
    }
}


void R2Image::
Harris(double sigma)
{
	
	R2Image tempImgX = *this;
	R2Image tempImgY = *this;

	tempImgX.SobelX();
	tempImgY.SobelY();

	R2Image tempImg1(width, height);
	R2Image tempImg2(width, height);
	R2Image tempImg3(width, height);
	
	for (int i = 1; i < width - 1; i++)
	{
	
		for (int j = 1; j < height - 1; j++)
		{
			//	tempImgX.SobelX();
			
			//	tempImgX.Pixel(i, j) *=tempImgX.Pixel(i, j);
			tempImg1.Pixel(i, j) = tempImgX.Pixel(i, j) *tempImgX.Pixel(i, j);
		//	tempImgY.Pixel(i, j) *= tempImgY.Pixel(i, j);
			
			//tempImg1.Pixel(i,j) += R2gray_pixel;
			tempImg2.Pixel(i, j) = tempImgY.Pixel(i, j)*tempImgY.Pixel(i, j);
		//	tempImg2.Pixel(i, j) += R2gray_pixel;
			
			tempImg3.Pixel(i, j) = tempImgX.Pixel(i, j)*tempImgY.Pixel(i, j);
			//tempImg3.Pixel(i, j) += R2gray_pixel;

		}
	}
	tempImg1.Blur(sigma);
	tempImg2.Blur(sigma);
	tempImg3.Blur(sigma);
	for (int i = 1; i < width - 1; i++)
	{
		for (int j = 1; j < height - 1; j++)
		{
		
			Pixel(i, j) = tempImg1.Pixel(i, j)*tempImg2.Pixel(i, j) - tempImg3.Pixel(i, j)*tempImg3.Pixel(i, j) -
			(0.04)*((tempImg1.Pixel(i,j)+tempImg2.Pixel(i,j))*(tempImg1.Pixel(i,j)+tempImg2.Pixel(i,j)));
			R2Pixel R2gray_pixel(0.5, 0.5, 0.5, 1.0);
			Pixel(i, j) += R2gray_pixel;
			
			Pixel(i, j).Clamp();
		}
	}
	
	
	int numFeat = 150;
	double threshold = sigma;//TODO
	int counter = 0;//number of marked location
	bool redOver, greenOver, blueOver;//booleans

	R2Image tempImg = *this;

	while (counter < numFeat){

		//loop thru image
		for (int i = 1; i < width - 1; i++){
			for (int j = 1; j < height - 1; j++){

				redOver = tempImg.Pixel(i, j).Red() >= threshold;
				greenOver = tempImg.Pixel(i, j).Green() >= threshold;
				blueOver = tempImg.Pixel(i, j).Blue() >= threshold;

				if (redOver && greenOver && blueOver){

					for (int x = -10; x <=10; x++){
						for (int y = -10; y <= 10; y++){

							//color red square on output image
							if ((x >=-1 && x <= 1) && (y >= -1&& y <= 1)){
								SetPixel(i + x, j + y, R2green_pixel);
							}

							//black out 10 by 10 area around red mark (where possible)
							if ((i + x >= 1) && (i + x < width - 1) && (j + y >= 1) && (j + y < height - 1)){
								tempImg.SetPixel(i + x, j + y, R2black_pixel);
							}

						}
					}
					counter++;
				}

				if (counter >= numFeat)
					break;

			}//inner-for

			if (counter >= numFeat)
				break;

		}//outer-for

		threshold -= 0.01;//decrement threshold to pick up more corners
	}

    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges
  

}


void R2Image::
Sharpen()
{
  // Sharpen an image using a linear filter. Use a kernel of your choosing.
       int k = 16;
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
        // using edge enhancement (e.g. 3x3 pixels with the following translations:
        //              [-k/8           -k/8            -k/8]
        //              [-k/8                   k                       -k/8]
        //              [-k/8                   -k/8            -k/8]
        //                                                                                                      )
 
        R2Image tempIMG(width, height);
        for (int i = 1; i < width-1; i++) {
            for (int j = 1; j < height-1; j++) {
 
                    double sumRed = Pixel(i, j).Red() * k;
                    double sumGreen = Pixel(i, j).Green() * k;
                    double sumBlue = Pixel(i, j).Blue() * k;
 
                    double redWeight = k;
                    double greenWeight = k;
                    double blueWeight = k;
 
                    int l[3] = { -1, 0, 1 };
                    for (int m = 0; m < 3; m++){
                        for (int n = 0; n < 3; n++) {
                            if (l[m] + i != i && l[n] + j != j) {
                                sumRed = sumRed + Pixel(l[m] + i, l[n] + j).Red() * (-k / 8);
                                sumBlue = sumBlue + Pixel(l[m] + i, l[n] + j).Blue() * (-k / 8);
                                sumGreen = sumGreen + Pixel(l[m] + i, l[n] + j).Green() * (-k / 8);
 
                                redWeight = redWeight + (-k / 8);
                                blueWeight = blueWeight + (-k / 8);
                                greenWeight = greenWeight + (-k / 8);
                            }
                        }
                    }
 
                    sumRed = sumRed / redWeight;
                    sumGreen = sumGreen / greenWeight;
                    sumBlue = sumBlue / blueWeight;
 
                    tempIMG.Pixel(i, j).SetRed(sumRed);
                    tempIMG.Pixel(i, j).SetBlue(sumBlue);
                    tempIMG.Pixel(i, j).SetGreen(sumGreen);
                }
        }
               
        for (int i = 1; i < width - 1; i++) {
            for (int j = 1; j < height - 1; j++) {
                            Pixel(i, j).SetRed(tempIMG.Pixel(i, j).Red());
                            Pixel(i, j).SetBlue(tempIMG.Pixel(i, j).Blue());
                            Pixel(i, j).SetGreen(tempIMG.Pixel(i, j).Green());
                            Pixel(i, j).Clamp();
            }
        }
 }
               
   



void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage" 
	// into this image with a 50% opacity.

	int numFeat = 150;
	//double sigma;
	double threshold = 2.0;

	Harris(1.0);
	otherImage->Harris(1.0);

	int numFeat2 = 150;
	double threshold2 = 1.0;//TODO
	int counter = 0;//number of marked location
	bool redOver, greenOver, blueOver;//booleans
	double r, g, b;
	double sRange = 7.0;


	R2Image tempImg = *this;
	tempImg.Harris(1.0);

	double difference = 0.0;
	double leastDifference;
	double xCord = 0;
	double yCord = 0;


	while (counter < numFeat){

		//loop thru image
		for (int i = 1; i < width - 1; i++){
			for (int j = 1; j < height - 1; j++){

				r = tempImg.Pixel(i, j).Red();
				g = tempImg.Pixel(i, j).Green();
				b = tempImg.Pixel(i, j).Blue();

				redOver = r >= threshold;
				greenOver = g >= threshold;
				blueOver = b >= threshold;

				if (redOver && greenOver && blueOver){
					double x;
					double y;
					for (x = -sRange; x <= sRange; x++){
						for (y = -sRange; y <= sRange; y++){

							difference += pow(r - otherImage->Pixel(i, j).Red(), 2);
							difference += pow(g - otherImage->Pixel(i, j).Green(), 2);
							difference += pow(b - otherImage->Pixel(i, j).Blue(), 2);

							leastDifference = 3.0;
							if (difference < leastDifference){
								leastDifference = difference;
								xCord = x + i;
								yCord = y + j;
							}

							//color red square on output image
							if ((x >= -1 && x <= 1) && (y >= -1 && y <= 1)){
								SetPixel(i + x, j + y, R2green_pixel);
								//otherImage->line(xCord, x, yCord, y, 1.0, 0.0, 1.0);
							}

							//black out 10 by 10 area around red mark (where possible)
							if ((i + x >= 1) && (i + x < width - 1) && (j + y >= 1) && (j + y < height - 1)){
								tempImg.SetPixel(i + x, j + y, R2black_pixel);


							}


						}
					}
					
					///("%d", line(xCord, x, yCord, y, 1.0, 0.0, 1.0));
					otherImage->line(xCord, x, yCord, y, 1.0, 0.0, 0.0);
					

					counter++;
				}

				if (counter >= numFeat)
					break;

			}//inner-for

			if (counter >= numFeat)
				break;

		}//outer-for

		threshold -= 0.01;//decrement threshold to pick up more corners
	}

}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
	fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
	return;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


	

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 75, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
	
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






