#include <cstdlib>
#include <ctime>
#include <iostream>

#include "distributions/distribution.h"
#include "distributions/variant.h"
#include "dataset/distr_array.h"
#include "distributions/joint_histogram.h"
#include "bmp_image.h"
#include <time.h>

using namespace edda;
using namespace std;
using namespace edda::dist;

//#define DEFAULT_TEST

DistrArray * make_JointHistogram_array() {
	shared_ary<JointHistogram> array(new JointHistogram[10], 10);
	DistrArray * abstract_array = new JointDistrArray<JointHistogram>(array);
	return abstract_array;
}

int main(int argc, char* argv[]) {
#ifdef DEFAULT_TEST
	// test the default constructor
	DistrArray * array = make_JointHistogram_array();
	int i;
	cout << "array: " << endl;
	for (i=0; i<10; i++) {
		vector<Real> v = array->getVector(i);
		cout << i << ": " << array->getDistr(i) <<
			": sample = " << v[0] << " " << v[1] << " " << v[2] << endl;
	}
#else
	// test with the image processing example
	if(argc!=2) {
		printf("[usage]:./test_jointHistogram [bmp image name]\n");
		return 0;
	}
	BMPImage image(argv[1]);
	BMPImage image1(argv[1]);
	BMPImage image2(argv[1]);
	BMPImage image3(argv[1]);
	BMPImage image4(argv[1]);
	BMPImage image5(argv[1]);
	BMPImage image6(argv[1]);
		
	int blockSize = 20;	// divide image into these many blocks
	int dsUs = image.width/blockSize;
	int dsVs = image.height/blockSize;
	shared_ary<JointHistogram> array(new JointHistogram[dsUs*dsVs], dsUs*dsVs);
	
	// test on marginalization 
	std::unordered_set<int> marVars;
	marVars.insert(1);// marginalized dimension
	marVars.insert(2);// marginalized dimension
	
	// test on conditional histogram
	std::unordered_set<int> vars;
	vars.insert(0); vars.insert(1);
	std::vector<int> condVars;
	condVars.push_back(2);
	std::vector<std::pair<int, int>> bin_range;
	std::pair<int, int> rg(0, 29);
	bin_range.push_back(rg);

	clock_t t1, t2;
	float elastTime = 0;
	// go through each block to get the joint histogram
	for (int dsV = 0; dsV < dsVs; dsV++){//row
		for (int dsU = 0; dsU < dsUs; dsU++){//col
			//printf("[test_jointHistogram debug] working on block: %d %d\n", dsV, dsU);
			// prepare raw data
			Real* var1 = (Real*)malloc(sizeof(Real)*blockSize*blockSize);
			Real* var2 = (Real*)malloc(sizeof(Real)*blockSize*blockSize);
			Real* var3 = (Real*)malloc(sizeof(Real)*blockSize*blockSize);
			std::vector<Real> min_val(3, FLT_MAX);
			std::vector<Real> max_val(3, -FLT_MAX);
			std::vector<int> bins(3, 30);
			int cnt = 0;
			for(int v=dsV*blockSize; v<(dsV+1)*blockSize; v++) {//row
				for(int u=dsU*blockSize; u<(dsU+1)*blockSize; u++) {//col
					var1[cnt] = (Real)(image.bitmapImage[(v*image.width+u)*3+0]); 
					var2[cnt] = (Real)(image.bitmapImage[(v*image.width+u)*3+1]); 
					var3[cnt] = (Real)(image.bitmapImage[(v*image.width+u)*3+2]); 
					if(var1[cnt]<min_val[0]) min_val[0] = var1[cnt];
					if(var1[cnt]>max_val[0]) max_val[0] = var1[cnt];
					if(var2[cnt]<min_val[1]) min_val[1] = var2[cnt];
					if(var2[cnt]>max_val[1]) max_val[1] = var2[cnt];
					if(var3[cnt]<min_val[2]) min_val[2] = var3[cnt];
					if(var3[cnt]>max_val[2]) max_val[2] = var3[cnt];
					cnt++;
				}
			}
			std::vector<Real*> dataAry;
			dataAry.push_back(var1);
			dataAry.push_back(var2);
			dataAry.push_back(var3);
			array[dsV*dsUs+dsU] = eddaComputeJointHistogram(dataAry, cnt, min_val, max_val, bins);
			// safe to free data, after constructing the distribution
			free(var1);
			free(var2);
			free(var3);
								
			// test getJointMean
			for(int v=dsV*blockSize; v<(dsV+1)*blockSize; v++) {//row
				for(int u=dsU*blockSize; u<(dsU+1)*blockSize; u++) {//col
					std::vector<Real> mean = getJointMean(array[dsV*dsUs+dsU]); 
					image1.bitmapImage[(v*image1.width+u)*3+0] = mean[0]; 
					image1.bitmapImage[(v*image1.width+u)*3+1] = mean[1]; 
					image1.bitmapImage[(v*image1.width+u)*3+2] = mean[2]; 
				}
			}
			
			// test getJointSample
			t1 = clock();
			for(int v=dsV*blockSize; v<(dsV+1)*blockSize; v++) {//row
				for(int u=dsU*blockSize; u<(dsU+1)*blockSize; u++) {//col
					std::vector<Real> sample = getJointSample(array[dsV*dsUs+dsU]); 
					image2.bitmapImage[(v*image2.width+u)*3+0] = sample[0]; 
					image2.bitmapImage[(v*image2.width+u)*3+1] = sample[1]; 
					image2.bitmapImage[(v*image2.width+u)*3+2] = sample[2]; 
				}
			}
			t2 = clock();
			elastTime += ((float)t2-(float)t1)/CLOCKS_PER_SEC;

			// test marginalization
			JointHistogram hist = array[dsV*dsUs+dsU].marginalization(marVars);
			for(int v=dsV*blockSize; v<(dsV+1)*blockSize; v++) {//row
				for(int u=dsU*blockSize; u<(dsU+1)*blockSize; u++) {//col
					std::vector<Real> mean = getJointMean(hist); 
					// the order below decides by the elements order in the marVars
					image3.bitmapImage[(v*image3.width+u)*3+0] = 0;
					image3.bitmapImage[(v*image3.width+u)*3+1] = mean[1]; 
					image3.bitmapImage[(v*image3.width+u)*3+2] = mean[0]; 
				}
			}
			for(int v=dsV*blockSize; v<(dsV+1)*blockSize; v++) {//row
				for(int u=dsU*blockSize; u<(dsU+1)*blockSize; u++) {//col
					std::vector<Real> sample = getJointSample(hist); 
					// the order below decides by the elements order in the marVars
					image4.bitmapImage[(v*image4.width+u)*3+0] = 0; 
					image4.bitmapImage[(v*image4.width+u)*3+1] = sample[1]; 
					image4.bitmapImage[(v*image4.width+u)*3+2] = sample[0]; 
				}
			}

			// test conditionalHistogram
			JointHistogram chist = array[dsV*dsUs+dsU].conditionalHist(vars, condVars, bin_range);
			for(int v=dsV*blockSize; v<(dsV+1)*blockSize; v++) {//row
				for(int u=dsU*blockSize; u<(dsU+1)*blockSize; u++) {//col
					std::vector<Real> mean = getJointMean(chist); 
					image5.bitmapImage[(v*image5.width+u)*3+0] = mean[0];
					image5.bitmapImage[(v*image5.width+u)*3+1] = mean[1]; 
					image5.bitmapImage[(v*image5.width+u)*3+2] = 0; 
				}
			}
			for(int v=dsV*blockSize; v<(dsV+1)*blockSize; v++) {//row
				for(int u=dsU*blockSize; u<(dsU+1)*blockSize; u++) {//col
					std::vector<Real> sample = getJointSample(chist); 
					image6.bitmapImage[(v*image6.width+u)*3+0] = sample[0]; 
					image6.bitmapImage[(v*image6.width+u)*3+1] = sample[1]; 
					image6.bitmapImage[(v*image6.width+u)*3+2] = 0; 
				}
			}
		}
	}
	// write back bmp images and print statistics on screen
	printf("[test_jointHistogram] total time for samping: %f\n", elastTime);
	printf("[test_jointHistogram] total blocks: %d, %d\n", dsVs, dsUs);
	printf("[test_jointHistogram] time for samping each block: %f\n", elastTime/(dsVs*dsUs));
	image1.writeImage(std::string("mean.bmp").c_str());
	image2.writeImage(std::string("dist.bmp").c_str());
	image3.writeImage(std::string("marMean.bmp").c_str());
	image4.writeImage(std::string("marDist.bmp").c_str());
	image5.writeImage(std::string("condMean.bmp").c_str());
	image6.writeImage(std::string("condDist.bmp").c_str());
	printf("[test_jointHistogram] results of getJointMean() and getJointSample() using original joint histogram are in mean.bmp and dist.bmp\n");
	printf("[test_jointHistogram] results of getJointMean() and getJointSample() using marginalized joint histogram are in marMean.bmp and marDist.bmp\n");
	printf("[test_jointHistogram] results of getJointMean() and getJointSample() using conditional joint histogram are in condMean.bmp and condDist.bmp\n");
#endif
}
