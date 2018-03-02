// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef ESTIMATE_GMM_H
#define ESTIMATE_GMM_H


#include <iostream>
#include <cmath>
#include <vector>

//edda includes
#include <core/interpolator.h>
#include <distributions/gaussian_mixture.h>
#include <distributions/gmm.h>

namespace edda
{

///////////////////////////////////////////////////////////////////////////////////////////
///// Traditional EM algorithm section

///////////////////////////////////////////
/// \brief EM helper function
///
inline int check_convergence(double fmax, double fmin, double ftol)
{
    double EPS = 1e-10;
    double delta = fabs(fmax - fmin);
    double accuracy = (fabs(fmax) + fabs(fmin)) * ftol;
    return (delta < (accuracy + EPS));
}

///////////////////////////////////////////
/// \brief Evaluate probability of a value from a Gaussian
/// \mean mean of a Gaussian
/// \sigma sigma of a Gaussian
/// \val given sample value
inline double eval_gaussian_density(double mean, double sigma, double val)
{
    double EPS = 1e-10;
    double prob=0.0;
    double exponent=0.0;

    if(sigma>0)
    {
        exponent = (val-mean)/sigma;
        exponent = -(exponent*exponent)/2.0;
        prob =  ( 1.0 / (sigma*sqrt( 2.0*M_PI))) * exp(exponent);
    }
    else
    {
        if(fabs(val-mean) < EPS)
            prob=1.0;
        else
            prob=0.0;
    }

    return prob;
}

///////////////////////////////////////////
/// \brief Upate parameter of a Gaussian
inline void update_parameters(int n, double * data, int k, double * prob, double * mean, double * sd, double ** class_prob)
{
    double EPS = 1e-6;

    //Update weights first
    for (int j = 0; j < k; j++)
    {
        prob[j] = 0.0;
        for (int i = 0; i < n; i++)
            prob[j] += class_prob[i][j];

        prob[j] /= n;
    }

    //update mean
    for (int j = 0; j < k; j++)
    {
        mean[j] = 0.0;
        for (int i = 0; i < n; i++)
            mean[j] += data[i] * class_prob[i][j];

        mean[j] /= n * prob[j] + EPS;
    }

    //update standard deviation
    for (int j = 0; j < k; j++)
    {
        sd[j] = 0.0;
        for (int i = 0; i < n; i++)
            sd[j] += (data[i] - mean[j])*(data[i] - mean[j]) * class_prob[i][j];
        sd[j] /= (n * prob[j] + EPS);
        sd[j] = sqrt(sd[j]);
    }
}

///////////////////////////////////////////
/// EM helper function
inline double classprob(int j, double x, int k, double *prob, double *mean, double *sd)
{
    double num = prob[j]*eval_gaussian_density(mean[j],sd[j],x);

    double denom=0;
    for(int i=0;i<k;i++)
        denom += prob[i]*eval_gaussian_density(mean[i],sd[i],x);

    return num/denom;
}

///////////////////////////////////////////
/// EM helper function
inline void update_class_prob(int n, double * data, int k, double * prob, double * mean, double * sd, double ** class_prob)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
            class_prob[i][j] = classprob(j, data[i],k, prob, mean, sd);
}

///////////////////////////////////////////
/// EM helper function
inline double computeLogLikelihood(int n, double* data,int k, double* prob,double* mean,double* sd)
{
    double llk=0;

    for(int p=0;p<n;p++)
    {
        double val=0;
        for(int q=0;q<k;q++)
        {
            val += prob[q]*eval_gaussian_density(mean[q],sd[q],data[p]);
        }

        llk += log(val);
    }

    return llk/n;
}


//////////////////////////////////////////////
/// Main EM function
///
inline dist::GMM eddaComputeGMM(double *dataArray, int nSamples, int nComps)
{
	double eps = 1e-6;
    int GMMs = nComps;

	if (nSamples<GMMs)
	{
		std::cout << "Number of samples must be larger than number of clusters..." << std::endl;
		dist::GMM new_gmm = dist::GMM(nComps);
		return new_gmm;
	}
	else
	{
		double llk = 0, prev_llk = 0;
		double *mean, *sd, *weight;
		double **class_prob;

		//Allocate memories for computation
		/////////////////////////////////////////////////////////////
		class_prob = (double **)malloc(sizeof(double *)*nSamples);
		for (int i = 0; i<nSamples; i++)
			class_prob[i] = (double *)malloc(sizeof(double)* GMMs);

		mean = (double *)malloc(sizeof(double)* GMMs);
		sd = (double *)malloc(sizeof(double)* GMMs);
		weight = (double *)malloc(sizeof(double)* GMMs);

		//initial estimate of parameters
		/////////////////////////////////////////////
		double mean1 = 0.0, sd1 = 0.0;

		for (int i = 0; i < nSamples; i++)
			mean1 += dataArray[i];
		mean1 /= nSamples;

		for (int i = 0; i < nSamples; i++)
			sd1 += (dataArray[i] - mean1)*(dataArray[i] - mean1);
		sd1 = sqrt(sd1 / nSamples);

		for (int j = 0; j < GMMs; j++)
		{
			weight[j] = 1.0 / GMMs;
			mean[j] = dataArray[rand() % nSamples];
			sd[j] = sd1;
		}

		//Do while loop for iterative estimation (limit the maximum interations as 100)
        /////////////////////////////////////////////////
        int iters = 0;
		do{
			//save prev likelihood
			prev_llk = llk;

			//update probabilities
			update_class_prob(nSamples, dataArray, GMMs, weight, mean, sd, class_prob);

			//update the parameters with newly estimated probabilities
			update_parameters(nSamples, dataArray, GMMs, weight, mean, sd, class_prob);

			//compute new likelihood
            llk = computeLogLikelihood(nSamples, dataArray, GMMs, weight, mean, sd);
            
            iters ++;
		} while (!check_convergence(llk, prev_llk, eps) && iters <100);

		//Update the gmm object with estimated values
		/////////////////////////////////////////////////
		dist::GMM new_gmm = dist::GMM(nComps);
		for (int m = 0; m<GMMs; m++)
		{
			new_gmm.models[m].m = mean[m];
			new_gmm.models[m].v = sd[m];
			new_gmm.models[m].w = weight[m];
		}

		//Clean up
		free(mean);
		free(sd);
		free(weight);

		for (int i = 0; i < nSamples; i++)
			free(class_prob[i]);
		free(class_prob);

		return new_gmm;
	}
}

//////////////////////////////////////////////
/// Main EM function
///
template <int GMMs>
void eddaComputeEM(double *samples, int numSamples, dist::GaussianMixture<GMMs>* new_gmm)
{
    double eps =  1e-6;

    if(numSamples<GMMs)
    {
        std::cout<<"Number of samples must be larger than number of clusters..."<<std::endl;
        return;
    }
    else
    {
        double llk = 0, prev_llk = 0;
        double *mean,*sd,*weight;
        double **class_prob;
        
        //Allocate memories for computation
        /////////////////////////////////////////////////////////////
        class_prob = (double **) malloc(sizeof(double *)*numSamples);
        for(int i=0;i<numSamples;i++)
            class_prob[i] = (double *) malloc(sizeof(double) * GMMs);

        mean = (double *) malloc(sizeof(double) * GMMs);
        sd = (double *) malloc(sizeof(double) * GMMs);
        weight = (double *) malloc(sizeof(double) * GMMs);

        //initial estimate of parameters
        /////////////////////////////////////////////
        double mean1 = 0.0, sd1 = 0.0;
        
        for (int i = 0; i < numSamples; i++)
            mean1 += samples[i];
        mean1 /= numSamples;
        
        for (int i = 0; i < numSamples; i++)
            sd1 += (samples[i] - mean1)*(samples[i] - mean1);
        sd1 = sqrt(sd1 / numSamples);
        
        for (int j = 0; j < GMMs; j++)
        {
            weight[j] = 1.0 / GMMs;
            mean[j] = samples[rand() % numSamples];
            sd[j] = sd1;
        }
        
        //Do while loop for iterative estimation (limit the maximum interations as 100)
        /////////////////////////////////////////////////
        int iters = 0;
        do{
            //save prev likelihood
            prev_llk = llk;

            //update probabilities
            update_class_prob(numSamples, samples, GMMs, weight, mean, sd, class_prob);

            //update the parameters with newly estimated probabilities
            update_parameters(numSamples, samples, GMMs, weight, mean, sd, class_prob);

            //compute new likelihood
            llk = computeLogLikelihood(numSamples, samples, GMMs, weight, mean, sd);

            iters ++;
        } while (!check_convergence(llk, prev_llk, eps) && iters < 100 );
        
        //Update the gmm object with estimated values
        /////////////////////////////////////////////////
        for(int m=0;m<GMMs;m++)
        {
            new_gmm->models[m].m = mean[m];
            new_gmm->models[m].v = sd[m];
            new_gmm->models[m].w = weight[m];
        }

        //Clean up
        free(mean);
        free(sd);
        free(weight);

        for (int i = 0; i < numSamples; i++)
            free(class_prob[i]);
        free(class_prob);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
///// Incremental section

///////////////////////////////////////////
/// Incremental update helper function
template <int GMMs>
bool testforMatch(double val, dist::GaussianMixture<GMMs>* currentGMM, std::vector<float> *valList, std::vector<float> *distList)
{
    float sigmaTh = 2.0;
    int match=0;
    valList->clear();
    distList->clear();

    for(int i=0;i<GMMs;i++)
    {
        float threshold[2];
        threshold[0] = currentGMM->models[i].m - sigmaTh*sqrt(currentGMM->models[i].v);
        threshold[1] = currentGMM->models[i].m + sigmaTh*sqrt(currentGMM->models[i].v);

        distList->push_back(fabs(currentGMM->models[i].m - val));

        if((val>= threshold[0]) && (val<= threshold[1]))
        {
            valList->push_back(fabs(currentGMM->models[i].m - val));
            match++;
        }
        else
        {
            valList->push_back(9999999);
        }
    }

    return (match);
}

//////////////////////////////////////////////
/// Main Incremental GMM update function
///
template <int GMMs>
void eddaUpdateGMMIncremental(double *samples, int numSamples, dist::GaussianMixture<GMMs>* new_gmm)
{
    float alpha = 0.2;//Mixing rate in updation

    float initVar = new_gmm->models[0].v;
    for(int i=0;i<GMMs;i++)
    {
        if(new_gmm->models[i].v>initVar)
            initVar =new_gmm->models[i].v;
    }

    for(int i=0;i<numSamples;i++)
    {
        std::vector<float> valList;
        std::vector<float> distList;
        int minIndex=0;
        float minProb=0.0;
        bool match;
        float sumWeight=0.0;

        match = testforMatch(samples[i], new_gmm,&valList, &distList);

        //Find the max prob Gaussian from the current GMM
        minProb = valList[0]; minIndex=0;
        for(int k=0;k<valList.size();k++)
        {
            if(minProb>valList[k])
            {
                minProb = valList[k];
                minIndex = k;
            }
        }

        //updation for the matched gaussian
        if(match)
        {
            //Update gaussians accordingly
            for(int k=0;k<valList.size();k++)
            {
                if(k==minIndex)
                {
                    //Update weight
                    new_gmm->models[k].w += alpha*(1-new_gmm->models[k].w);

                    //Update Mean
                    new_gmm->models[k].m = (1-alpha)*new_gmm->models[k].m + alpha*samples[i];

                    //Update standard Deviation
                    float dev = new_gmm->models[k].v;
                    float temp = (samples[i] - new_gmm->models[k].m);
                    dev = (1-alpha)*dev + alpha*temp*temp;

                    new_gmm->models[k].v = dev;
                }
                //updation for the other gaussians
                else
                {
                    //Update weight
                    new_gmm->models[k].w = (1-alpha)*new_gmm->models[k].w;
                }
            }

        }
        //Update if no match found
        else
        {
            int mid = 0;
            if(valList.size() > GMMs-1)
            {
                //find gaussian with min weight
                float minw = new_gmm->models[0].w;
                for(int k=0;k<valList.size();k++)
                {
                    if(new_gmm->models[k].w < minw)
                    {
                        mid=k;
                        minw = new_gmm->models[k].w;
                    }
                }
            }

            new_gmm->models[mid].m = samples[i];
            new_gmm->models[mid].v = initVar;
            new_gmm->models[mid].w = 0.00001;
        }

        //Renormalize weights
        sumWeight=0.0;
        for(int k=0;k<valList.size();k++)
            sumWeight += new_gmm->models[k].w;

        for(int k=0;k<valList.size();k++)
            new_gmm->models[k].w = new_gmm->models[k].w/sumWeight;

        valList.clear();
    }
}

}

#endif  // ESTIMATE_GMM
