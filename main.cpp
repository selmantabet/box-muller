/* Selman Tabet (@selmantabet - https://selman.io/) - UIN 724009859
ECEN-449 Project - Part 2: Pseudo Random Number Generator (PRNG).

This program generates samples from a Gaussian distribution that is based on
two separate uniform distributionswith specified mean (mu) and standard
deviation (sigma). The program takes 2000 samples at a rate of 10Hz, and the
mean and standard deviation would be recomputed based on all samples collected.
The results would then be printed through the Serial port.
The error margin between the specified mean and the current calculated mean
are calculated and displayed as an integer value on the SPI 
interface (SPI LEDs on the board).

Developed using the Mbed IDE. Tested on an EA LPC4088 QuickStart Board. */

#include "mbed.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include <time.h>

Ticker tick; //Define ticker object
Serial pc(USBTX, USBRX); //Serial channel over HDK USB interface
SPI shifter(p5, p6, p7); //MOSI, MISO, SCLK
DigitalOut chip_select(p30); //Chip select

const int sample_size = 2000; //Can be dynamically adjusted per spec.
double samples[sample_size]; //Array of samples.

int counter = 0; //Keeps count for sigma calculations and for CLI print.
double sample; //Single sample.
double mean = 0; //Current mean value.
double std_deviation = 0; //Sigma (Standard Deviation AKA Variance^0.5)
int margin; //Integer percentage margin of error of the mean.

double rand_normal(double mu, double sigma) { //Box-Muller PRNG.

    /*Mersenne Twister engine (mt19937) is a more widely used PRNG (for good
    reasons) but the handout forces us to use rand() and Box-Muller. So here
    it is, for the sake of fulfilling the specified requirements.
    
    The implementation is based on the R. Knop revision of the method
    shown on the Wikipedia article linked in the handout. The z0 and z1
    formulae used here are available under the "Polar form" section
    of the Box-Muller Wikipedia article.*/
    
    static double z_1 = 0.0;
    static int z_1_cached = 0; //For temporarily storing z1 between calls.
    if (!z_1_cached)
    {
        double u, v, s; // u : x-component v : y-component  s : R^2 = u^2 + v^2
        do
        {
            u = 2.0*rand()/RAND_MAX - 1; //Uniform distribution 1.
            v = 2.0*rand()/RAND_MAX - 1; //Uniform distribution 2.

            s = u*u + v*v; //Polar form set expression.
        }
        while (s == 0.0 || s > 1.0); //Re-gen u,v pair if s-condition is false
        {
            double d = sqrt(-2.0*log(s)/s); //B-M intermediate.
            double z_0 = u*d;
            z_1 = v*d;
            double result = z_0*sigma + mu; //Scaling to defined parameters.
            z_1_cached = 1; //Use z1 on next call.
            return result;
        }
    }
    else
    {
        z_1_cached = 0; //Discard for pair regeneration on next call.
        return z_1*sigma + mu;
    }
}

void generate(){
    if (counter >= sample_size){
        pc.printf("Final mean: %f         Final standard deviation: %f\n",
        mean, std_deviation);
        exit(0);   
    }
    else {
        sample = rand_normal(10.0, 1.0); //Mean = 10, Sigma = 1.
        samples[counter] = sample;
        //Recompute mean and sigma.
        mean = ((mean * counter) + sample)/(counter + 1);
        double tmp = 0.0; //Squared-diff intermediate.
        for (int i = 0; i <= counter; i++){
            tmp += pow(abs(samples[i] - mean), 2.0); //Squared-diff summation.
        }
        counter++;
        std_deviation = sqrt(tmp / counter); //Final sigma calculation step.
        
        
        pc.printf("Sample #%i: %f    ( Mean: %f    Standard deviation: %f )\n",
        counter, sample, mean, std_deviation);
        
        margin = (int)((abs(10 - mean)*10) + 0.5); //Rounded via typecast.

        chip_select = 0; shifter.write(margin); chip_select = 1;
        
    }
}

int main(){
    chip_select = 0; //Select the device by setting chip select low
    shifter.write(0); //Clear LEDs
    chip_select = 1; //Deselect the device
    srand(time(NULL)); //Seed PRNG.
    tick.attach(&generate, 0.1); //Sampling rate set at 10Hz.
    while(1); //Run forever.
}