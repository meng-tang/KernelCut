#ifndef _GPUFILTERING_H_
#define _GPUFILTERING_H_
#include <iostream>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <ctime>

#define BLOCK_SIZE 768
using namespace std;

int gpufilteringdemo(int N);
float * gpufiltering(float * X, float * values, float * sigmas, size_t N, char flag);
#endif
