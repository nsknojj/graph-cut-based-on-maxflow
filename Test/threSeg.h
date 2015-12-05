#ifndef THRESEG_H
#define THRESEG_H
#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/legacy/legacy.hpp>
#include <vector>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace cv;

Mat threSeg(const Mat &Image); //全局阈值分割

#endif // THRESEG_H
