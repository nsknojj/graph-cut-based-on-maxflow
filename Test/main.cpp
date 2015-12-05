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
#include "maxflow.h"
#include "threSeg.h"

using namespace std;
using namespace cv;

const int maxn = 1510, maxm = 1510;

int R[maxn][maxn], G[maxn][maxm], B[maxn][maxm];
double a[maxn][maxm], b[maxn][maxm], p[maxn][maxm][2];
int n, m;
extern int S, T;
extern bool f[];
const float INF = 100000000;

// colorReduce: an example of dealing image 
void colorReduce(Mat &image, int div = 64) {
	int nr = image.rows; // number of rows
	int nc = image.cols * image.channels(); // total number of elements per line
	cout << image.channels() << endl;
	for (int j = 0; j<nr; j++) {
		uchar* data = image.ptr<uchar>(j);
		for (int i = 0; i<nc; i++) {
			*data = *data / div*div + div / 2;
			if (i % 3 != 2) *data = 0;		// 0 blue 1 green 2 red
			data++;
		} // end of row                 
	}
}
// get RGB mat
void getRGB(Mat &image) {
	int nr = image.rows; // number of rows
	int nc = image.cols * image.channels(); // total number of elements per line
	//cout << image.channels() << endl;
	n = nr;
	m = nc / 3;
	for (int j = 0; j<nr; j++) {
		uchar* data = image.ptr<uchar>(j);
		for (int i = 0; i<nc; i++) {
			if (i % 3 == 0) B[j][i / 3] = *data;
			if (i % 3 == 1) G[j][i / 3] = *data;
			if (i % 3 == 2) R[j][i / 3] = *data;
			data++;
		} // end of row                 
	}
}

float getdiff(int i1, int j1, int i2, int j2){
	int R1, G1, B1, R2, G2, B2;
	float Y1, U1, V1, Y2, U2, V2;
	R1 = R[i1][j1];
	G1 = G[i1][j1];
	B1 = B[i1][j1];
	R2 = R[i2][j2];
	G2 = G[i2][j2];
	B2 = B[i2][j2];
	Y1 = R1;
	U1 = G1;
	V1 = B1;
	Y2 = R2;
	U2 = G2;
	V2 = B2;
	int u1 = 0.299*R1 + 0.587*G1 + 0.114*B1;
	int u2 = 0.299*R2 + 0.587*G2 + 0.114*B2;
	Y1 -= u1;
	U1 -= u1;
	V1 -= u1;
	Y2 -= u2;
	U2 -= u2;
	V2 -= u2;

	/*
	Y1 = 0.299*R1 + 0.587*G1 + 0.114*B1;
	U1 = 0.436*(B1 - Y1) / (1 - 0.114) + 128;
	V1 = 0.615*(R1 - Y1) / (1 - 0.299) + 128;
	Y2 = 0.299*R2 + 0.587*G2 + 0.114*B2;
	U2 = 0.436*(B2 - Y2) / (1 - 0.114) + 128;
	V2 = 0.615*(R2 - Y2) / (1 - 0.299) + 128;
	*/
	return pow(((Y1 - Y2)*(Y1 - Y2) + (U1 - U2)*(U1 - U2) + (V1 - V2)*(V1 - V2)), 0.5) + 1.5 * abs(u1 - u2);
}

bool cmp(uchar a, uchar b){
	return abs(a - b)>10;
}

int d[maxn][maxm], from[maxn][maxm][2];
int q[maxn*maxm][2];
const int D[4][2] = { {0,1},{1,0},{-1,0},{0,-1} };
double avebgcolor[3], avefgcolor[3];
int sum1, sum2;

float sqr(float x){ return x*x; }

float colordis(double R1, double G1, double B1, double R2, double G2, double B2) {
	double Y1, U1, V1, Y2, U2, V2;
	Y1 = R1;
	U1 = G1;
	V1 = B1;
	Y2 = R2;
	U2 = G2;
	V2 = B2;
	/*
	Y1 = 0.299*R1 + 0.587*G1 + 0.114*B1;
	U1 = 0.436*(B1 - Y1) / (1 - 0.114) + 128;
	V1 = 0.615*(R1 - Y1) / (1 - 0.299) + 128;
	Y2 = 0.299*R2 + 0.587*G2 + 0.114*B2;
	U2 = 0.436*(B2 - Y2) / (1 - 0.114) + 128;
	V2 = 0.615*(R2 - Y2) / (1 - 0.299) + 128;
	*/
	int u1 = 0.299*R1 + 0.587*G1 + 0.114*B1;
	int u2 = 0.299*R2 + 0.587*G2 + 0.114*B2;
	Y1 -= u1;
	U1 -= u1;
	V1 -= u1;
	Y2 -= u2;
	U2 -= u2;
	V2 -= u2;
	return pow(((Y1 - Y2)*(Y1 - Y2) + (U1 - U2)*(U1 - U2) + (V1 - V2)*(V1 - V2)), 0.5) + 1.5 * abs(u1 - u2);
}

void calcAB(Mat image, Mat image2, Mat image3) {
	int nr = image.rows; // number of rows
	int nc = image.cols * image.channels(); // total number of elements per line
	int l, r;
	l = 1; r = 0;
	int maxd = 0;
	sum1 = sum2 = 0;
	memset(d, 0, sizeof d);
	for (int j = 0; j<nr; j++) {
		uchar* data = image.ptr<uchar>(j);
		uchar* data2 = image2.ptr<uchar>(j);
		uchar* data3 = image3.ptr<uchar>(j);
		for (int i = 0; i<nc; i += 3) {
			if (cmp(*data, *data3) || cmp(*(data + 1), *(data3 + 1)) || cmp(*(data + 2), *(data3 + 2))) {
				avebgcolor[0] += *data;
				avebgcolor[1] += *(data + 1);
				avebgcolor[2] += *(data + 2);
				sum2++;
				b[j][i/3] = INF;
				a[j][i/3] = 0;
			}
			if (cmp(*data,*data2) || cmp(*(data + 1),*(data2 + 1)) || cmp(*(data + 2),*(data2 + 2))) {
				b[j][i/3] = 0;
				a[j][i/3] = INF;
				d[j][i / 3] = 1;
				avefgcolor[0] += *data;
				avefgcolor[1] += *(data+1);
				avefgcolor[2] += *(data + 2);
				sum1++;
				q[++r][0] = j;
				q[r][1] = i/3;
				from[j][i / 3][0] = j;
				from[j][i / 3][1] = i / 3;
			}
			data = data + 3;
			data2 = data2 + 3;
			data3 = data3 + 3;
		}              
	}
	avefgcolor[0] /= sum1;
	avefgcolor[1] /= sum1;
	avefgcolor[2] /= sum1;
	avebgcolor[0] /= sum2;
	avebgcolor[1] /= sum2;
	avebgcolor[2] /= sum2;
	
	for (; l <= r; l++){
		int i = q[l][0], j = q[l][1];
		maxd = max(d[i][j], maxd);
		for (int dd = 0; dd < 4; dd++) {
			int tx = i + D[dd][0], ty = j + D[dd][1];
			if (tx >= 0 && tx < n&&ty >= 0 && ty < m) {
				if (!d[tx][ty]) {
					from[tx][ty][0] = from[i][j][0];
					from[tx][ty][1] = from[i][j][1];
					d[tx][ty] = d[i][j] + 1;
					q[++r][0] = tx;
					q[r][1] = ty;
				}
			}
		}
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++){
			if (a[i][j] >= INF-1 || b[i][j] >= INF-1){}
			else {
				float db = (colordis(R[i][j], G[i][j], B[i][j], avebgcolor[0], avebgcolor[1], avebgcolor[2]));
				float df = (colordis(R[i][j], G[i][j], B[i][j], avefgcolor[0], avefgcolor[1], avefgcolor[2]));
				a[i][j] = db / (db + df);
				b[i][j] = df / (db + df);
			}
		}

}

void calcPij() {
	for (int i = 0; i < n;i++)
	for (int j = 0; j < m; j++) {
		double lamda = 100.0, eps = 0.01;
		p[i][j][0] = lamda/(eps+(getdiff(i, j, i + 1, j)));
		p[i][j][1] = lamda/(eps+(getdiff(i, j, i, j + 1)));
	}
}

int getnum(int i, int j) {
	return i * m + j + 3;
}

int geti(int i) {
	return (i - 3) / m;
}

int getj(int i){
	return (i - 3) % m;
}

/*
	build maxflow graph
	a[i] foreground probability of pixel
	b[i] background probability
	p[i][j] penalty of i,j belonging to different sets
	S T
	edge:
		S->i, a[i]
		i->T, b[i]
		i->j, p[i][j]
*/
void build() {
	S = 1;
	T = 2;
	for (int i = 0; i < n;i++)
		for (int j = 0; j < m; j++){
			addedge(S, getnum(i, j), a[i][j]);
			addedge(getnum(i, j), S, 0);
			addedge(getnum(i, j), T, b[i][j]);
			addedge(T, getnum(i, j), 0);
			if (i < n) {
				addedge(getnum(i, j), getnum(i + 1, j), p[i][j][0]);
				addedge(getnum(i + 1, j), getnum(i, j), p[i][j][0]);
			}
			if (j < m) {
				addedge(getnum(i, j), getnum(i, j + 1), p[i][j][1]);
				addedge(getnum(i, j + 1), getnum(i, j), p[i][j][1]);
			}
		}
}

void setimage(Mat image) {
	int nr = image.rows; // number of rows
	int nc = image.cols * image.channels(); // total number of elements per line
	//cout << image.channels() << endl;
	for (int j = 0; j<nr; j++) {
		uchar* data = image.ptr<uchar>(j);
		for (int i = 0; i<nc; i++) {
			if (!f[getnum(j, i / 3)]) *data = 0;
			data++;
		}
	}
}

void addimg(Mat &image, Mat image2){
	int nr = image.rows; // number of rows
	int nc = image.cols * image.channels(); // total number of elements per line
	for (int j = 0; j<nr; j++) {
		uchar* data = image.ptr<uchar>(j);
		uchar* data2 = image2.ptr<uchar>(j);
		for (int i = 0; i<nc; i++) {
			*data = min((int)*data2+(int)*data, 255);
			data++;
			data2++;
		}
	}
}

int main()
{
	const char* imagename = "..\\Images\\in.bmp";
	const char* imagename2 = "..\\Images\\fg.bmp";
	const char* imagename3 = "..\\Images\\bg.bmp";

	//从文件中读入图像
	Mat img = imread(imagename);
	Mat img2 = imread(imagename2);
	Mat img3 = imread(imagename3);
	//Mat img3 = threSeg(img);
	//Mat img4;
	//cvtColor(img3, img4, CV_GRAY2RGB);
	//addimg(img4, img);
	//imshow("image", img);
	//imshow("image2",img4);
	getRGB(img);
	calcAB(img, img2, img3);
	calcPij();
	build();
	getflow();
	cut();
	setimage(img);
	blur(img, img, Size(2, 2));
	imshow("image", img);
	cvWaitKey();
	return 0;
}