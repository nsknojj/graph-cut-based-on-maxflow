#include "threSeg.h"

Mat threSeg(const Mat &Image)
{
	Mat ans;
	cvtColor(Image, ans, CV_RGB2GRAY);
	int *pixelNum = (int*)calloc(256, sizeof(int));                      //图象直方图，共256个点
	int n, n1, n2;
	int total;                              //total为总和，累计值
	double m1, m2, sum, csum, fmax, sb;     //sb为类间方差，fmax存储最大方差值
	int k, t, q;
	int threshValue = 1;                    // 阈值
	Mat_<uchar>::iterator it = ans.begin<uchar>(), end = ans.end<uchar>();
	for (; it != end; ++it)
		++pixelNum[*it];
	//直方图平滑化
	for (k = 0; k <= 255; k++)
	{
		total = 0;
		for (t = -2; t <= 2; t++)              //与附近2个灰度做平滑化，t值应取较小的值
		{
			q = k + t;
			if (q < 0)                     //越界处理
				q = 0;
			if (q > 255)
				q = 255;
			total = total + pixelNum[q];    //total为总和，累计值
		}
		pixelNum[k] = (int)((float)total / 5.0 + 0.5);    //平滑化，左边2个+中间1个+右边2个灰度，共5个，所以总和除以5，后面加0.5是用修正值
	}
	//求阈值
	sum = csum = 0.0;
	n = 0;
	//计算总的图象的点数和质量矩，为后面的计算做准备
	for (k = 0; k <= 255; k++)
	{
		sum += (double)k * (double)pixelNum[k];     //x*f(x)质量矩，也就是每个灰度的值乘以其点数（归一化后为概率），sum为其总和
		n += pixelNum[k];                       //n为图象总的点数，归一化后就是累积概率
	}
	fmax = -1.0;                          //类间方差sb不可能为负，所以fmax初始值为-1不影响计算的进行
	n1 = 0;
	for (k = 0; k <= 255; k++)                  //对每个灰度（从0到255）计算一次分割后的类间方差sb
	{
		n1 += pixelNum[k];                //n1为在当前阈值遍前景图象的点数
		if (n1 == 0)
		{
			continue;    //没有分出前景后景
		}
		n2 = n - n1;                        //n2为背景图象的点数
		if (n2 == 0)
		{
			break;    //n2为0表示全部都是后景图象，与n1=0情况类似，之后的遍历不可能使前景点数增加，所以此时可以退出循环
		}
		csum += (double)k * pixelNum[k];    //前景的“灰度的值*其点数”的总和
		m1 = csum / n1;                     //m1为前景的平均灰度
		m2 = (sum - csum) / n2;               //m2为背景的平均灰度
		sb = (double)n1 * (double)n2 * (m1 - m2) * (m1 - m2);   //sb为类间方差
		if (sb - fmax>0.001)                  //如果算出的类间方差大于前一次算出的类间方差
		{
			fmax = sb;                    //fmax始终为最大类间方差（otsu）
			threshValue = k;              //取最大类间方差时对应的灰度的k就是最佳阈值
		}
	}
	for (it = ans.begin<uchar>(); it != end; ++it) //根据阈值二值化图像
	{
		if (*it>threshValue)
			*it = 255;
		else
			*it = 0;
	}
	return ans;
}