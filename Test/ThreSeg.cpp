#include "threSeg.h"

Mat threSeg(const Mat &Image)
{
	Mat ans;
	cvtColor(Image, ans, CV_RGB2GRAY);
	int *pixelNum = (int*)calloc(256, sizeof(int));                      //ͼ��ֱ��ͼ����256����
	int n, n1, n2;
	int total;                              //totalΪ�ܺͣ��ۼ�ֵ
	double m1, m2, sum, csum, fmax, sb;     //sbΪ��䷽�fmax�洢��󷽲�ֵ
	int k, t, q;
	int threshValue = 1;                    // ��ֵ
	Mat_<uchar>::iterator it = ans.begin<uchar>(), end = ans.end<uchar>();
	for (; it != end; ++it)
		++pixelNum[*it];
	//ֱ��ͼƽ����
	for (k = 0; k <= 255; k++)
	{
		total = 0;
		for (t = -2; t <= 2; t++)              //�븽��2���Ҷ���ƽ������tֵӦȡ��С��ֵ
		{
			q = k + t;
			if (q < 0)                     //Խ�紦��
				q = 0;
			if (q > 255)
				q = 255;
			total = total + pixelNum[q];    //totalΪ�ܺͣ��ۼ�ֵ
		}
		pixelNum[k] = (int)((float)total / 5.0 + 0.5);    //ƽ���������2��+�м�1��+�ұ�2���Ҷȣ���5���������ܺͳ���5�������0.5��������ֵ
	}
	//����ֵ
	sum = csum = 0.0;
	n = 0;
	//�����ܵ�ͼ��ĵ����������أ�Ϊ����ļ�����׼��
	for (k = 0; k <= 255; k++)
	{
		sum += (double)k * (double)pixelNum[k];     //x*f(x)�����أ�Ҳ����ÿ���Ҷȵ�ֵ�������������һ����Ϊ���ʣ���sumΪ���ܺ�
		n += pixelNum[k];                       //nΪͼ���ܵĵ�������һ��������ۻ�����
	}
	fmax = -1.0;                          //��䷽��sb������Ϊ��������fmax��ʼֵΪ-1��Ӱ�����Ľ���
	n1 = 0;
	for (k = 0; k <= 255; k++)                  //��ÿ���Ҷȣ���0��255������һ�ηָ�����䷽��sb
	{
		n1 += pixelNum[k];                //n1Ϊ�ڵ�ǰ��ֵ��ǰ��ͼ��ĵ���
		if (n1 == 0)
		{
			continue;    //û�зֳ�ǰ����
		}
		n2 = n - n1;                        //n2Ϊ����ͼ��ĵ���
		if (n2 == 0)
		{
			break;    //n2Ϊ0��ʾȫ�����Ǻ�ͼ����n1=0������ƣ�֮��ı���������ʹǰ���������ӣ����Դ�ʱ�����˳�ѭ��
		}
		csum += (double)k * pixelNum[k];    //ǰ���ġ��Ҷȵ�ֵ*����������ܺ�
		m1 = csum / n1;                     //m1Ϊǰ����ƽ���Ҷ�
		m2 = (sum - csum) / n2;               //m2Ϊ������ƽ���Ҷ�
		sb = (double)n1 * (double)n2 * (m1 - m2) * (m1 - m2);   //sbΪ��䷽��
		if (sb - fmax>0.001)                  //����������䷽�����ǰһ���������䷽��
		{
			fmax = sb;                    //fmaxʼ��Ϊ�����䷽�otsu��
			threshValue = k;              //ȡ�����䷽��ʱ��Ӧ�ĻҶȵ�k���������ֵ
		}
	}
	for (it = ans.begin<uchar>(); it != end; ++it) //������ֵ��ֵ��ͼ��
	{
		if (*it>threshValue)
			*it = 255;
		else
			*it = 0;
	}
	return ans;
}