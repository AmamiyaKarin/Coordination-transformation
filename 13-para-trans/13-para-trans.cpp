// 13-para-trans.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <string>
#include <cmath>
using namespace Eigen;
using namespace std;

//定义一个点的结构体数组
struct Data {
	string name;
	double X;
	double Y;
	double Z;
}KnownData[17], TargetData[17];

//下面是读取文件的函数
void read(string FilePath, struct Data a[17]) {
	ifstream infile(FilePath, ios::in);
	int order = 0;//order变量用于向后遍历字符串
	int num = 0;//num变量用于检测这是这一行的第几个逗号，因此可以区分这一段数据应该赋给谁
	string s;
	cout << "reading from the file" << endl;
	while (getline(infile, s)) //逐行读取文件信息
	{
		num = 0;
		//定义两个字符串，用于之后的读取
		string string1 = s;
		string string2 = "";
		for (int i = 0; i < string1.length(); i++)//提取出一行后，一个个往后遍历，由于以逗号分隔，因此需要做if判断
		{
			if (string1[i] == ',')//如果检测到了逗号，停止读取，处理之前的一段数据
			{
				if (num == 0)//这是第一个逗号，因此前面一段应该是点名
				{
					a[order].name = string2;
				}
				if (num == 1)//这是第二个逗号，因此前面一段应该是x坐标
				{
					a[order].X = atof(string2.c_str());
				}
				if (num == 2)//这是第三个逗号，因此前面一段应该是y坐标
				{
					a[order].Y = atof(string2.c_str());
				}
				string2 = "";
				num++;
			}
			if (string1[i] != ',')//如果该字符不是逗号，则将其添加到string2中去，检测到逗号时停止读取，将之前读取的一段转化成相应的数据类型然后赋值
			{
				string2 += string1[i];
			}
		}
		if (string2 != "")
			a[order].Z = atof(string2.c_str());
		order++;
	}
	cout << "succeed!" << endl;
	infile.close();
	order = 0;//order清零 便于下次文件读取
}

double sum(Matrix<double,13,1>a)
{
	double sum=0;
	for (int i = 0; i < 13; i++)
	{
		sum += abs(a(i,0));
	}
	return sum;
}
int main()
{
	//首先进行文件的读取
	read("XYZ_origin_3.xyz", KnownData);
	read("XYZ_target_3.xyz", TargetData);
	//建立矩阵，要将限制条件方程和误差方程组合，形成新的系数矩阵
	Matrix<double, 39, 13>B;
	Matrix<double, 39, 1>l;
	Matrix<double, 39, 39>P;
	Matrix<double, 6, 13>C;
	Matrix<double, 13, 1>BTPl;
	MatrixXd NBB;
	MatrixXd xhat;
	int m = 0;
	double DX = 0, DY = 0, DZ = 0, a11 = 1, a12 = 0, a13 = 0, a21 = 0, a22 = 1, a23 = 0, a31 = 0, a32 = 0, a33 = 1, k = 1;
	//设定矩阵初值，便于之后对特定位置赋值
	for (int i = 0; i < 39; i++)
	{
		for (int j = 0; j < 13; j++)
		{
			B(i, j) = 0;
		}
		l(i, 0) = 0;
	}
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 13; j++)
		{
			C(i, j) = 0;
		}
		BTPl(i, 0) = 0;
	}
	for (int i = 0; i < 39; i++)
	{
		if(i<33)P(i, i) = 1;
		else P(i, i) = 10000000;
	}
	//接下来进行迭代
	do {
		//先对B矩阵进行初始化赋值
		for (int i = 0; i < 33; i++)
		{
			if (i % 3 == 0)
			{
				B(i, 0) = 1;
				B(i, 1) = 0;
				B(i, 2) = 0;
				B(i, 3) = a11 * KnownData[i / 3].X + a12 * KnownData[i / 3].Y + a13 * KnownData[i / 3].Z;
				B(i, 4) = k * KnownData[i / 3].X;
				B(i, 5) = k * KnownData[i / 3].Y;
				B(i, 6) = k * KnownData[i / 3].Z;
			}
			else if (i % 3 == 1)
			{
				B(i, 0) = 0;
				B(i, 1) = 1;
				B(i, 2) = 0;
				B(i, 3) = a21 * KnownData[i / 3].X + a22 * KnownData[i / 3].Y + a23 * KnownData[i / 3].Z;
				B(i, 7) = k * KnownData[i / 3].X;
				B(i, 8) = k * KnownData[i / 3].Y;
				B(i, 9) = k * KnownData[i / 3].Z;
			}
			else
			{
				B(i, 0) = 0;
				B(i, 1) = 0;
				B(i, 2) = 1;
				B(i, 3) = a31 * KnownData[i / 3].X + a32 * KnownData[i / 3].Y + a33 * KnownData[i / 3].Z;
				B(i, 10) = k * KnownData[i / 3].X;
				B(i, 11) = k * KnownData[i / 3].Y;
				B(i, 12) = k * KnownData[i / 3].Z;
			}
		}
		//由于是将条件方程也计入误差方程，因此将B矩阵增广收纳进C矩阵，这里要对C矩阵的部分单独赋值
		B(33, 4) = 2 * a11; B(33, 5) = 2 * a12; B(33, 6) = 2 * a13;
		B(34, 7) = 2 * a21; B(34, 8) = 2 * a22; B(34, 9) = 2 * a23;
		B(35, 10) = 2 * a31; B(35, 11) = 2 * a32; B(35, 12) = 2 * a33;
		B(36, 4) = a12; B(36, 5) = a11; B(36, 7) = a22; B(36, 8) = a21; B(36, 10) = a32; B(36, 11) = a31;
		B(37, 4) = a13; B(37, 6) = a11; B(37, 7) = a23; B(37, 9) = a21; B(37, 10) = a33; B(37, 12) = a31;
		B(38, 5) = a13; B(38, 6) = a12; B(38, 8) = a23; B(38, 9) = a22; B(38, 11) = a33; B(38, 12) = a32;

		//然后对l矩阵初始化
		for (int i = 0; i < 33; i++)
		{
			if (i % 3 == 0) l(i, 0) = TargetData[i / 3].X - k * (a11*KnownData[i / 3].X + a12 * KnownData[i / 3].Y + a13 * KnownData[i / 3].Z) - DX;
			else if (i % 3 == 1) l(i, 0) = TargetData[i / 3].Y - k * (a21*KnownData[i / 3].X + a22 * KnownData[i / 3].Y + a23 * KnownData[i / 3].Z) - DY;
			else l(i, 0) = TargetData[i / 3].Z - k * (a31*KnownData[i / 3].X + a32 * KnownData[i / 3].Y + a33 * KnownData[i / 3].Z) - DZ;

		}
		//同理，l矩阵与W矩阵也进行了合并，这里赋值W矩阵的部分
		l(33, 0) = -(a11 * a11 + a12 * a12 + a13 * a13 - 1);
		l(34, 0) = -(a21 * a21 + a22 * a22 + a23 * a23 - 1);
		l(35, 0) = -(a31 * a31 + a32 * a32 + a33 * a33 - 1);
		l(36, 0) = -(a11 * a12 + a21 * a22 + a31 * a32);
		l(37, 0) = -(a11 * a13 + a21 * a23 + a31 * a33);
		l(38, 0) = -(a12 * a13 + a22 * a23 + a32 * a33);


		//进行矩阵的计算
		NBB = B.transpose()*P*B;
		BTPl = B.transpose()*P*l;
		xhat = NBB.inverse()*BTPl;

		//参数的重新赋值，作为下一次迭代的新值
		DX += xhat(0, 0);
		DY += xhat(1, 0);
		DZ += xhat(2, 0);
		k += xhat(3, 0);
		a11 += xhat(4, 0);
		a12 += xhat(5, 0);
		a13 += xhat(6, 0);
		a21 += xhat(7, 0);
		a22 += xhat(8, 0);
		a23 += xhat(9, 0);
		a31 += xhat(10, 0);
		a32 += xhat(11, 0);
		a33 += xhat(12, 0);
		m++;//m用以记录迭代次数
	} while (sum(xhat) > 0.00000001);//当精度符合要求时，跳出循环
	if (m > 2000)cout << "failed!" << endl;
	cout << "一共迭代了" << m << "次" << endl;

	MatrixXd V = B * xhat - l;
	MatrixXd sigma2 = V.transpose()*V;
	//构建新矩阵对检核点做运算，计算误差
	Matrix<double, 3, 3>B2;
	Matrix<double, 18, 1>V2;
	Matrix<double, 3, 1>l2;
	l2 << DX, DY, DZ;
	B2 << a11, a12, a13, a21, a22, a23, a31, a32, a33;
	for (int i = 0; i < 18; i += 3)
	{
		Matrix<double, 3, 1>x1;
		x1 << KnownData[i / 3 + 11].X, KnownData[i / 3 + 11].Y, KnownData[i / 3 + 11].Z;
		MatrixXd x2 = k * B2*x1 + l2;
		V2(i, 0) = TargetData[i / 3 + 11].X - x2(0, 0);
		V2(i + 1, 0) = TargetData[i / 3 + 11].Y - x2(1, 0);
		V2(i + 2, 0) = TargetData[i / 3 + 11].Z - x2(2, 0);
	}





	//最后进行输出
	cout.precision(11);
	cout << "下面输出条件数和冗余数：" << endl;
	cout << "n=33       r=n-t=26" << endl;
	cout << "下面输出参数: "<<endl;
	cout << "DX=" << DX << endl;
	cout << "DY=" << DY << endl;
	cout << "DZ=" << DZ << endl;
	cout << "k=" << k << endl;
	cout << "a11=" << a11 << endl;
	cout << "a12=" << a12 << endl;
	cout << "a13=" << a13 << endl;
	cout << "a21=" << a21 << endl;
	cout << "a22=" << a22 << endl;
	cout << "a23=" << a23 << endl;
	cout << "a31=" << a31 << endl;
	cout << "a32=" << a32 << endl;
	cout << "a33=" << a33 << endl;
	cout << "验后单位权中误差为" << sqrt(sigma2(0, 0)/26.00) << endl;
	cout << "下面输出公共点坐标残差:" << endl;
	for (int i = 0; i < 33; i += 3)
	{
		cout << KnownData[i / 3].name << "," << "Vx=" << V(i, 0) << " Vy=" << V(i + 1, 0) << " Vz=" << V(i + 2, 0) << endl;
	}
	cout << "下面输出检核点坐标残差:" << endl;
	for (int i = 0; i < 18; i += 3)
	{
		cout << KnownData[i / 3 + 11].name << "," << " Vx=" << V2(i, 0) << " Vy=" << V2(i + 1, 0) << " Vz=" << V2(i + 2, 0) << endl;
	}
	system("pause");
	return 0;
}

