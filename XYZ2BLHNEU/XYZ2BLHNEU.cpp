// XYZ2BLHNEU.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include "pch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <string>
#include <cmath>
const double pi = acos(-1);
using namespace Eigen;
using namespace std;


struct Data {
	string name;
	double X;
	double Y;
	double Z;
}XYZ[14435],NEU[14435],BLH[14435];

//下面是读取文件的函数
void read(string FilePath, struct Data a[6]) {
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
//下面是输出文件的函数
void out(string FilePath, struct Data a[14435])
{
	ofstream outfile;
	outfile.open(FilePath);
	cout << "outputting to the file" << endl;
	for (int i = 0; i < 14435; i++)
	{
		outfile << a[i].name << "," << a[i].X << "," << a[i].Y << "," << a[i].Z  << endl;
	}
	cout << "succeed" << endl;
	outfile.close();
}

double Rad2Deg(double rad)
{
	double deg = rad * 180 / pi;
	return deg;
}
double Deg2Rad(double deg)
{
	double rad = deg * pi / 180;
	return rad;
}
double ComputeE2()
{
	double f = 1 / 298.257222101;
	double e2 = 2 * f - f * f;
	return e2;
}
double ComputeN(double B)
{
	double a = 6378137;
	double e2 = ComputeE2();
	double sinB = sin(Deg2Rad(B));
	double W = sqrt(1 - e2 * sinB * sinB);
	double N = a / W;
	return N;
}
double ComputeH(double X, double Y, double B)
{
	double N = ComputeN(B);
	double H = (sqrt(X * X + Y * Y) / cos(Deg2Rad(B))) - N;
	return H;
}
double ComputeB(double B0, double X, double Y, double Z, double H)
{
	double e2 = ComputeE2();
	double N = ComputeN(B0);
	double upper = Z * (N + H);
	double under = sqrt(X * X + Y * Y) * (N * (1 - e2) + H);
	double B = atan(upper / under);
	B = Rad2Deg(B);
	return B;
}
double ComputeL(double X, double Y)
{
	double L = 0;
	if (abs(X) < 1e-15 && Y > 0)
		L = 0;
	else if (abs(X) < 1e-15 && Y < 0)
		L = 180;
	else
	{
		L = atan(Y / X);
		L = Rad2Deg(L);
		if (Y < 0 && X > 0)
			L = L + 360;
		if (Y > 0 && X < 0)
			L += 180;
		if (X < 0 && Y < 0)
			L += 180;
	}
	return L;
}

int main()
{
	read("XYZ2BLHNEU.xyz", XYZ);
	for (int i = 0; i < 14435; i++)
	{
		double B0 = XYZ[i].Z/(sqrt(XYZ[i].X * XYZ[i].X + XYZ[i].Y * XYZ[i].Y));
		double deltaB = 90;
		BLH[i].name = XYZ[i].name;
		BLH[i].Y = ComputeL(XYZ[i].X, XYZ[i].Y);
		do
		{
			BLH[i].Z = ComputeH(XYZ[i].X, XYZ[i].Y, B0);
			BLH[i].X = ComputeB(B0, XYZ[i].X, XYZ[i].Y, XYZ[i].Z, BLH[i].Z);
			deltaB = BLH[i].X - B0;
			B0 = BLH[i].X;
		} while (abs(deltaB) > 1e-10);
	}
	Matrix<double, 3, 3>S;
	Matrix<double, 3, 1>delta;
	Matrix<double, 3, 1>X;
	for (int i = 0; i < 14435; i++)
	{
		
		S << -sin(Deg2Rad(BLH[0].Y)), cos(Deg2Rad(BLH[0].Y)), 0,
			-sin(Deg2Rad(BLH[0].X))*cos(Deg2Rad(BLH[0].Y)), -sin(Deg2Rad(BLH[0].X))*sin(Deg2Rad(BLH[0].Y)), cos(Deg2Rad(BLH[0].X)),
			cos(Deg2Rad(BLH[0].X))*cos(Deg2Rad(BLH[0].Y)), cos(Deg2Rad(BLH[0].X))*sin(Deg2Rad(BLH[0].Y)), sin(Deg2Rad(BLH[0].X));
		delta << (XYZ[i].X - XYZ[0].X), (XYZ[i].Y - XYZ[0].Y), (XYZ[i].Z - XYZ[0].Z);
		NEU[i].name = XYZ[i].name;
		X = S * delta;
		NEU[i].X = X(1, 0);
		NEU[i].Y = X(0, 0);
		NEU[i].Z = X(2, 0);
	}
	out("BLH.xyz", BLH);
	out("NEU.xyz", NEU);
	system("pause");
	return 0;
}

