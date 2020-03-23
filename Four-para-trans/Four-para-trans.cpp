// Four-para-trans.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
#include "pch.h"
#include <iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<Eigen/Dense>
#include<cmath>
using namespace Eigen;
using namespace std;

//定义一个点的结构体数组
struct Data {
	string name;
	double X;
	double Y;
	double Z;
}KnownData[8],TargetData[8];

//下面是读取文件的函数
void read( string FilePath, struct Data a[8]) {
	ifstream infile(FilePath, ios::in);
	int order = 0;//order变量用于向后遍历字符串
	int num= 0;//num变量用于检测这是这一行的第几个逗号，因此可以区分这一段数据应该赋给谁
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

int main()
{
	//分别读取两个文件作为初始数据
	read("XYZ_origin_1.xyz", KnownData);
	read("XYZ_target_1.xyz", TargetData);

	//定义矩阵，
	Matrix<double, 8 , 4>B;
	Matrix < double , 8, 1 > l;

	//矩阵初始化
	for (int i = 0; i < 8; i++)
	{
		if (i % 2 == 0)
		{
			B(i, 0) = 1;
			B(i, 1) = 0;
			B(i, 2) = KnownData[i / 2].Y;
			B(i, 3) = KnownData[i / 2].X;
			l(i, 0) = TargetData[i / 2].X - KnownData[i / 2].X;
		}
		else
		{
			B(i, 0) = 0;
			B(i, 1) = 1;
			B(i, 2) = -KnownData[i / 2].X;
			B(i, 3) = KnownData[i / 2].Y;
			l(i, 0) = TargetData[i / 2].Y - KnownData[i / 2].Y;
		}
	}

	//下面进行矩阵计算，并进行内符合指标的计算，此处P为单位矩阵，故略去
	MatrixXd BTB = B.transpose()*B;
	MatrixXd W = B.transpose()*l;
	MatrixXd Para = BTB.inverse()*W;
	MatrixXd V = B * Para - l;
	MatrixXd sigma2 = V.transpose()*V;

	//定义坐标转换用的新矩阵
	Matrix < double , 8, 4 > MatB;
	Matrix<double, 8, 1>l0;//origindata
	Matrix<double, 8, 1>l2;//targetdata
	//下面进行系数矩阵MatB，l0和l2矩阵的初始化
	for (int i = 0; i < 8; i++)
	{
		if (i % 2 == 0)
		{
			MatB(i, 0) = 1;
			MatB(i, 1) = 0;
			MatB(i, 2) = KnownData[i / 2 + 4].Y;
			MatB(i, 3) = KnownData[i / 2 + 4].X;
			l0(i, 0) = KnownData[(i / 2) + 4].X;
			l2(i, 0) = TargetData[(i / 2) + 4].X;
		}
		else
		{
			MatB(i, 0) = 0;
			MatB(i, 1) = 1;
			MatB(i, 2) = -KnownData[i / 2 + 4].X;
			MatB(i, 3) = KnownData[i / 2 + 4].Y;
			l0(i, 0) = KnownData[(i / 2) + 4].Y;
			l2(i, 0) = TargetData[(i / 2) + 4].Y;
		}
	}
	MatrixXd x1 = MatB * Para + l0;
	MatrixXd v2 = x1 - l2;
	Matrix < double , 2, 4 > B3;
	Matrix<double, 2, 1>l3;
	B3 << 1, 0.000, 
		118337.826,115337.237, 
		0.000, 1, 
		-115337.237, 118337.826;
	l3 << 115337.237, 118337.826;
	MatrixXd x3 = B3 * Para + l3;

	//输出部分
	cout.precision(11);
	cout << "下面输出条件数与冗余度：" << endl;
	cout << "n=8   r=n-t=4" << endl;
	cout << "下面输出计算得出的参数：" << endl;
	cout << "DX=" << Para(0, 0) << endl;
	cout << "DY=" << Para(1, 0) << endl;
	cout << "DR=" << Para(2, 0) << endl;
	cout << "DK=" << Para(3, 0) << endl;
	cout << "验后单位权中误差=" << sqrt(sigma2(0, 0) / 4) << endl << endl;
	cout << "下面输出公共点的坐标残差：" << endl;
	for (int i = 0; i < 8; i += 2)
	{
		cout << KnownData[i / 2].name << "," << "Vx=" << ":" << V(i, 0) << " Vy=" << ":" << V(i + 1) << endl;
	}
	cout << "下面输出检核点的坐标残差：" << endl;
	for (int i = 0; i < 8; i += 2)
	{
		cout << KnownData[i / 2 + 4].name << "," << "Vx= " << v2(i, 0) << "," << "Vy= " << v2(i + 1, 0) << endl;
	}
	cout<< "下面输出由学号生成的坐标转换结果：" << endl;
	cout<<"("<<x3(0, 0) << "," << x3(1, 0) <<")"<< endl;
	system("pause");
	return 0;
}