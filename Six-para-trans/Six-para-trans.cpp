// Six-para-trans.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
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
}KnownData[6], TargetData[6];

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

int main()
{	
	//首先进行文件的读取
	read("XYZ_origin_2.xyz", KnownData);
	read("XYZ_target_2.xyz", TargetData);
	//然后初始化B和l矩阵，计算参数
	Matrix<double, 12, 6>B;
	Matrix<double, 12, 1>l;
	for (int i = 0; i < 12; i++)
	{
		if (i % 3 == 0)
		{
			B(i, 0) = 1;
			B(i, 1) = 0;
			B(i, 2) = 0;
			B(i, 3) = 0;
			B(i, 4) = -KnownData[i / 3].Z;
			B(i, 5) = KnownData[i / 3].Y;
			l(i, 0) = TargetData[i / 3].X - KnownData[i / 3].X;
		}
		else if (i % 3 == 1)
		{
			B(i, 0) = 0;
			B(i, 1) = 1;
			B(i, 2) = 0;
			B(i, 3) = KnownData[i / 3].Z;
			B(i, 4) = 0;
			B(i, 5) = -KnownData[i / 3].X;
			l(i, 0) = TargetData[i / 3].Y - KnownData[i / 3].Y;
		}
		else
		{
			B(i, 0) = 0;
			B(i, 1) = 0;
			B(i, 2) = 1;
			B(i, 3) = -KnownData[i / 3].Y;
			B(i, 4) = KnownData[i / 3].X;
			B(i, 5) = 0;
			l(i, 0) = TargetData[i / 3].Z - KnownData[i / 3].Z;
		}
	}
	MatrixXd BTB = B.transpose()*B;
	MatrixXd W = B.transpose()*l;
	MatrixXd Para = BTB.inverse()*W;
	MatrixXd V = B * Para - l;
	MatrixXd sigma2 = V.transpose()*V;

	//定义坐标转换用的新矩阵
	Matrix<double, 6, 6>MatB;
	Matrix<double, 6, 1>l0;//origin矩阵
	Matrix<double, 6, 1>l2;//target矩阵
	for (int i = 0; i < 6; i++)
	{
		if (i % 3 == 0)
		{
			MatB(i, 0) = 1;
			MatB(i, 1) = 0;
			MatB(i, 2) = 0;
			MatB(i, 3) = 0;
			MatB(i, 4) = -KnownData[i / 3 + 4].Z;
			MatB(i, 5) = KnownData[i / 3 + 4].Y;
			l0(i, 0) = KnownData[(i / 3) + 4].X;
			l2(i, 0) = TargetData[(i / 3) + 4].X;
		}
		else if (i % 3 == 1)
		{
			MatB(i, 0) = 0;
			MatB(i, 1) = 1;
			MatB(i, 2) = 0;
			MatB(i, 3) = KnownData[i / 3 + 4].Z;
			MatB(i, 4) = 0;
			MatB(i, 5) = -KnownData[i / 3 + 4].X;
			l0(i, 0) = KnownData[(i / 3) + 4].Y;
			l2(i, 0) = TargetData[(i / 3) + 4].Y;
		}
		else
		{
			MatB(i, 0) = 0;
			MatB(i, 1) = 0;
			MatB(i, 2) = 1;
			MatB(i, 3) = -KnownData[i / 3 + 4].Y;
			MatB(i, 4) = KnownData[i / 3 + 4].X;
			MatB(i, 5) = 0;
			l0(i, 0) = KnownData[(i / 3) + 4].Z;
			l2(i, 0) = TargetData[(i / 3) + 4].Z;
		}
	}
	MatrixXd x1 = MatB * Para + l0;
	MatrixXd V2 = x1 - l2;
	Matrix<double, 3, 6>B3;
	Matrix<double, 3, 1>l3;
	B3 << 1, 0, 0, 0, -2894337.6030, 5496337.0138,
		0, 1, 0, 2894337.6030, 0, 2100337.2849,
		0, 0, 1, -5496337.0138, -2100337.2849, 0;
	l3 << -2100337.2849, 5496337.0137, 2894337.6030;
	MatrixXd x3 = B3 * Para + l3;

	//输出部分
	cout.precision(10);
	cout << "下面输出条件数和冗余数：";
	cout << "n=12    r=n-t=6" << endl;
	cout << "下面输出参数：" << endl;
	cout << "DX=" << Para(0, 0) << endl;
	cout << "DY=" << Para(1, 0) << endl;
	cout << "DR=" << Para(2, 0) << endl;
	cout << "RX=" << Para(3, 0) << endl;
	cout << "RY=" << Para(4, 0) << endl;
	cout << "RZ=" << Para(5, 0) << endl;
	cout << "验后单位权中误差=" << sqrt(sigma2(0, 0) / 6) << endl << endl;
	cout << "下面输出每个公共点坐标残差：" << endl;
	for (int i = 0; i < 12; i += 3)
	{
		cout << KnownData[i / 3].name << "," << "Vx=" << V(i, 0) << " Vy=" << V(i + 1, 0) << " Vz=" << V(i + 2, 0) << endl;
	}
	cout << "下面输出检核点坐标残差：" << endl;
	cout << "Vx= " << V2(0, 0) << " Vy= " << V2(1, 0) << " Vz= " << V2(2, 0) << endl;
	cout << "Vx= " << V2(3, 0) << " Vy= " << V2(4, 0) << " Vz= " << V2(5, 0) << endl;
	cout.precision(11);
	cout << endl << "学号相关的坐标  x:" << x3(0, 0) << "  y:" << x3(1, 0) << "  z:" << x3(2, 0) << endl;
	system("pause");
	return 0;
}