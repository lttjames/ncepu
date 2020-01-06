#include<iostream>
#include<math.h>
#include<stdio.h>
#include<fstream>

using namespace std;

struct Node
{
	int node_num;
	double U,P,Q,theta;
};

int Read_data(double pos[])
{
	int i;
	int node_num;
	ifstream in("Input.txt");
	for(i=0;i<100;i++)
	{
		in>>pos[i];
	}
	node_num=pos[0];
	in.close();
	return(node_num);
}

void Form_Y(double ad[],int node)
{
	double G[10][10],B[10][10];
	int i,j;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			G[i][j]=0;
			B[i][j]=0;
		}
	}
	//初始化节点导纳矩阵
	G[0][1]=-ad[25]/(ad[25]*ad[25]+ad[26]*ad[26]);       //互导纳计算
	G[0][2]=-ad[30]/(ad[30]*ad[30]+ad[31]*ad[31]);
	G[1][0]=G[0][1];
	G[2][0]=G[0][2];
	G[1][1]=-G[0][1];                                                      //自导纳计算
	G[2][2]=-G[0][2];
	G[0][0]=-G[0][1]-G[0][2];
	B[0][1]=ad[26]/(ad[25]*ad[25]+ad[26]*ad[26]);
	B[0][2]=ad[31]/(ad[30]*ad[30]+ad[31]*ad[31]);
	B[1][0]=B[0][1];
	B[2][0]=B[0][2];
	B[1][1]=-B[0][1];
	B[2][2]=-B[0][2];
	B[0][0]=-B[0][1]-B[0][2];
	//形成节点导纳矩阵
	cout<<"电导矩阵"<<endl;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			cout<<G[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<"电纳矩阵"<<endl;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			cout<<B[i][j]<<" ";
		}
		cout<<endl;
	}
}

double Calculate_Unbalanced_Para()      //计算PQ节点不平衡分量deltaPi、deltaQi和PV节点的deltaPi、deltaUi^2
{
	return 0;
}

void main()
{
	double ad[100];
	int i,j,node;
	for(i=0;i<100;i++)
	{
		ad[i]=0;
	}
	//初始化输入数组
	node=Read_data(ad);
	//线路节点数
	Form_Y(ad,node);
	system("pause");
}