#include<iostream>
#include<math.h>
#include<stdio.h>
#include<fstream>

using namespace std;

int  Read_data(double pos[])
{
	int i;
	ifstream in("C:\\Users\\Administrator\\Documents\\Visual Studio 2010\\Projects\\test\\Input.txt");
	for(i=0;i<100;i++)
	{
		in>>pos[i];
	}
	for(i=0;i<100;i++)
	{
		if(pos[i]==100)
			break;
	}
	in.close();
	return(i);
}

void main()
{
	double ad[100];
	int n,i,j,node;
	for(i=0;i<100;i++)
	{
		ad[i]=100;
	}
	n=Read_data(ad);
	double G[10],B[10];
	for(i=0;i<n;i++)
	{
		G[0]=ad[25];
		G[1]=ad[30];
		B[0]=ad[26];
		B[1]=ad[31];
		node=ad[0];
	}
	double gb[10][10];
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			gb[i][j]=0;
		}
	}
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			gb[0][0]=G[0]+G[1];
			gb[1][1]=G[0];
			gb[2][2]=G[1];
		}
	}
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			cout<<gb[i][j]<<" ";
		}
		cout<<endl;
	}
	system("pause");
}