#include<iostream>
#include<math.h>
#include<stdio.h>
#include<fstream>

using namespace std;

int node_num,pq_num,pv_num;                                                                                        //����ڵ�����pq�ڵ�����pv�ڵ���
double Ac;                                                                                                                        //���徫��
double G[10][10],B[10][10];                                                                                               //����絼�͵��ɾ���
double c1,c2,d1,d2,alpha,M;                                                                                              //���帴���������
double Jacob[10][10],D[10],deltaU[10];                                                                              //�Ÿ��Ⱦ��󣬲�ƽ��������� 
struct Node
{
	int num,t;                                                                                                                       //numΪ�ڵ�ţ�tΪ�ڵ�����
	double P,Q,e,f,U,theta,delta_P,delta_Q,delta_e,delta_f,delta_U;                                        //�ڵ��й����޹����ʣ���ѹ�ݡ����������ѹģֵ����ѹ����Լ���ƽ�����
}Node[100];

void Read_data(double pos[])
{
	int i;
	ifstream in("Input.txt");
	for(i=0;i<100;i++)
	{
		in>>pos[i];
	}
	node_num=pos[0];
	pq_num=pos[2];
	pv_num=pos[3];
	Ac=pos[5];
	for(i=1;i<node_num;i++)
	{
		Node[i].num=pos[i*5+1];
		Node[i].t=pos[i*5+2];
		Node[i].P=pos[i*5+3];
		Node[i].Q=pos[i*5+4];
		Node[i].e=pos[i*5+5];
		Node[i].f=pos[i*5+6];
		Node[node_num].num=pos[node_num*5+1];
		Node[node_num].t=pos[node_num*5+2];
		Node[node_num].e=pos[node_num*5+3];
		Node[node_num].f=pos[node_num*5+4];
	}
	in.close();
}

void Form_Y(double ad[],int node)
{
	int i,j;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			G[i][j]=0;
			B[i][j]=0;
		}
	}
	//��ʼ���ڵ㵼�ɾ���
	G[0][1]=-ad[25]/(ad[25]*ad[25]+ad[26]*ad[26]);       //�����ɼ���
	G[0][2]=-ad[30]/(ad[30]*ad[30]+ad[31]*ad[31]);
	G[1][0]=G[0][1];
	G[2][0]=G[0][2];
	G[1][1]=-G[0][1];                                                      //�Ե��ɼ���
	G[2][2]=-G[0][2];
	G[0][0]=-G[0][1]-G[0][2];
	B[0][1]=ad[26]/(ad[25]*ad[25]+ad[26]*ad[26]);
	B[0][2]=ad[31]/(ad[30]*ad[30]+ad[31]*ad[31]);
	B[1][0]=B[0][1];
	B[2][0]=B[0][2];
	B[1][1]=-B[0][1];
	B[2][2]=-B[0][2];
	B[0][0]=-B[0][1]-B[0][2];
	//�γɽڵ㵼�ɾ���
	cout<<"�絼����"<<endl;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			cout<<G[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<"���ɾ���"<<endl;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			cout<<B[i][j]<<" ";
		}
		cout<<endl;
	}
}

 double mozhi(double a0,double b0)                      /*������ģֵ����*/
 { 	 
	 M=sqrt(a0*a0+b0*b0); 
	 return 0;
 } 
 double ji(double a1,double b1,double a2,double b2)     /*�����������*/
 {   c1=a1*a2-b1*b2;
	 d1=a1*b2+a2*b1;
	 return 0;
 }
 double shang(double a3,double b3,double a4,double b4)  /*�������̺���*/
 {   c2=(a3*a4+b3*b4)/(a4*a4+b4*b4);
     d2=(a4*b3-a3*b4)/(a4*a4+b4*b4);
	 return 0;
 }

double Calculate_Unbalanced_Para()      //����PQ�ڵ㲻ƽ�����delta_P��delta_Q��PV�ڵ��delta_P��deltaU^2
{
	int i,j,k,a,g;
	g=pq_num+pv_num;
	double Cc[10],Dd[10];
	for(i=1;i<=node_num;i++)                       /*����PQ�ڵ㲻ƽ����*/
	{
		if(Node[i].t==1)
		{
			k=Node[i].num;
		    Cc[k]=Dd[k]=0;
		    for(j=1;j<=node_num;j++)
			{
				for(a=1;a<=node_num;a++)		  
				{
					if(Node[a].num==j)
					    break; 
				}                                         
			    ji(G[k][j],-B[k][j],Node[a].e,-Node[a].f);
			    Cc[k]+=c1;
				Dd[k]+=d1;                  
			}	                                    
				ji(Node[i].e,Node[i].f,Cc[k],Dd[k]);
				Node[i].delta_P=Node[i].P-c1;
				Node[i].delta_Q=Node[i].Q-d1;      
		}
		if(Node[i].t==2)                      /*����PV�ڵ㲻ƽ����*/
		{
			k=Node[i].num;
		    Cc[k]=Dd[k]=0;
		    for(j=1;j<=node_num;j++)
			{
				for(a=1;a<=node_num;a++)		  
				{
					if(Node[a].num==j)
					    break; 
				}                                         
			    ji(G[k][j],-B[k][j],Node[a].e,-Node[a].f);
			    Cc[k]+=c1;
				Dd[k]+=d1;                  
			}	                                    
				ji(Node[i].e,Node[i].f,Cc[k],Dd[k]);
				Node[i].delta_P=Node[i].P-c1;
				Node[i].Q=d1;
				Node[i].delta_U=Node[i].U*Node[i].U-(Node[i].e*Node[i].e+Node[i].f*Node[i].f);    
		}
    for(i=1;i<=pq_num;i++)                                   /*�γɲ�ƽ��������D[M]*/
    {   		
		D[2*i-1]=Node[i].delta_P;
	    D[2*i]=Node[i].delta_Q; 
	}
    for(i=pq_num+1;i<=g;i++)
    {  
		D[2*i-1]=Node[i].delta_P;
	    D[2*i]=Node[i].delta_U; 
	}
	for(i=1;i<=pq_num;i++)
	{	
		k=Node[i].num;
		cout<<"delta_P["<<i<<"]="<<D[2*k-1]<<endl;
		cout<<"delta_Q["<<i<<"]="<<D[2*k]<<endl;
	} 
	for(i=pq_num+1;i<=g;i++)
	{	
		k=Node[i].num;
		cout<<"delta_P["<<i<<"]="<<D[2*k-1]<<endl;
		cout<<"delta_U["<<i<<"]="<<D[2*k]<<endl;
	} 
}
	return 0;
}

void main()
{
	double ad[100];
	int i,j;
	for(i=0;i<100;i++)
	{
		ad[i]=0;
	}
	//��ʼ����������
	Read_data(ad);
	//��·�ڵ���
	Form_Y(ad,node_num);
	//����ţ��-����ѷ���еĲ�ƽ�����
	Calculate_Unbalanced_Para();
	system("pause");
}