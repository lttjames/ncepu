#include<iostream>
#include<math.h>
#include<stdio.h>
#include<fstream>

using namespace std;

ofstream out("Output.txt");

int node_num,branch_num,pq_num,pv_num;                                                                     //定义节点数，支路数，pq_num节点数，pv节点数
double Ac,max1;                                                                                                                //定义精度
double G[10][10],B[10][10];                                                                                                //定义电导和电纳矩阵
double c1,c2,d1,d2,alpha,M;                                                                                              //定义复数计算变量
double Jacob[100][100],D[10],deltaU[10];                                                                          //雅各比矩阵，不平衡分量矩阵 
int g;                                                                                                                                  //定义除去平衡节点的节点总数
struct Node
{
	int num,t;                                                                                                                       //num为节点号，t为节点类型
	double P,Q,e,f,U,theta,delta_P,delta_Q,delta_e,delta_f,delta_U;                                        //节点有功、无功功率，电压纵、横分量，电压模值，电压相角以及不平衡分量
}Node[10];
struct branch
{
	int numb;                                                                                                                       //numb为支路号
	int z1,z2;                                                                                                                        //z1、z2为支路连接的节点号
	double R,X;                                                                                                                    //R、X分别为支路电阻和电抗
}branch[10];

void Read_data(double pos[])
{
	int i;
	ifstream in("Input.txt");
	for(i=0;i<100;i++)
	{
		in>>pos[i];
	}
	node_num=pos[0];
	branch_num=pos[1];
	pq_num=pos[2];
	pv_num=pos[3];
	g = pq_num + pv_num;
	Ac=pos[5];
	for(i=1;i<node_num;i++)                                                                                        //定义节点号、节点类型、节点有功无功等
	{
		Node[i].num=pos[i*6];
		Node[i].t=pos[i*6+1];
		Node[i].P=pos[i*6+2];
		Node[1].Q=pos[9];
		Node[2].Q = 0;
		Node[2].U = pos[15];
		Node[i].e=pos[i*6+4];
		Node[i].f=pos[i*6+5];
		Node[node_num].num=pos[node_num*6];
		Node[node_num].t=pos[node_num*6+1];
		Node[node_num].e=pos[node_num*6+2];
		Node[node_num].f=pos[node_num*6+3];
	}
	for(i=1;i<=branch_num;i++)
	{
		branch[1].numb=pos[22];
		branch[1].z1=pos[23];
		branch[1].z2=pos[24];
		branch[1].R=pos[25];
		branch[1].X=pos[26];
		branch[2].numb=pos[27];
		branch[2].z1=pos[28];
		branch[2].z2=pos[29];
		branch[2].R=pos[30];
		branch[2].X=pos[31];
	}
	in.close();
}

void Form_Y(double ad[],int node)
{
	out<<"--------------------潮流计算过程-------------------"<<endl;
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
	out<<"电导矩阵"<<endl;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			out<<G[i][j]<<" ";
		}
		out<<endl;
	}
	out<<"电纳矩阵"<<endl;
	for(i=0;i<node;i++)
	{
		for(j=0;j<node;j++)
		{
			out<<B[i][j]<<" ";
		}
		out<<endl;
	}
}

void mozhi(double a0,double b0)                      /*复数求模值函数*/
 { 	 
	 M=sqrt(a0*a0+b0*b0); 
 } 
void ji(double a1,double b1,double a2,double b2)     /*复数求积函数*/
 {   c1=a1*a2-b1*b2;
	 d1=a1*b2+a2*b1;
 }
void shang(double a3,double b3,double a4,double b4)  /*复数求商函数*/
 {   c2=(a3*a4+b3*b4)/(a4*a4+b4*b4);
     d2=(a4*b3-a3*b4)/(a4*a4+b4*b4);
 }

void Calculate_Unbalanced_Para()      //计算pq_num节点不平衡分量delta_P、delta_Q和PV节点的delta_P、deltaU^2
{
	int i,j,k,a;
	double Cc[10],Dd[10];
	for (i = 1; i <= node_num; i++)                       /*计算PQ节点不平衡量*/
	{
		if (Node[i].t == 1)
		{
			k = Node[i].num;
			Cc[k] = Dd[k] = 0;
			for (j = 1; j <= node_num; j++)
			{
				for (a = 1; a <= node_num; a++)
				{
					if (Node[a].num == j)
						break;
				}
				ji(G[k - 1][j - 1], -B[k - 1][j - 1], Node[a].e, -Node[a].f);
				Cc[k] += c1;
				Dd[k] += d1;
			}
			ji(Node[i].e, Node[i].f, Cc[k], Dd[k]);
			Node[i].delta_P = Node[i].P - c1;
			Node[i].delta_Q = Node[i].Q - d1;
		}
		if (Node[i].t == 2)                              /*计算PV节点不平衡量*/
		{
			k = Node[i].num;
			Cc[k] = Dd[k] = 0;
			for (j = 1; j <= node_num; j++)
			{
				for (a = 1; a <= node_num; a++)
				{
					if (Node[a].num == j)
						break;
				}
				ji(G[k - 1][j - 1], -B[k - 1][j - 1], Node[a].e, -Node[a].f);
				Cc[k] += c1;
				Dd[k] += d1;
			}
			ji(Node[i].e, Node[i].f, Cc[k], Dd[k]);
			Node[i].delta_P = Node[i].P - c1;
			Node[i].Q = d1;
			Node[i].delta_U = Node[i].U * Node[i].U - (Node[i].e * Node[i].e + Node[i].f * Node[i].f);
		}
	}
    for(i=1;i<=pq_num;i++)                                   /*形成不平衡量矩阵D[ ]*/
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
		out<<"delta_P["<<i<<"]="<<D[2*k-1]<<endl;
		out<<"delta_Q["<<i<<"]="<<D[2*k]<<endl;
	} 
	for(i=pq_num+1;i<=g;i++)
	{	
		k=Node[i].num;
		out<<"delta_P["<<i<<"]="<<D[2*k-1]<<endl;
		out<<"delta_U["<<i<<"]="<<D[2*k]<<endl;
	} 
}

void Form_J_Matrix()
{
	int i, j, k;
	double Aa[10], Bb[10];
	for (i = 1; i <= pq_num; i++)            /*形成PQ节点子阵*/
	{
		for (j = 1; j < node_num; j++)
		{
			int i1 = Node[i].num;
			int j1 = Node[j].num;
			double ei = Node[i].e;
			double ej = Node[j].e;
			double fi = Node[i].f;
			double fj = Node[j].f;

			if (i != j)                                              /*求i!=j时的H、N、J、L*/
			{
				Jacob[2 * i - 1][2 * j - 1] = -B[i1 - 1][j1 - 1] * ei + G[i1 - 1][j1 - 1] * fi;      /* H */
				Jacob[2 * i - 1][2 * j] = G[i1 - 1][j1 - 1] * ei + B[i1 - 1][j1 - 1] * fi;         /* N */
				Jacob[2 * i][2 * j - 1] = -G[i1 - 1][j1 - 1] * ei - B[i1 - 1][j1 - 1] * fi;        /* J */
				Jacob[2 * i][2 * j] = -B[i1 - 1][j1 - 1] * ei + G[i1 - 1][j1 - 1] * fi;          /* L */
			}                                                                                                                                                                                                                                                    
			else                                                   /*求i=j时的H、N、J、L*/
			{
				Aa[i] = 0; Bb[i] = 0;
				for (k = 1; k <= node_num; k++)
				{
					int k1 = Node[k].num;
					ji(G[i1 - 1][k1 - 1], B[i1 - 1][k1 - 1], Node[k].e, Node[k].f);
					Aa[i] += c1;
					Bb[i] += d1;
				}
				Jacob[2 * i - 1][2 * j - 1] = -B[i1 - 1][i1 - 1] * ei + G[i1 - 1][i1 - 1] * fi + Bb[i];      /*H*/
				Jacob[2 * i - 1][2 * j] = G[i1 - 1][i1 - 1] * ei + B[i1 - 1][i1 - 1] * fi + Aa[i];         /*N*/
				Jacob[2 * i][2 * j - 1] = -G[i1 - 1][i1 - 1] * ei - B[i1 - 1][i1 - 1] * fi + Aa[i];        /*J*/
				Jacob[2 * i][2 * j] = -B[i1 - 1][i1 - 1] * ei + G[i1 - 1][i1 - 1] * fi - Bb[i];          /*L*/
			}
		}
	}

	for (i = pq_num + 1; i <= g; i++)                                        /*形成PV节点子阵*/
		for (j = 1; j < node_num; j++)
		{
			int i1 = Node[i].num;
			int j1 = Node[j].num;
			double ei = Node[i].e;
			double ej = Node[j].e;
			double fi = Node[i].f;
			double fj = Node[j].f;
			if (i != j)                                           /*求i!=j时的H、N*/
			{
				Jacob[2 * i - 1][2 * j - 1] = -B[i1 - 1][j1 - 1] * ei + G[i1 - 1][j1 - 1] * fi;     /*H*/
				Jacob[2 * i - 1][2 * j] = G[i1 - 1][j1 - 1] * ei + B[i1 - 1][j1 - 1] * fi;        /*N*/
				Jacob[2 * i][2 * j - 1] = Jacob[2 * i][2 * j] = 0;                  /*R、S*/
			}

			else                                               /*求i=j时的H、N、R、S*/
			{
				Aa[i] = 0; Bb[i] = 0;
				for (k = 1; k <= node_num; k++)
				{
					int k1 = Node[k].num;
					ji(G[i1 - 1][k1 - 1], B[i1 - 1][k1 - 1], Node[k].e, Node[k].f);
					Aa[i] = Aa[i] + c1;
					Bb[i] = Bb[i] + d1;
				}
				Jacob[2 * i - 1][2 * j - 1] = -B[i1 - 1][i1 - 1] * ei + G[i1 - 1][i1 - 1] * fi + Bb[i];      /*H*/
				Jacob[2 * i - 1][2 * j] = G[i1 - 1][i1 - 1] * ei + B[i1 - 1][i1 - 1] * fi + Aa[i];         /*N*/
				Jacob[2 * i][2 * j - 1] = 2 * fi;                                /*R*/
				Jacob[2 * i][2 * j] = 2 * ei;                                  /*S*/
			}
		}
	out << "雅各比矩阵为：" << endl;
	for ( i = 1; i <= 2*g; i++)
	{
		for (j = 1; j <= 2 * g; j++)
		{
			out << Jacob[i][j] << " ";
		}
		out << endl;
	}
}

void Solve_Equations()                              /* 求解修正方程组 (用列主元消去法)*/
{
	double rr, tt, t;
	int i,j, k, l;
	for (i = 1; i <= 2 * g; i++)                             /*把函数残量矩阵编入修正方程*/
		Jacob[i][2 * g + 1] = D[i];
	k = 1;
	do
	{
		rr = Jacob[k][k];
		l = k;
		i = k + 1;
		do
		{
			if (fabs(Jacob[i][k]) > fabs(rr))                  /*列选主元*/
			{
				rr = Jacob[i][k];
				l = i;
			}
			i++;
		} while (i <= 2 * g);
		if (l != k)
		{
			for (j = k; j <= 2 * g + 1; j++)
			{
				tt = Jacob[l][j];
				Jacob[l][j] = Jacob[k][j];
				Jacob[k][j] = t;
			}
		}
		for (j = k + 1; j <= 2 * g + 1; j++)                           /*消元*/
			Jacob[k][j] /= Jacob[k][k];
		for (i = k + 1; i <= 2 * g; i++)
			for (j = k + 1; j <= 2 * g + 1; j++)
				Jacob[i][j] -= Jacob[i][k] * Jacob[k][j];
		k++;
	} while (k <= 2 * g);

	if (k != 1)
	{
		for (i = 2 * g; i >= 1; i--)                               /*回代*/
		{
			tt = 0;
			for (j = i + 1; j <= 2 * g; j++)
				tt += Jacob[i][j] * D[j];
			D[i] = Jacob[i][2 * g + 1] - tt;
		}
	}

	for (i = 1; i <= g; i++)
	{
		Node[i].e += D[2 * i];
		Node[i].f += D[2 * i - 1];
	}
	max1 = fabs(D[1]);
	for (i = 1; i <= 2 * (pq_num + pv_num); i++)
		if (fabs(D[i]) > max1)
			max1 = fabs(D[i]);

}

void Newton_Raphson()
{
	int z = 1,a,c;
	do
	{
		out<<"*****************************************************************"<<endl;
		max1 = 1;
		if ((z < 100) && (max1 >= Ac))
		{
			out << "迭代次数" << z-1 << endl;
		}

		Calculate_Unbalanced_Para();
		Form_J_Matrix();
		Solve_Equations();


		out << "输出delta_e、delta_f" << endl;
		for (c = 1; c <= node_num; c++)
		{
			for (a = 1; a <= node_num; a++)
			{
				if (Node[a].num == c)
					break;
			}
			out << "节点" << c << " " << "delta_f=" << D[2 * a - 1] << " " << "delta_e=" << D[2 * a] << endl;
		}
		out << "输出迭代过程中的电压值: " << endl;
		for (c = 1; c <= node_num; c++)
		{
			for (a = 1; a <= node_num; a++)
			{
				if (Node[a].num == c)
					break;
			}
			out << "U[" << c << "]=" << Node[a].e ;
			if (Node[a].f >= 0)
				out << "+j" << Node[a].f << endl;
			else
				out << "-j" << -Node[a].f << endl;
		}
		z++;
	} while ((z < 100) && (max1 >=Ac));

}

void Powerflow_Result()
{
	int n1=Node[node_num].num,a,b,c,i,j;
	double Aa[10],Bb[10];
	double rr,tt;
	out<<"--------------------潮流计算结果-------------------"<<endl;
	out<<"各节点电压："<<endl;                                                                           //计算各节点电压
        for(c=1;c<=node_num;c++)
		{
			for(a=1;a<=node_num;a++)		  
			{
				if(Node[a].num==c)
					break; 
			}   
			out<<"U["<<c<<"]="<<Node[a].e;
			if(Node[a].f>=0)
				out<<"+j"<<Node[a].f<<endl;
			else
				out<<"-j"<<-Node[a].f<<endl;
		}
	
	rr=tt=0;	
	for(i=1;i<=node_num;i++)
	{
		int i1=Node[i].num;
		ji(G[n1-1][i1-1],-B[n1-1][i1-1],Node[i].e,-Node[i].f);
		rr+=c1;
		tt+=d1;
	}
	ji(Node[node_num].e,Node[node_num].f,rr,tt);
	out<<"各节点注入功率："<<endl;                                                    //计算各节点注入功率
  for(i=1;i<=pq_num;i++)
  { 
	  int i1=Node[i].num;
	  out<<"PQ节点："<<"节点"<<i1<<"     "<<"S["<<i1<<"]="<<Node[i].P;
	  if(Node[i].Q>=0)
		  out<<"+j"<<Node[i].Q<<endl;
	  else
		  out<<"-j"<<-Node[i].Q<<endl;
  }
  for(i=pq_num+1;i<=g;i++)
  { 
	  int i1=Node[i].num;
	  out<<"PV节点："<<"节点"<<i1<<"     "<<"S["<<i1<<"]="<<Node[i].P;
	  if(Node[i].Q>=0)
		  out<<"+j"<<Node[i].Q<<endl;
	  else
		  out<<"-j"<<-Node[i].Q<<endl;
  }
  
  out<<"平衡节点:   "<<"节点"<<Node[node_num].num<<"     "<<"S["<<n1<<"]="<<c1;
	if(d1>=0)
		out<<"+j"<<d1<<endl;
	else
		out<<"-j"<<-d1<<endl;

	out<<"线路功率："<<endl;                                                            //计算线路功率
	rr=tt=0;
    for(i=1;i<=branch_num;i++)
	{
		int i1=branch[i].z1;
		int j1=branch[i].z2;
		Aa[i]=Bb[i]=0;
		for(a=1;a<=node_num;a++)			  
		{
			if(Node[a].num==i1)
				break;
		} 
		for(b=1;b<=node_num;b++)			  
		{
			if(Node[b].num==j1)
				break;
		}
		mozhi(branch[i].R,branch[i].X);
		if(M==0)
			continue;
		shang(1,0,branch[i].R,branch[i].X);
		ji(Node[a].e-Node[b].e,-Node[a].f+Node[b].f,c2,-d2);
		ji(Node[a].e,Node[a].f,c1,d1);
		out<<"线路"<<i<<"     "<<i1<<"--"<<j1<<"的功率为:   "<<c1;
		if(d1>=0)
			out<<"+j"<<d1<<endl;
	    else
			out<<"-j"<<-d1<<endl;
		Aa[i]+=c1;
		Bb[i]+=d1;
		ji(Node[b].e-Node[a].e,-Node[b].f+Node[a].f,c2,-d2);
		ji(Node[b].e,Node[b].f,c1,d1);
		out<<"线路"<<i<<"     "<<j1<<"--"<<i1<<"的功率为:   "<<c1;
		if(d1>=0)
			out<<"+j"<<d1<<endl;
	    else
			out<<"-j"<<d1<<endl;
        Aa[i]+=c1;
		Bb[i]+=d1;
		rr+=Aa[i];
		tt+=Bb[i];

	} 
	out<<"线路损耗功率:"<<endl;                                                                //计算线路损耗的功率
	for(i=1;i<=branch_num;i++)
	{
		int i1=branch[i].z1;
		int j1=branch[i].z2;
		out<<"线路"<<i<<"损耗的功率为:   "<<Aa[i];
		if(Bb[i]>=0)
			out<<"+j"<<Bb[i]<<endl;
	    else
			out<<"-j"<<-Bb[i]<<endl;
	}
			
	out<<"网络总损耗功率为:   "<<rr;                                                       //计算网络总损耗功率
		if(tt>=0)
			out<<"+j"<<tt<<endl;
	    else
			out<<"-j"<<-tt<<endl;
	out<<"********************************************************************"<<endl;
	out<<"----------潮流计算结束----------"<<endl;
}

void main()
{
	double ad[100];
	int i;
	for(i=0;i<100;i++)
	{
		ad[i]=0;
	}
	//初始化输入数组
	Read_data(ad);
	out<<"                       潮流上机计算"<<endl;
	out<<"题号：45       学号：201705000307      姓名：李腾      班级：电力英1703"<<endl;
	Form_Y(ad,node_num);
	Newton_Raphson();
	Powerflow_Result();
	out.close();
	cout<<"潮流计算成功！"<<endl;
	cout<<"请在Output.txt中查看结果。"<<endl;
	system("pause");
}