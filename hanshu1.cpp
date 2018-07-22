#include<iostream>
using namespace std;
int gy(int a,int b)
{
	int c,t;
	if(a>b)
	{
		t=a;
		a=b;
		b=t;
	}
	for(c=a;c>0;c--)
		if(a%c==0&&b%c==0)
			break;
	return c;
}
int gb(int a,int b)
{
	int c,t;
    if(a>b)
	{
		t=a;
		a=b;
		b=t;
	}
	for(c=b;c<=a*b;c++)
		if(c%a==0&&c%b==0)
			break;
	return c;
}
int main()
{
	int a,b,c,d;
	cin>>a>>b;
	c=gy(a,b);
	d=gb(a,b);
	cout<<c<<" "<<d<<endl;
	return 0;
}

