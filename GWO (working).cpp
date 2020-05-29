
//Merhaba hocam programin çalismasi icin asagidaki ayari yapmaniz gerekmektedir.
//Tools -> Compiler Options... -> Settings -> Code Generation -> Language Standart -> ISO C++ 11 olmali


#include<iostream>

#include<string>

#include<algorithm>

#include<stdio.h>

#include<math.h>

#include<cstdint>

#include<numeric>

#include <iomanip>

#include <random>

#include<limits>

#include <fstream>


#include <bits/stdc++.h> 



using namespace std;
int len(int x[])
{
	return (sizeof(x)/sizeof(x[0]));
}
int lenDouble(double x[])
{
	return (sizeof(x)/sizeof(x[0]));
}
double *power(double x[], int y)
{
	int lx=lenDouble(x);
	double *xp= new double[lx]; 
	
	for(int i=0;i<lx;i++)
	{
		xp[i]=pow(x[i],y);
	}
	return (xp);
}
double prod(double a[])
{
	int l=lenDouble(a);
	double result;
	for (int i=0;i<l;i++)
	{
		result=a[i]*result;
	}
	return result;
}
double arraySum(double a[])  
{
	int n = lenDouble(a); 
    double sum  = 0;
    for(int i=0;i<n;i++)
    {
    	sum=sum+a[i];
    }  
    return (sum); 
}
double F1(double x[])
{
	double fit=arraySum(power(x,2));
	return fit;
}
double F2(double x[])
{
	
	double fit=arraySum(x)+prod(x);
	return fit;
}


int GWO(int lb,int ub,int dim,int SearchAgents_no,int Max_iter)
{
	double Alpha_pos[dim];
	double Alpha_score = std::numeric_limits<double>::infinity();


	double Beta_pos[dim];
	double Beta_score = std::numeric_limits<double>::infinity();

	double Delta_pos[dim];
	double Delta_score = std::numeric_limits<double>::infinity();
	
	double Positions[SearchAgents_no][dim]; 
	
	
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> distribution(lb, ub);
    double array[SearchAgents_no][dim];

    for(int i = 0; i < SearchAgents_no; i++)
    {
    	for (int j=0;j<dim;j++)
    	{
        	double d = distribution(mt);
        	array[i][j] = d;
        }
    }

	double productArray[SearchAgents_no][dim];
	for (int i=0;i<SearchAgents_no;i++)
	{
		for(int j=0;j<dim;j++)
		{
			productArray[i][j]=array[i][j]*(ub-lb);
		}
	}
	
	double addArray[SearchAgents_no][dim];
	for (int i=0; i<SearchAgents_no;i++)
	{
		for (int j=0;j<dim;j++)
		{
			addArray[i][j]=productArray[i][j]+lb;
		}
	}

	for (int i=0;i<SearchAgents_no;i++) 
	{																			
		for(int j=0;j<dim;j++)
		{
			Positions[i][j]=addArray[i][j];																
		}
	}																			

	double Convergence_curve[Max_iter];
	

	
	
	for(int l=0;l<Max_iter;l++)
	{
	
		for(int i=0;i<SearchAgents_no;i++)
		{	
			
			double fitness=F1(Positions[i]);
		
			if(fitness<Alpha_score)
			{
				Alpha_score=fitness;  
				
				for(int z=0;z<SearchAgents_no;z++)
				{
					Alpha_pos[z]=Positions[i][z];
				}
			}

			if(fitness>Alpha_score && fitness<Beta_score)
			{
				Beta_score=fitness;
			
				for(int z=0;z<SearchAgents_no;z++)
				{
					Beta_pos[z]=Positions[i][z];
				}
			}

			if(fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score)
			{
				Delta_score=fitness; 
				
				for(int z=0;z<SearchAgents_no;z++)
				{
					Delta_pos[z]=Positions[i][z];
				}
			}
		}
		int a=2-1*((2)/Max_iter); 
		
		for (int i=0;i<SearchAgents_no;i++)
		{


			for (int j=0;j<dim;j++)
			{
			
				int r1=rand()%2;  
				int r2=rand()%2;   
				int A1=2*a*r1-a; 
				int C1=2*r2; 
				int D_alpha=abs(C1*Alpha_pos[j]-Positions[i][j]);
				int X1=Alpha_pos[j]-A1*D_alpha; 

				r1=rand()%2;
				r2=rand()%2;
				int A2=2*a*r1-a; 
				int C2=2*r2; 
				int D_beta=abs(C2*Beta_pos[j]-Positions[i][j]); 
				int X2=Beta_pos[j]-A2*D_beta; 

				r1=rand()%2;
				r2=rand()%2;
				int A3=2*a*r1-a; 
				int C3=2*r2; 
				int D_delta=abs(C2*Delta_pos[j]-Positions[i][j]); 
				int X3=Delta_pos[j]-A3*D_delta; 

				Positions[i][j]=(X1+X2+X3)/3; 
			}
		}
		

		cout<< "At iteration " << l << "the best fitness is " << Alpha_score<<"\n";
	}
}

int main(int args, char* arg[])
{

	int NumOfRuns=2;
	int PopulationSize=50;
	int Iterations=100;

	for (int j=0;j<23;j++)
	{
		for(int k=0;k<NumOfRuns;k++)
		{
			GWO(-100,100,30,PopulationSize,Iterations);
		}
	}
}
