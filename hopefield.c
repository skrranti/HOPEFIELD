#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"


#define N 30
#define P 157
#define seme 3845
#define T 0.0001
#define montsteps 50
#define comienzo 1
#define def 1

gsl_rng *rn;

double minimo(double w, double v);
int swap(int elemento);

int main()
{
		FILE *f1, *f2, *f3, *f4, *f5, *f6;

		f1=fopen("memory.txt","r");
		f2=fopen("red.read_ini.txt","r");
		f3=fopen("red.write_ini.txt","w");
		f4=fopen("red.write_fin.txt","w");
		f5=fopen("solapamiento.txt","w");
		f6=fopen("prueba.txt","w");


		extern gsl_rng *rn;
		int semilla=seme;
		
		rn=gsl_rng_alloc(gsl_rng_taus);
		gsl_rng_set(rn,semilla);		

		int i,j,k,l,m;
		int x,y;
		double z;
		int crtl;
		double sum;
		
		int red[N][N];
		int chi[P][N][N];
		double a[P];
		double w[N][N][N][N];
		double theta[N][N];
		double solapamiento[P];
		int steps;

		double p;
		double DeltaE;
		double epsilon;

		int inicio=comienzo;

		steps=montsteps*N*N;

		


		for(i=0; i<(P); i++)  //Se generan P patrones aleatoriamente
		{
			for(j=0; j<N; j++)
			{
				for(k=0; k<N; k++)
				{
					z=gsl_rng_uniform(rn);			

					if(z>0.5) chi[i][j][k]=1;
					else chi[i][j][k]=0;
					fprintf(f6,"%i\t", chi[i][j][k]); //Aqui se escriben los patrones que se han generado
				}
				fprintf(f6,"\n");	
			}
			fprintf(f6,"\n");
		}


		for(i=0; i<P; i++)  //Se calculan las a que son los valores medios de cada patron
		{
			a[i]=0;	
			for(j=0; j<N; j++)
			{
				for(k=0; k<N; k++)
				{
					a[i]=a[i]+chi[i][j][k];
				}
			}
			a[i]=a[i]*1./(N*N);
		}


		for(i=0; i<N; i++)  //Se calculan las a que son los valores medios de cada patron
		{	
			for(j=0; j<N; j++)
			{
				for(k=0; k<N; k++)
				{
					for(l=0; l<N; l++)
					{	
						w[i][j][k][l]=0;	
						if((i!=k)||(j!=l))
						{
							for(m=0; m<P; m++)
							{
								w[i][j][k][l]=w[i][j][k][l]+((chi[m][i][j]-a[m])*(chi[m][k][l]-a[m]));	
							}
						}
						w[i][j][k][l]=w[i][j][k][l]*1./(N*N);
					}
				}
			}
		}


		for(i=0; i<N; i++)  //Se calculan las a que son los valores medios de cada patron
		{	
			for(j=0; j<N; j++)
			{
				theta[i][j]=0;
				for(k=0; k<N; k++)
				{
					for(l=0; l<N; l++)
					{	
						theta[i][j]=theta[i][j]+w[i][j][k][l];
					}
				}
				theta[i][j]=theta[i][j]*0.5;
			}
		}
		
		

		if(inicio==0)
		{
			for(i=0; i<N; i++)  //Inicializamos la red con valores aleatorios
			{
				for(j=0; j<N; j++)
				{
					z=gsl_rng_uniform(rn);			

					if(z>0.5) red[i][j]=1;
					else red[i][j]=0;
					fprintf(f3,"%i\t", red[i][j]);
				}
				fprintf(f3,"\n");
			}
		}
		
		else
		{
			for(i=0; i<N; i++)  //Inicializamos la red leyendo el primer patron
			{
				for(j=0; j<N; j++)
				{

					red[i][j]=chi[0][i][j]; 	
				}
			}


			for(i=0; i<(def*N); i++)  //Ese patrón que se ha leido se deforma ahora 
			{
				x=gsl_rng_uniform_int(rn,N);
				y=gsl_rng_uniform_int(rn,N);
				red[x][y]=swap(red[x][y]);
			}

			for(i=0; i<N; i++)  //Aqui se escribe el patron inicial que queda tras la deformacion
			{
				for(j=0; j<N; j++)
				{
					fprintf(f3,"%i\t", red[i][j]);	
				}
				fprintf(f3,"\n");
			}

			

		}




		for(j=0; j<1; j++)  //CALCULO DEL SOLAPAMIENTO INICIAL (solamente el del primer patron)
		{
			solapamiento[j]=0;
			for(k=0; k<N; k++)
			{
				for(l=0; l<N; l++)
				{	
					solapamiento[j]=solapamiento[j]+((chi[j][k][l]-a[j])*(red[k][l]-a[j]));	
				}
			}
			solapamiento[j]=solapamiento[j]*1./(N*N*a[j]*(1-a[j]));
			fprintf(f5, "%lf\t", solapamiento[j]);
		}
		fprintf(f5, "0\n");	




	
		for(i=0; i<steps; i++)
		{
			x=gsl_rng_uniform_int(rn,N);
			y=gsl_rng_uniform_int(rn,N);

			sum=0;
			for(k=0; k<N; k++)  //Inicializamos la red leyendo de un fichero
			{
				for(l=0; l<N; l++)
				{
					sum=sum+(red[k][l]*w[x][y][k][l]);	
				}
			}


			DeltaE=((red[x][y]-swap(red[x][y]))*sum)-((red[x][y]-swap(red[x][y]))*theta[x][y]);

			p=minimo(1,exp(-DeltaE*1./T));
			epsilon=gsl_rng_uniform(rn);		
			
			if(epsilon<p) red[x][y]=swap(red[x][y]);

			crtl=(i+1)%(N*N);		

			if(crtl==0)
			{
				for(j=0; j<1; j++)  
				{
					solapamiento[j]=0;
					for(k=0; k<N; k++)
					{
						for(l=0; l<N; l++)
						{	
							solapamiento[j]=solapamiento[j]+((chi[j][k][l]-a[j])*(red[k][l]-a[j]));	
						}
					}
					solapamiento[j]=solapamiento[j]*1./(N*N*a[j]*(1-a[j]));
					fprintf(f5, "%lf\t", solapamiento[j]);
				}
				fprintf(f5, "%i\n", (i+1)/(N*N));	
			}
			
		}

		for(k=0; k<N; k++)
		{
			for(l=0; l<N; l++)
			{
				fprintf(f4, "%i\t", red[k][l]);
			}
			fprintf(f4, "\n");
		}
		
		

		fclose(f1);
		fclose(f2);
		fclose(f3);
		fclose(f4);
		fclose(f5);

		return 0;
}

double minimo(double w, double v)
{
		double f;
		if(w<=v) f=w;
		else f=v;

		return f; 
}

int swap(int elemento)
{
	int f;
	if(elemento==0) f=1;
	else f=0;

	return f;
}
