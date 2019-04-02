#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Nx 256 // Numero de cuadro en la direccion x
#define Ny 256 // Numero de cuadro en la direccion y
#define Nx1 (Nx+1)
#define Ny1 (Ny+1)
#define L (Ny+1) // Ancho de la cavidad
#define Q 9 		// Numero de velocidades discretas
#define rho0 1.0  // Densidad Inicial
#define ux0 0.0   // Velocidad inicial en la componente x
#define uy0 0.0   // Velocidad inicial en la componente y 
#define uw 0.1
#define Re 400.0

int cx[Q]={0, 1, 0, -1, 0, 1, -1, -1, 1};
int cy[Q]={0, 0, 1, 0, -1, 1, 1, -1, -1};
double f[Ny1][Nx1][Q]; //Arreglo de las funciones de distribucion
double f_post[Ny1][Nx1][Q]; // Arreglo de las funciones de distribucion luego de la colision
double rho[Ny1][Nx1], ux[Ny1][Nx1], uy[Ny1][Nx1];
// Arreglo de la densidad, velocidad en x e y 
double tau; // Tiempo de relajacion en el modelo BGK
double w[Q]={4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/36 ,1.0/36,
1.0/36,1.0/36}; // Pesos 
//int rc[Q]={0,3,4,1,2,7,8,5,6}; // index of reversed velocity
void Init_Eq(void); //Funcion de Initializacion
double feq(double RHO, double U, double V, int k);
// Funcionde distribucion de equilibrio 
void Coll_BGK(void); // BGK colision
void Streaming(void); // Streaming (Transmision)
void Den_Vel(void); // Variables macroscopicas
void Bounce_back(void); // Bounce-back Condiciones de frontera
double Err(void); // Funcion del error para parar el ciclo 
double u0[Ny1][Nx1],v0[Ny1][Nx1]; // definicion de matrices de condiciones iniciales
void Data_Output(void); // Funcion que escribe los datos

//=========================================================
//=========================================================
int main()
{
	int k,M2,N2;
	double err;
	M2=Ny/2; N2=Nx/2;
	k=0;
	err=1.0;
	tau=3*L*uw/Re+0.5; // Tiempo de relajacion para BGK
	Init_Eq();
	while(err>1.0e-2)
	{
		k++;
		Coll_BGK(); //BGK colision
		Streaming(); // Streaming (Transmision)
		Bounce_back(); // Condiciones de frontera
		Den_Vel(); // Variables macroscopicas de fluido 
		if(k%1000==0)
		{
			err=Err(); // Diferencia entre las velocidades cada mil pasos
			printf("err=%e ux_center=%e uy_center=%e k=%d\n",err,ux[M2][N2],uy[M2][N2], k); 
		}
	}
	Data_Output(); //Escribir los pasos cuando acabe la iteracion
}


// Funcion que inicializa las matrices
void Init_Eq()
{
	int j, i, k;
	for(j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)
	{
		rho[j][i]=rho0;
		ux[j][i]=ux0;
		uy[j][i]=uy0;
		for(k=0;k<Q;k++)
		f[j][i][k]=feq(rho[j][i],ux[j][i],uy[j][i],k);
	}
}


// Calculo de la distribucion de equilibrio 
double feq(double RHO, double U, double V, int k)
{
	double cu, U2;
	cu=cx[k]*U+cy[k]*V; // c k*u
	U2=U*U+V*V; // u*u; norma al cuadrado
	return w[k]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


// Funcion que hace la colision BGK
void Coll_BGK()
{
	int j, i, k;
	double FEQ;
	for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++) for(k=0;k<Q;k++)
	{
		FEQ=feq(rho[j][i],ux[j][i],uy[j][i],k); // EDF
		f_post[j][i][k] = f[j][i][k]-(f[j][i][k]-FEQ)/tau;// Post-collision funciones de distribucion
	}
}


// Streaming (Transmision de la informacion)
void Streaming()
{
	int j, i, jd, id, k;
	for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++) for(k=0;k<Q;k++)
	{
	jd=j-cy[k]; id=i-cx[k]; 
	if(jd>=0 && jd<=Ny && id>=0 && id<=Nx) 
		f[j][i][k]=f_post[jd][id][k]; // streaming
	}
}


//Bounce Back
void Bounce_back()
{
	int i,j;
	// j=Ny: Parte superior
	for(i=0;i<=Nx;i++)
	{
		f[Ny][i][4]=f_post[Ny][i][2];
		f[Ny][i][7]=f_post[Ny][i][5]+6*rho[Ny][i]*w[7]*cx[7]*uw;
		f[Ny][i][8]=f_post[Ny][i][6]+6*rho[Ny][i]*w[8]*cx[8]*uw;
	}
	// j=0: Parte inferior
	for(i=0;i<=Nx;i++)
	{
		f[0][i][2]=f_post[0][i][4];
		f[0][i][5]=f_post[0][i][7];
		f[0][i][6]=f_post[0][i][8];
	}
	// i=0: Pared izquierda
	

	for(j=0;j<=Ny;j++)
	{
		f[j][0][1]=f_post[j][0][3];
		f[j][0][5]=f_post[j][0][7];
		f[j][0][8]=f_post[j][0][6];
	}
	// i=Nx: Pared derecha
	for(j=0;j<=Ny;j++)
	{
		f[j][Nx][3]=f_post[j][Nx][1];
		f[j][Nx][7]=f_post[j][Nx][5];
		f[j][Nx][6]=f_post[j][Nx][8];
	}
}



//Calculo de las variables macroscopicas
void Den_Vel()
{
	int j, i;
	for(j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)
	{
	rho[j][i]=f[j][i][0]+f[j][i][1]+f[j][i][2]+f[j][i][3]
	+f[j][i][4]+f[j][i][5]+f[j][i][6]+f[j][i][7]+
	f[j][i][8];
	ux[j][i]=(f[j][i][1]+f[j][i][5]+f[j][i][8]-f[j][i][3]-
	f[j][i][6]-f[j][i][7])/rho[j][i];
	uy[j][i]=(f[j][i][5]+f[j][i][6]+f[j][i][2]-f[j][i][7]-
	f[j][i][8]-f[j][i][4])/rho[j][i];
	}
}


//Error
//=========================================================
double Err() // Calculo del error relativo cada mil pasos
{
	int j, i;
	double e1,e2;
	e1=e2=0.0;
	for(j=1;j<Ny;j++) for(i=0;i<Nx;i++)
	{
		e1+=sqrt((ux[j][i]-u0[j][i])*(ux[j][i]-u0[j][i])
		+(uy[j][i]-v0[j][i])*(uy[j][i]-v0[j][i]));
		e2+=sqrt(ux[j][i]*ux[j][i]+uy[j][i]*uy[j][i]);
		u0[j][i]=ux[j][i];v0[j][i]=uy[j][i];
	}
	return e1/e2;
}



void Data_Output() // Datos de salida
{
	int i,j;
	FILE *fp;
	fp=fopen("pruebax.dat","w+");
	for(i=0;i<=Nx;i++) fprintf(fp,"%e \n", (i+0.5)/L);
	fclose(fp);
	fp=fopen("pruebay.dat","w+");
	for(j=0;j<=Ny;j++) fprintf(fp,"%e \n", (j+0.5)/L);
	fclose(fp);
	fp=fopen("pruebaux.dat","w");
	for(j=0;j<=Ny;j++) {
		for (i=0; i<=Nx; i++) 
		{
			fprintf(fp,"%e ",ux[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("pruebauy.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",uy[j][i]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("pruebarho.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",rho[j][i]);
	fprintf(fp,"\n");
	}
	fclose(fp);
}