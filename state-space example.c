#include <stdio.h>
#include <stdlib.h>




int main()
{

    float x1,x1p,x2,x2p,x3,x3p,x4,x4p,y,u;
    int k;
    int N=30;

    //Open an Input File
    FILE *inputFile;
    inputFile=fopen("input.txt", "r");

    //Open an Output File
    FILE *outputFile;
    outputFile=fopen("output.txt","w");

    //Initialize the State Variable to Zero
    x1=0;
    x2=0;
    x3=0;
    x4=0;


    for(k=0;k<N;k++)
    {
        fscanf(inputFile,"%f,",&u);
        //State-Update Equations
		x1p=0.9499*x1+0.04910*x2+0*x3+0*x4+(-7.1089)*u;
		x2p=(-0.0491)*x1+0.9499*x2+0*x3+0*x4+64.5389*u;
		x3p=0*x1+0*x2+0.3001*x3+0.4*x4+4.4020*u;
		x4p=0*x1+0*x2+(-0.4)*x3+0.3001*x4+2.2418*u;
		//Output Equation
		y=(-0.02)*x1+0.0146*x2+0.111*x3+(-.2568)*x4+0*u;

		//Remember the New States (move from k to k+1)
		x1=x1p;
		x2=x2p;
		x3=x3p;
		x4=x4p;

        fprintf(outputFile,"%f, ", y);
        printf("%5.4f   %5.4f\n",u,y);
    }


    fclose(inputFile);
    fclose(outputFile);
}
