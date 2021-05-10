/*
 * TestMat3.cpp
 *
 *  Created on: Nov 5, 2020
 *      Author: shadowcage72
 */






#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <Eigen/Dense>
#include <iomanip>
using namespace std;
using namespace Eigen;


int index2(double y,double z,int NG)
{
	if(z>NG-1)
	{
		z-=NG;
	}
	if(z<0)
	{
		z+=NG;
	}
	if(y>NG-1)
	{
		y-=NG;
	}
	if(y<0)
	{
		y+=NG;
	}
	return y*NG+z;
}








int main()
{
	string path="C:\\Users\\shadowcage72\\Eclipse\\TestMat3\\Data Numbered\\";


	double starttime = clock();
	double e0=1;//Epsilon Naught
	double k=.000000000168637205;//Boltzman
	double meterconv=0.00000000000003511282;//Convert from our units to meters
	double L=9/meterconv;//Length
	double DT=L/120;//L/60;//DeltaTs

	for(int number=37;number++;number<40)
	{
	double separation=.32*L+.04*L*(number-34);


	int NT=600000;//Iterations
	int NG=64;//Number of Grid Points
	double NG1=NG;
	double ndensity=10000000*meterconv*meterconv*meterconv;
	double totalparticles=2*ndensity*L*L;//Total number of electrons in simulation


	//At least 25 per cell is ideal

	int N=150000;//Number of Particles
	double me=totalparticles/N;//Mass of electron
	double mp=1*totalparticles/N;//mass of proton
	double e=totalparticles/N;//electron charge
	double Nin=N;
	int Nin1=N;
	double dx=L/NG;//Delta X
	double acc=1;//Number close to 1 as to make the Poisson Matrix nonsingular

	double dE=0;//Initialization of Energy Graph Value
	double dKE=0;
	double dUE=0;
	double dxmom=0;
	double dymom=0;
	double lostenergy=0;
	double esphere1=0;
	double esphere2=0;



	// initial loading for Particle Positions
	VectorXd xp(2*N);
	VectorXd yp(2*N);
	int n;
	for(n=0;n<N;n++)
	{
	xp(n)=L*((double) rand() / (RAND_MAX));
	yp(n)=L*((double) rand() / (RAND_MAX));

	xp(n+N)=0;
	yp(n+N)=0;
	};

	//Mass Matrix
	VectorXd massmat(2*N);
	for(n=0;n<N/2;n++)
		{
	massmat(n)=me;//Make all into ions
	massmat(n+N/2)=mp;
	massmat(n+N)=0;
	massmat(n+N+N/2)=0;
		}

	//Charge Matrix
	VectorXd emat(2*N);
		for(n=0;n<N/2;n++)
			{
		emat(n)=-e;//Make all into ions
		emat(n+N/2)=-e;
		emat(n+N)=0;
		emat(n+N+N/2)=0;
			}

	double Tempin=90*11700;//Initial Temperature
	double Temp=Tempin;//Temperature that will change in time

	VectorXd vxp(2*N);
	VectorXd vyp(2*N);
	double pi=3.141592653589793;
	double vrand;
	double anglerand;
	VectorXd InitialVelocity(2*N);
	VectorXd Pot(2*N);
	for(n=0;n<N;n++)
	{
		vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
		while(isinf(vrand))
		{
			vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
		}
		anglerand=2*pi*((double) rand() / (RAND_MAX));
		vxp(n)=vrand*cos(anglerand);
		vxp(n+N)=0;
		vyp(n)=vrand*sin(anglerand);//Make all into ions
		vyp(n+N)=0;
		InitialVelocity(n)=vrand;
	}
MatrixXd invPoi(NG*NG,NG*NG);
	for(int n=0;n<NG*NG;n++)
	{
		for(int m=0;m<NG*NG;m++)
		{
			invPoi(n,m)=0;
		}
	}
	for(n=0;n<NG*NG;n++)
	{
		invPoi(n,n)=-4.0001;
		if(n>0)
		{
			invPoi(n,n-1)=1;
		}
		if(n<NG*NG-2)
		{
			invPoi(n,n+1)=1;
		}
		if(n>NG-1)
		{
			invPoi(n,n-NG)=1;
		}
		if(n<NG*NG-NG-1)
		{
			invPoi(n,n+NG)=1;
		}
		if(n<NG-1)
		{
			invPoi(n,NG*NG-NG+n)=1;
		}
		if(n<NG)
		{
			invPoi(NG*NG-NG+n,n)=1;
		}
	}
	for(n=1;n<NG;n++)
	{
		invPoi(NG*n-1,NG*n)=0;
		invPoi(NG*n,NG*n-1)=0;

		invPoi(NG*n-1,NG*n-NG)=1;
		invPoi(NG*n-NG,NG*n-1)=1;
	}
	invPoi(NG*NG-1,NG*NG-NG)=1;
	invPoi(NG-1,NG*NG-1)=1;
	invPoi(NG*NG-2,NG*NG-1)=1;
	invPoi(NG*NG-NG,NG*NG-1)=1;
	invPoi(NG*NG-1-NG,NG*NG-1)=1;
//	invPoi(NG*NG-1,NG-1)=0;
//	invPoi(NG*NG-1,NG*NG-2)=0;
//	invPoi(NG*NG-1,NG*NG-NG)=0;
//	invPoi(NG*NG-1,NG*NG-1-NG)=0;


cout << "\nPreparing PIC Simulation";
invPoi=invPoi.inverse();

//Preallocate memory
	VectorXd Ex(NG*NG);
	VectorXd Ey(NG*NG);
	VectorXd xmom(NT);
	VectorXd ymom(NT,1);
	VectorXd ForceXSphere1(NT,1);
	VectorXd ForceYSphere1(NT,1);
	VectorXd ForceXSphere2(NT,1);
	VectorXd ForceYSphere2(NT,1);
	VectorXd totEn(NT,1);
	VectorXd totUEn(NT,1);
	VectorXd totKEn(NT,1);
	VectorXd totMomX(NT,1);
	VectorXd totMomY(NT,1);
	VectorXd UdUE(NT,1);
	VectorXd LdUE(NT,1);
	VectorXd Volt1(NT,1);
	VectorXd Volt2(NT,1);
	VectorXd Charge2(NT,1);
	VectorXd Charge1(NT,1);
	MatrixXd rho1allit(NG,NG);
	MatrixXd Phi1allit(NG,NG);
	for(n=0;n<NG;n++)
	{
		for(int m=0;m<NG;m++)
		{
			rho1allit(n,m)=0;
			Phi1allit(n,m)=0;
		}
	}
	VectorXd TotalMomentumAbsorbed1X(NT,1);
	VectorXd TotalMomentumAbsorbed1Y(NT,1);
	VectorXd TotalMomentumAbsorbed2X(NT,1);
	VectorXd TotalMomentumAbsorbed2Y(NT,1);
	VectorXd ElectricForceSphere1X(NT,1);
	VectorXd ElectricForceSphere1Y(NT,1);
	VectorXd ElectricForceSphere2X(NT,1);
	VectorXd ElectricForceSphere2Y(NT,1);
	VectorXd Temp1(NT,1);
	VectorXd TotalCharge1(NT,1);
	VectorXd totEnWLost(NT,1);
	VectorXd ForceYmean(NT,1);
	MatrixXd rho1(NG,NG);
	MatrixXd Phi1(NG,NG);
	MatrixXd Ex1(NG,NG);
	MatrixXd Ey1(NG,NG);


	double LD =sqrt(k*Temp/ndensity);


	//sphere1
	double sxc1=.5*L;//Sphere 1 x center
	double syc1=.25*L;//Sphere 1 y center
	double sr1=.02*L;//Sphere 1 radius


	//sphere2
	double sxc2=.5*L;//Sphere 2 x center
	double syc2=.75*L;//Sphere 2 y center
	double sr2=.02*L;//Sphere 2 radius

	for(n=0;n<N;n++)
	{
		if(((xp(n)-sxc1)*(xp(n)-sxc1)+(yp(n)-syc1)*(yp(n)-syc1))<sr1*sr1)
		{
			emat(n)=0;
		}
		if(((xp(n)-sxc2)*(xp(n)-sxc2)+(yp(n)-syc2)*(yp(n)-syc2))<sr2*sr2)
		{
			emat(n)=0;
		}
	}

	//Capacity Matrix

	//Put a particle on every grid point
	int i=0;
	int j=0;
	VectorXd xpcap(NG*NG);
	VectorXd ypcap(NG*NG);
	for(int n=0;n<NG*NG;n++)
	{
		xpcap(n)=dx*(j+.5);
		ypcap(n)=dx*(i+.5);
		if(i==NG-1)
		{
			i=0;
			j++;
		}
		else
		{
			i++;
		}
	}

	//Now a particle is on each point




	VectorXd ematc1(NG*NG);
	VectorXd ematc2(NG*NG);
	VectorXd ematc3(NG*NG);
	int lengthin1=0;
	int lengthin2=0;
	int lengthin3=0;
	for(n=0;n<NG*NG;n++)
	{
		ematc3(n)=0;
		if(((xpcap(n)-sxc1)*(xpcap(n)-sxc1)+(ypcap(n)-syc1)*(ypcap(n)-syc1))<sr1*sr1)
		{
			ematc1(n)=1;
			ematc3(n)=1;
			lengthin1++;
			lengthin3++;
		}
		else
		{
			ematc1(n)=0;
		}

		if(((xpcap(n)-sxc2)*(xpcap(n)-sxc2)+(ypcap(n)-syc2)*(ypcap(n)-syc2))<sr2*sr2)
		{
			ematc2(n)=1;
			ematc3(n)=1;
			lengthin2++;
			lengthin3++;
		}
		else
		{
			ematc2(n)=0;
		}
	}
//	cout <<"\nlengthin2="<<lengthin2;

	VectorXd rhotemp1(NG*NG);
	VectorXd Phitemp1(NG*NG);
	MatrixXd Bprime1(lengthin1,lengthin1);
	MatrixXd Bprime2(lengthin2,lengthin2);
	MatrixXd Bprime3(lengthin3,lengthin3);
	for(n=0;n<NG*NG;n++)
	{
		rhotemp1(n)=0;
	}
	i=0;
	j=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc1(n)==1)
		{
			rhotemp1(n)=1/dx/dx;
			Phitemp1=-invPoi*rhotemp1*dx*dx/e0;
			for (int m=0;m<NG*NG;m++)
			{
				if(ematc1(m)==1)
				{
					Bprime1(i,j)=Phitemp1(m);
					i++;
				}
			}
			i=0;
			j++;
			rhotemp1(n)=0;
		}
	}
	MatrixXd Cap1=Bprime1.inverse();
	i=0;
	j=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc2(n)==1)
		{
			rhotemp1(n)=1/dx/dx;
			Phitemp1=-invPoi*rhotemp1*dx*dx/e0;
			for (int m=0;m<NG*NG;m++)
			{
				if(ematc2(m)==1)
				{
					Bprime2(i,j)=Phitemp1(m);
					i++;
				}
			}
			i=0;
			j++;
			rhotemp1(n)=0;
		}
	}
	MatrixXd Cap2=Bprime2.inverse();

	i=0;
	j=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc3(n)==1)
		{
			rhotemp1(n)=1/dx/dx;
			Phitemp1=-invPoi*rhotemp1*dx*dx/e0;
			for (int m=0;m<NG*NG;m++)
			{
				if(ematc3(m)==1)
				{
					Bprime3(i,j)=Phitemp1(m);
					i++;
				}
			}
			i=0;
			j++;
			rhotemp1(n)=0;
		}
	}
	MatrixXd Cap3=Bprime3.inverse();



	//Circular Matrix representing 1 particle around Sphere 1
		int NCirc=64;//Number of particles for circle
		double NCirc1=NCirc;
		VectorXd xpc(NCirc);
		VectorXd ypc(NCirc);
		for (n=0;n<NCirc;n++)
		{
			xpc(n)=sxc1+sr1*sin(2*pi/NCirc*n);
			ypc(n)=syc1+sr1*cos(2*pi/NCirc*n);
		}
		double charge=1/NCirc1;

		MatrixXd gc(2*NCirc,2);
		for (n=0;n<NCirc;n++)
		{
			gc(n,0)=floor(xpc(n)/dx-.5);
			gc(n,1)=floor(ypc(n)/dx-.5);
			gc(NCirc+n,0)=(floor(xpc(n)/dx-.5)+1);
			gc(NCirc+n,1)=(floor(ypc(n)/dx-.5)+1);
		}


		MatrixXd frazc(2*NCirc,2);
		for (n=0;n<NCirc;n++)
		{
			frazc(n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(1-(ypc(n)/dx-.5-gc(n,1)));
			frazc(n,1)=(xpc(n)/dx-.5-gc(n,0))*(1-(ypc(n)/dx-.5-gc(n,1)));
			frazc(NCirc+n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(ypc(n)/dx-.5-gc(n,1));
			frazc(NCirc+n,1)=(xpc(n)/dx-.5-gc(n,0))*(ypc(n)/dx-.5-gc(n,1));
		}


		//Make the 1 particle sphere matrix
		VectorXd rho_sphere1(NG*NG);
		for (n=0;n<NG*NG;n++)
		{
			rho_sphere1(n)=0;
		}
		for(n=0;n<NCirc;n++)
		{
			rho_sphere1(index2(gc(n,0),gc(n,1),NG))+=charge*frazc(n,0)/dx/dx;
			rho_sphere1(index2(gc(n,0),gc(NCirc+n,1),NG))+=charge*frazc(n,1)/dx/dx;
			rho_sphere1(index2(gc(NCirc+n,0),gc(n,1),NG))+=charge*frazc(NCirc+n,0)/dx/dx;
			rho_sphere1(index2(gc(NCirc+n,0),gc(NCirc+n,1),NG))+=charge*frazc(NCirc+n,1)/dx/dx;
		}



		//Circular Matrix representing 1 particle around Sphere 2
			for (n=0;n<NCirc;n++)
			{
				xpc(n)=sxc2+sr2*sin(2*pi/NCirc*n);
				ypc(n)=syc2+sr2*cos(2*pi/NCirc*n);
			}

			for (n=0;n<NCirc;n++)
			{
				gc(n,0)=floor(xpc(n)/dx-.5);
				gc(n,1)=floor(ypc(n)/dx-.5);
				gc(NCirc+n,0)=(floor(xpc(n)/dx-.5)+1);
				gc(NCirc+n,1)=(floor(ypc(n)/dx-.5)+1);
			}


			for (n=0;n<NCirc;n++)
			{
				frazc(n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(1-(ypc(n)/dx-.5-gc(n,1)));
				frazc(n,1)=(xpc(n)/dx-.5-gc(n,0))*(1-(ypc(n)/dx-.5-gc(n,1)));
				frazc(NCirc+n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(ypc(n)/dx-.5-gc(n,1));
				frazc(NCirc+n,1)=(xpc(n)/dx-.5-gc(n,0))*(ypc(n)/dx-.5-gc(n,1));
			}

			//Make the 1 particle sphere matrix
					VectorXd rho_sphere2(NG*NG);
					for (n=0;n<NG*NG;n++)
					{
						rho_sphere2(n)=0;
					}
					for(n=0;n<NCirc;n++)
					{
						rho_sphere2(index2(gc(n,0),gc(n,1),NG))+=charge*frazc(n,0)/dx/dx;
						rho_sphere2(index2(gc(n,0),gc(NCirc+n,1),NG))+=charge*frazc(n,1)/dx/dx;
						rho_sphere2(index2(gc(NCirc+n,0),gc(n,1),NG))+=charge*frazc(NCirc+n,0)/dx/dx;
						rho_sphere2(index2(gc(NCirc+n,0),gc(NCirc+n,1),NG))+=charge*frazc(NCirc+n,1)/dx/dx;
					}



//cout <<"lengthin1="<<lengthin1;
//
	//color
	VectorXd c(2*N);
	for (n=0;n<N/2;n++)
		{
		c(n)=1;
		c(n+N/2)=10;
		c(n+N)=0;
		c(n+N+N/2)=0;
		}

cout << "\nOpening Files";

	ofstream myfile1;
		myfile1.open(path+to_string(number)+"Volt1.txt");

		ofstream myfile2;
		myfile2.open(path+to_string(number)+"Volt2.txt");

		ofstream myfile3;
		myfile3.open(path+to_string(number)+"Charge1.txt");

		ofstream myfile4;
		myfile4.open(path+to_string(number)+"Charge2.txt");


		ofstream myfile5;
		myfile5.open(path+to_string(number)+"rho1allit.txt");
		myfile5 << rho1allit;
		myfile5.close();

		ofstream myfile6;
		myfile6.open(path+to_string(number)+"Phi1allit.txt");
		myfile6 << Phi1allit;
		myfile6.close();

		ofstream myfile7;
		myfile7.open(path+to_string(number)+"Temp1.txt");

		ofstream myfile8;
		myfile8.open(path+to_string(number)+"velocitysquared.txt");
		myfile8 << vxp.array().square()+vyp.array().square();
		myfile8.close();

		ofstream myfile9;
		myfile9.open(path+to_string(number)+"totUEn.txt");

		ofstream myfile10;
		myfile10.open(path+to_string(number)+"totKEn.txt");

		ofstream myfile11;
		myfile11.open(path+to_string(number)+"totEn.txt");

		ofstream myfile12;
		myfile12.open(path+to_string(number)+"Phi1.txt");
		myfile12 << Phi1;
		myfile12.close();

		ofstream myfile13;
		myfile13.open(path+to_string(number)+"rho1.txt");
		myfile13 << rho1;
		myfile13.close();

		ofstream myfile14;
		myfile14.open(path+to_string(number)+"emat.txt");
		myfile14 << emat;
		myfile14.close();
//		cout <<"\nCheckpoint6";
		ofstream myfile15;
		myfile15.open(path+to_string(number)+"xp.txt");
		myfile15 << xp;
		myfile15.close();
//		cout <<"\nCheckpoint7";
		ofstream myfile16;
		myfile16.open(path+to_string(number)+"yp.txt");
		myfile16 << yp;
		myfile16.close();
//		cout <<"\nCheckpoint8";
		ofstream myfile17;
		myfile17.open(path+to_string(number)+"vxp.txt");
		myfile17 << vxp;
		myfile17.close();
//		cout <<"\nCheckpoint9";
		ofstream myfile18;
		myfile18.open(path+to_string(number)+"vyp.txt");
		myfile18 << vyp;
		myfile18.close();
//		cout <<"\nCheckpoint10";
		ofstream myfile19;
		myfile19.open(path+to_string(number)+"invPoi88.txt");
		//myfile19 << invPoi;
		myfile19.close();
//		cout <<"\nCheckpoint11";
		ofstream myfile20;
		myfile20.open(path+to_string(number)+"Cap1.txt");
		myfile20 << Cap1;
		myfile20.close();
//		cout <<"\nCheckpoint12";
		ofstream myfile21;
		myfile21.open(path+to_string(number)+"Cap2.txt");
		myfile21 << Cap2;
		myfile21.close();
//		cout <<"\nCheckpoint13";
		ofstream myfile22;
		myfile22.open(path+to_string(number)+"Cap3.txt");
		myfile22 << Cap3;
		myfile22.close();
//		cout <<"\nCheckpoint14";
		ofstream myfile23;
		myfile23.open(path+to_string(number)+"TotalMomentumAbsorbed1X.txt");
		ofstream myfile24;
		myfile24.open(path+to_string(number)+"TotalMomentumAbsorbed1Y.txt");
		ofstream myfile25;
		myfile25.open(path+to_string(number)+"TotalMomentumAbsorbed2X.txt");
		ofstream myfile26;
		myfile26.open(path+to_string(number)+"TotalMomentumAbsorbed2Y.txt");
		ofstream myfile27;
		myfile27.open(path+to_string(number)+"ElectricForceSphere1X.txt");
		ofstream myfile28;
		myfile28.open(path+to_string(number)+"ElectricForceSphere1Y.txt");
		ofstream myfile29;
		myfile29.open(path+to_string(number)+"ElectricForceSphere2X.txt");
		ofstream myfile30;
		myfile30.open(path+to_string(number)+"ElectricForceSphere2Y.txt");
		ofstream myfile31;
		myfile31.open(path+to_string(number)+"TotalMomentumX.txt");
		ofstream myfile32;
		myfile32.open(path+to_string(number)+"TotalMomentumY.txt");

		  if (myfile30.is_open())
		  {
		    cout << "\nStarting Simulation";
		  }
		  else
		  {
			  cout <<"\nFailed to Open Files";
		  }
	// Main computational cycle
	for(int it=0;it<NT;it++)
	{

		  // Get the system time.
		   unsigned seed = time(0);

		   // Seed the random number generator.
		   srand(seed);


	// update xp
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
					{
						xp(n)+=vxp(n)*DT;
						yp(n)+=vyp(n)*DT;
					}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
					{
						xp(n)+=vxp(n)*DT;
						yp(n)+=vyp(n)*DT;
					}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
					{
						xp(n)+=vxp(n)*DT;
						yp(n)+=vyp(n)*DT;
					}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
					{
						xp(n)+=vxp(n)*DT;
						yp(n)+=vyp(n)*DT;
					}
			}
		}


	//Capacity Matrix

	//Particle Delete if in a bordering grid space or out of bounds
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
				{
					if(xp(n)<0)
					{
						xp(n)+=L;
					}
					if(xp(n)>L)
					{
						xp(n)-=L;
					}
					if(yp(n)<0)
					{
						yp(n)+=L;
					}
					if(yp(n)>L)
					{
						yp(n)-=L;
					}
				}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
				{
					if(xp(n)<0)
					{
						xp(n)+=L;
					}
					if(xp(n)>L)
					{
						xp(n)-=L;
					}
					if(yp(n)<0)
					{
						yp(n)+=L;
					}
					if(yp(n)>L)
					{
						yp(n)-=L;
					}
				}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
				{
					if(xp(n)<0)
					{
						xp(n)+=L;
					}
					if(xp(n)>L)
					{
						xp(n)-=L;
					}
					if(yp(n)<0)
					{
						yp(n)+=L;
					}
					if(yp(n)>L)
					{
						yp(n)-=L;
					}
				}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
				{
					if(xp(n)<0)
					{
						xp(n)+=L;
					}
					if(xp(n)>L)
					{
						xp(n)-=L;
					}
					if(yp(n)<0)
					{
						yp(n)+=L;
					}
					if(yp(n)>L)
					{
						yp(n)-=L;
					}
				}
			}
		}



	cout  << "\nN="<<N;




	//Add Fresh Plasma
	double NGnew=4*NG-4;//Border Grid spaces getting replaced
	double Nnew=floor(NGnew*Nin/NG1/NG1+.5);//Put in a number of particles to keep N about the same
	double randmat;//Which side will they appear?




	MatrixXd g(4*Nin1,2);
#pragma omp parallel sections num_threads(4)
				{
		#pragma omp section
					{
						for(n=0;n<N/4;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
		#pragma omp section
					{
						for(n=N/4;n<N/2;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
		#pragma omp section
					{
						for(n=N/2;n<N/4*3;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
		#pragma omp section
					{
						for(n=N/4*3;n<N;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
				}
	MatrixXd fraz(4*Nin1,2);
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
		}
	    //Put the weights into a 3d matrix, where each row in the third
	    //dimension is a new particle. Adding all allements on each row should
	    //equal 1 in 4 adjacent indices

	//Compute the charge density
	VectorXd rho(NG*NG);
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<NG*NG/4;n++)
				{
						rho(n)=0;
					}
			}
#pragma omp section
			{
				for(n=NG*NG/4;n<NG*NG/2;n++)
				{
						rho(n)=0;
					}
			}
#pragma omp section
			{
				for(n=NG*NG/2;n<NG*NG/4*3;n++)
				{
						rho(n)=0;
					}
			}
#pragma omp section
			{
				for(n=NG*NG/4*3;n<NG*NG;n++)
				{
						rho(n)=0;
					}
			}
		}

#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
		}

	rho=rho+esphere1*rho_sphere1+esphere2*rho_sphere2;


	VectorXd Phiseg1(NG*NG/4);
	VectorXd Phiseg2(NG*NG/4);
	VectorXd Phiseg3(NG*NG/4);
	VectorXd Phiseg4(NG*NG/4);

	double PhiC1;
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				Phiseg1=-invPoi.block(0,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg2=-invPoi.block(NG*NG/4,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg3=-invPoi.block(NG*NG/2,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg4=-invPoi.block(3*NG*NG/4,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
		}


	VectorXd Phiin1(lengthin1);
	VectorXd Phitemp1(lengthin1);
	double PhiC2;
	VectorXd Phiin2(lengthin2);
	VectorXd Phitemp2(lengthin2);

	VectorXd PhiC1C2(lengthin3);

	VectorXd drho3=Cap3.transpose()*PhiC1C2/dx/dx;
	VectorXd drho1(NG*NG);
	VectorXd drho2(NG*NG);
int m=0;

	MatrixXd rho1(NG,NG);
	for(n=0;n<NG;n++)
	{
	    for(m=0;m<NG;m++)
	    {
	    	rho1(n,m)=rho(n*NG+m);
	    }
	}

	if(it>19999)
	{
		rho1allit=(it-20000)*rho1allit/(it-19999)+rho1/(it-19999);
	}



	MatrixXd Phi1(NG,NG);
	VectorXd Phi(NG*NG);

	for(n=0;n<NG;n++)
	{
	    for(m=0;m<NG;m++)
	    {
			if(n*NG+m<NG*NG/4)
			{
		    	Phi(n*NG+m)=Phiseg1(n*NG+m);
				Phi1(n,m)=Phiseg1(n*NG+m);
			}
			else if(n*NG+m<NG*NG/2)
			{
		    	Phi(n*NG+m)=Phiseg2(n*NG+m-NG*NG/4);
		    	Phi1(n,m)=Phiseg2(n*NG+m-NG*NG/4);
			}
			else if(n*NG+m<3*NG*NG/4)
			{
		    	Phi(n*NG+m)=Phiseg3(n*NG+m-NG*NG/2);
		    	Phi1(n,m)=Phiseg3(n*NG+m-NG*NG/2);
			}
			else
			{
		    	Phi(n*NG+m)=Phiseg4(n*NG+m-3*NG*NG/4);
		    	Phi1(n,m)=Phiseg4(n*NG+m-3*NG*NG/4);
			}
	    }
	}


if(it>19999)
{
	Phi1allit=(it-20000)*Phi1allit/(it-19999)+Phi1/(it-19999);
}





		for(n=0;n<N;n++)
				{
					Pot(n)=0;
					if(index2(g(n,0),g(n,1),NG)<NG*NG/4)
					{
						Pot(n)+=Phiseg1(index2(g(n,0),g(n,1),NG))*fraz(n,0);
					}
					else if(index2(g(n,0),g(n,1),NG)<NG*NG/2)
					{
						Pot(n)+=Phiseg2(index2(g(n,0),g(n,1),NG)-NG*NG/4)*fraz(n,0);
					}
					else if(index2(g(n,0),g(n,1),NG)<NG*NG*3/4)
					{
						Pot(n)+=Phiseg3(index2(g(n,0),g(n,1),NG)-NG*NG/2)*fraz(n,0);
					}
					else
					{
						Pot(n)+=Phiseg4(index2(g(n,0),g(n,1),NG)-NG*NG*3/4)*fraz(n,0);
					}

					if(index2(g(n,0),g(N+n,1),NG)<NG*NG/4)
					{
						Pot(n)+=Phiseg1(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0);
					}
					else if(index2(g(n,0),g(N+n,1),NG)<NG*NG/2)
					{
						Pot(n)+=Phiseg2(index2(g(n,0),g(N+n,1),NG)-NG*NG/4)*fraz(n+N,0);
					}
					else if(index2(g(n,0),g(N+n,1),NG)<NG*NG*3/4)
					{
						Pot(n)+=Phiseg3(index2(g(n,0),g(N+n,1),NG)-NG*NG/2)*fraz(n+N,0);
					}
					else
					{
						Pot(n)+=Phiseg4(index2(g(n,0),g(N+n,1),NG)-NG*NG*3/4)*fraz(n+N,0);
					}

					if(index2(g(N+n,0),g(n,1),NG)<NG*NG/4)
					{
						Pot(n)+=Phiseg1(index2(g(N+n,0),g(n,1),NG))*fraz(n,1);
					}
					else if(index2(g(N+n,0),g(n,1),NG)<NG*NG/2)
					{
						Pot(n)+=Phiseg2(index2(g(N+n,0),g(n,1),NG)-NG*NG/4)*fraz(n,1);
					}
					else if(index2(g(N+n,0),g(n,1),NG)<NG*NG*3/4)
					{
						Pot(n)+=Phiseg3(index2(g(N+n,0),g(n,1),NG)-NG*NG/2)*fraz(n,1);
					}
					else
					{
						Pot(n)+=Phiseg4(index2(g(N+n,0),g(n,1),NG)-NG*NG*3/4)*fraz(n,1);
					}

					if(index2(g(N+n,0),g(N+n,1),NG)<NG*NG/4)
					{
						Pot(n)+=Phiseg1(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1);
					}
					else if(index2(g(N+n,0),g(N+n,1),NG)<NG*NG/2)
					{
						Pot(n)+=Phiseg2(index2(g(N+n,0),g(N+n,1),NG)-NG*NG/4)*fraz(n+N,1);
					}
					else if(index2(g(N+n,0),g(N+n,1),NG)<NG*NG*3/4)
					{
						Pot(n)+=Phiseg3(index2(g(N+n,0),g(N+n,1),NG)-NG*NG/2)*fraz(n+N,1);
					}
					else
					{
						Pot(n)+=Phiseg4(index2(g(N+n,0),g(N+n,1),NG)-NG*NG*3/4)*fraz(n+N,1);
					}
				}

				double tempvel;
				double expvel;


//			    //Update Velocity
//
//
//			#pragma omp parallel sections num_threads(4)
//							{
//					#pragma omp section
//								{
//									for(n=0;n<N/4;n++)
//									{
//										g(n,0)=floor(xp(n)/dx);
//										g(n,1)=floor(yp(n)/dx);
//										g(N+n,0)=(floor(xp(n)/dx)+1);
//										g(N+n,1)=(floor(yp(n)/dx)+1);
//									}
//								}
//					#pragma omp section
//								{
//									for(n=N/4;n<N/2;n++)
//									{
//										g(n,0)=floor(xp(n)/dx);
//										g(n,1)=floor(yp(n)/dx);
//										g(N+n,0)=(floor(xp(n)/dx)+1);
//										g(N+n,1)=(floor(yp(n)/dx)+1);
//									}
//								}
//					#pragma omp section
//								{
//									for(n=N/2;n<N/4*3;n++)
//									{
//										g(n,0)=floor(xp(n)/dx);
//										g(n,1)=floor(yp(n)/dx);
//										g(N+n,0)=(floor(xp(n)/dx)+1);
//										g(N+n,1)=(floor(yp(n)/dx)+1);
//									}
//								}
//					#pragma omp section
//								{
//									for(n=N/4*3;n<N;n++)
//									{
//										g(n,0)=floor(xp(n)/dx);
//										g(n,1)=floor(yp(n)/dx);
//										g(N+n,0)=(floor(xp(n)/dx)+1);
//										g(N+n,1)=(floor(yp(n)/dx)+1);
//									}
//								}
//							}
//
//			#pragma omp parallel sections num_threads(4)
//					{
//			#pragma omp section
//						{
//							for(n=0;n<N/4;n++)
//							{
//								fraz(n,0)=(1-(xp(n)/dx-g(n,0)))*(1-(yp(n)/dx-g(n,1)));
//								fraz(n,1)=(xp(n)/dx-g(n,0))*(1-(yp(n)/dx-g(n,1)));
//								fraz(N+n,0)=(1-(xp(n)/dx-g(n,0)))*(yp(n)/dx-g(n,1));
//								fraz(N+n,1)=(xp(n)/dx-g(n,0))*(yp(n)/dx-g(n,1));
//							}
//						}
//			#pragma omp section
//						{
//							for(n=N/4;n<N/2;n++)
//							{
//								fraz(n,0)=(1-(xp(n)/dx-g(n,0)))*(1-(yp(n)/dx-g(n,1)));
//								fraz(n,1)=(xp(n)/dx-g(n,0))*(1-(yp(n)/dx-g(n,1)));
//								fraz(N+n,0)=(1-(xp(n)/dx-g(n,0)))*(yp(n)/dx-g(n,1));
//								fraz(N+n,1)=(xp(n)/dx-g(n,0))*(yp(n)/dx-g(n,1));
//							}
//						}
//			#pragma omp section
//						{
//							for(n=N/2;n<N/4*3;n++)
//							{
//								fraz(n,0)=(1-(xp(n)/dx-g(n,0)))*(1-(yp(n)/dx-g(n,1)));
//								fraz(n,1)=(xp(n)/dx-g(n,0))*(1-(yp(n)/dx-g(n,1)));
//								fraz(N+n,0)=(1-(xp(n)/dx-g(n,0)))*(yp(n)/dx-g(n,1));
//								fraz(N+n,1)=(xp(n)/dx-g(n,0))*(yp(n)/dx-g(n,1));
//							}
//						}
//			#pragma omp section
//						{
//							for(n=N/4*3;n<N;n++)
//							{
//								fraz(n,0)=(1-(xp(n)/dx-g(n,0)))*(1-(yp(n)/dx-g(n,1)));
//								fraz(n,1)=(xp(n)/dx-g(n,0))*(1-(yp(n)/dx-g(n,1)));
//								fraz(N+n,0)=(1-(xp(n)/dx-g(n,0)))*(yp(n)/dx-g(n,1));
//								fraz(N+n,1)=(xp(n)/dx-g(n,0))*(yp(n)/dx-g(n,1));
//							}
//						}
//					}




					//Electric Field in x direction
					int in1;
					int in2;
#pragma omp parallel sections num_threads(4)
							{
					#pragma omp section
								{
									for(n=0;n<NG*NG/4;n++)
									{
										        in1=n+NG;//Index of right node
										        in2=n-NG;//Index of left node
										        if(in1>NG*NG-1)
										        {
										            Ex(n)=(Phi(in1-NG*NG)-Phi(in2))/(2*dx);//Electric field in x, zero out edge
										        }
										        else if (in2<0)
												{
										            Ex(n)=(Phi(in1)-Phi(in2+NG*NG))/(2*dx);//Electric field in x, zero out edge
												}
										        else
										        {
										        	Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in x
										        }
										}
								}
					#pragma omp section
								{
									for(n=NG*NG/4;n<NG*NG/2;n++)
									{
										        in1=n+NG;//Index of right node
										        in2=n-NG;//Index of left node
										        if(in1>NG*NG-1)
										        {
										            Ex(n)=(Phi(in1-NG*NG)-Phi(in2))/(2*dx);//Electric field in x, zero out edge
										        }
										        else if (in2<0)
												{
										            Ex(n)=(Phi(in1)-Phi(in2+NG*NG))/(2*dx);//Electric field in x, zero out edge
												}
										        else
										        {
										        	Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in x
										        }
										}
								}
					#pragma omp section
								{
									for(n=NG*NG/2;n<NG*NG/4*3;n++)
									{
										        in1=n+NG;//Index of right node
										        in2=n-NG;//Index of left node
										        if(in1>NG*NG-1)
										        {
										            Ex(n)=(Phi(in1-NG*NG)-Phi(in2))/(2*dx);//Electric field in x, zero out edge
										        }
										        else if (in2<0)
												{
										            Ex(n)=(Phi(in1)-Phi(in2+NG*NG))/(2*dx);//Electric field in x, zero out edge
												}
										        else
										        {
										        	Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in x
										        }
										}
								}
					#pragma omp section
								{
									for(n=NG*NG/4*3;n<NG*NG;n++)
									{
										        in1=n+NG;//Index of right node
										        in2=n-NG;//Index of left node
										        if(in1>NG*NG-1)
										        {
										            Ex(n)=(Phi(in1-NG*NG)-Phi(in2))/(2*dx);//Electric field in x, zero out edge
										        }
										        else if (in2<0)
												{
										            Ex(n)=(Phi(in1)-Phi(in2+NG*NG))/(2*dx);//Electric field in x, zero out edge
												}
										        else
										        {
										        	Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in x
										        }
										}
								}
							}

						for(n=0;n<NG;n++)
						{
						    for(m=0;m<NG;m++)
						    {
						    	Ex1(n,m)=Ex(n*NG+m);
						    }
						}

					//Electric Field in Y Direction
				#pragma omp parallel sections num_threads(4)
						{
				#pragma omp section
							{
								for(n=0;n<NG*NG/4;n++)
								{
									        in1=n+1;//Index of above node
									        in2=n-1;//index of lower node
									        if (in1%NG==0)
									        {
									            Ey(n)=(Phi(in1-NG)-Phi(in2))/(2*dx);//Electric field in y, zere out edge
									        }
									        else if ((in2+NG)%NG==NG-1)
									        {
									            Ey(n)=(Phi(in1)-Phi(in2+NG))/(2*dx);//Electric field in y, zero out edge
									        }
									        else
									        {
									        	Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in y
									        }
									}
							}
				#pragma omp section
							{
								for(n=NG*NG/4;n<NG*NG/2;n++)
								{
									        in1=n+1;//Index of above node
									        in2=n-1;//index of lower node
									        if (in1%NG==0)
									        {
									            Ey(n)=(Phi(in1-NG)-Phi(in2))/(2*dx);//Electric field in y, zere out edge
									        }
									        else if ((in2+NG)%NG==NG-1)
									        {
									            Ey(n)=(Phi(in1)-Phi(in2+NG))/(2*dx);//Electric field in y, zero out edge
									        }
									        else
									        {
									        	Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in y
									        }
									}
							}
				#pragma omp section
							{
								for(n=NG*NG/2;n<NG*NG/4*3;n++)
								{
									        in1=n+1;//Index of above node
									        in2=n-1;//index of lower node
									        if (in1%NG==0)
									        {
									            Ey(n)=(Phi(in1-NG)-Phi(in2))/(2*dx);//Electric field in y, zere out edge
									        }
									        else if ((in2+NG)%NG==NG-1)
									        {
									            Ey(n)=(Phi(in1)-Phi(in2+NG))/(2*dx);//Electric field in y, zero out edge
									        }
									        else
									        {
									        	Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in y
									        }
									}
							}
				#pragma omp section
							{
								for(n=NG*NG/4*3;n<NG*NG;n++)
								{
									        in1=n+1;//Index of above node
									        in2=n-1;//index of lower node
									        if (in1%NG==0)
									        {
									            Ey(n)=(Phi(in1-NG)-Phi(in2))/(2*dx);//Electric field in y, zere out edge
									        }
									        else if ((in2+NG)%NG==NG-1)
									        {
									            Ey(n)=(Phi(in1)-Phi(in2+NG))/(2*dx);//Electric field in y, zero out edge
									        }
									        else
									        {
									        	Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in y
									        }
									}
							}
						}

							for(n=0;n<NG;n++)
							{
							    for(m=0;m<NG;m++)
							    {
							    	Ey1(n,m)=Ey(n*NG+m);
							    }
							}


				    //Put the weights into a 3d matrix, where each row in the third
				    //dimension is a new particle. Adding all allements on each row should
				    //equal 1 in 4 adjacent indices


						    //Update Velocity
					#pragma omp parallel sections num_threads(4)
							{
					#pragma omp section
								{
									for(n=0;n<N/4;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
					#pragma omp section
								{
									for(n=N/4;n<N/2;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
					#pragma omp section
								{
									for(n=N/2;n<N/4*3;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
					#pragma omp section
								{
									for(n=N/4*3;n<N;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
							}


	cout << "\nit="<<it;

	double Phiin1a;
	double Phiin2a;
	for (n=0;n<NG*NG;n++)
	{
		if(n%NG>2&&n%NG<NG-3)
		{
			in1=n-1;
			in2=n+1;
		if(in1<NG*NG/4)
		{
			Phiin1a=Phiseg1(in1);
		}
		else if(in1<NG*NG/2)
		{
			Phiin1a=Phiseg2(in1-NG*NG/4);
		}
		else if(in1<NG*NG*3/4)
		{
			Phiin1a=Phiseg3(in1-NG*NG/2);
		}
		else
		{
			Phiin1a=Phiseg4(in1-NG*NG*3/4);
		}
		if(in2<NG*NG/4)
		{
			Phiin2a=Phiseg1(in2);
		}
		else if(in2<NG*NG/2)
		{
			Phiin2a=Phiseg2(in2-NG*NG/4);
		}
		else if(in2<NG*NG*3/4)
		{
			Phiin2a=Phiseg3(in2-NG*NG/2);
		}
		else
		{
			Phiin2a=Phiseg4(in2-NG*NG*3/4);
		}
		ElectricForceSphere1X(it)+=(Phiin2a-Phiin1a)*dx*(esphere1*rho_sphere1(n)+drho1(n))/2;
		}
	}


	for (n=0;n<NG*NG;n++)
	{
		if(n>=2*NG&&n<NG*NG-2*NG)
		{
			in1=n-NG;
			in2=n+NG;
		if(in1<NG*NG/4)
		{
			Phiin1a=Phiseg1(in1);
		}
		else if(in1<NG*NG/2)
		{
			Phiin1a=Phiseg2(in1-NG*NG/4);
		}
		else if(in1<NG*NG*3/4)
		{
			Phiin1a=Phiseg3(in1-NG*NG/2);
		}
		else
		{
			Phiin1a=Phiseg4(in1-NG*NG*3/4);
		}
		if(in2<NG*NG/4)
		{
			Phiin2a=Phiseg1(in2);
		}
		else if(in2<NG*NG/2)
		{
			Phiin2a=Phiseg2(in2-NG*NG/4);
		}
		else if(in2<NG*NG*3/4)
		{
			Phiin2a=Phiseg3(in2-NG*NG/2);
		}
		else
		{
			Phiin2a=Phiseg4(in2-NG*NG*3/4);
		}
		ElectricForceSphere1Y(it)+=(Phiin2a-Phiin1a)*dx*(esphere1*rho_sphere1(n)+drho1(n))/2;
		}
	}

	for (n=0;n<NG*NG;n++)
	{
		if(n%NG>2&&n%NG<NG-3)
		{
			in1=n-1;
			in2=n+1;
		if(in1<NG*NG/4)
		{
			Phiin1a=Phiseg1(in1);
		}
		else if(in1<NG*NG/2)
		{
			Phiin1a=Phiseg2(in1-NG*NG/4);
		}
		else if(in1<NG*NG*3/4)
		{
			Phiin1a=Phiseg3(in1-NG*NG/2);
		}
		else
		{
			Phiin1a=Phiseg4(in1-NG*NG*3/4);
		}
		if(in2<NG*NG/4)
		{
			Phiin2a=Phiseg1(in2);
		}
		else if(in2<NG*NG/2)
		{
			Phiin2a=Phiseg2(in2-NG*NG/4);
		}
		else if(in2<NG*NG*3/4)
		{
			Phiin2a=Phiseg3(in2-NG*NG/2);
		}
		else
		{
			Phiin2a=Phiseg4(in2-NG*NG*3/4);
		}
		ElectricForceSphere2X(it)+=(Phiin2a-Phiin1a)*dx*(esphere1*rho_sphere2(n)+drho2(n))/2;
		}
	}


	for (n=0;n<NG*NG;n++)
	{
		if(n>=2*NG&&n<NG*NG-2*NG)
		{
			in1=n-NG;
			in2=n+NG;
		if(in1<NG*NG/4)
		{
			Phiin1a=Phiseg1(in1);
		}
		else if(in1<NG*NG/2)
		{
			Phiin1a=Phiseg2(in1-NG*NG/4);
		}
		else if(in1<NG*NG*3/4)
		{
			Phiin1a=Phiseg3(in1-NG*NG/2);
		}
		else
		{
			Phiin1a=Phiseg4(in1-NG*NG*3/4);
		}
		if(in2<NG*NG/4)
		{
			Phiin2a=Phiseg1(in2);
		}
		else if(in2<NG*NG/2)
		{
			Phiin2a=Phiseg2(in2-NG*NG/4);
		}
		else if(in2<NG*NG*3/4)
		{
			Phiin2a=Phiseg3(in2-NG*NG/2);
		}
		else
		{
			Phiin2a=Phiseg4(in2-NG*NG*3/4);
		}
		ElectricForceSphere2Y(it)+=(Phiin2a-Phiin1a)*dx*(esphere1*rho_sphere2(n)+drho2(n))/2;
		}
	}

//	cout << "\nForce 1 X:" << ElectricForceSphere1X(it);
//	cout << "\nForce 2 X:" << ElectricForceSphere2X(it);
//	cout << "\nForce 1 Y:" << ElectricForceSphere1Y(it);
//	cout << "\nForce 2 Y:" << ElectricForceSphere2Y(it);



	double Temp=0;
	for(n=0;n<N;n++)
	{
		Temp+=massmat(n)/me*(vxp(n)*vxp(n)+vyp(n)*vyp(n))/N/2/k;
	}

	totUEn(it)=0;
	for(n=0;n<NG*NG;n++)
	{
		totUEn(it)+=0.5*dx*dx*e0*(Ex(n)*Ex(n)+Ey(n)*Ey(n));
	}

	totKEn(it)=0;

	for(n=0;n<N;n++)
	{
		totKEn(it)+=.5*massmat(n)*(vxp(n)*vxp(n)+vyp(n)*vyp(n));
	}

	totEn(it)=totUEn(it)+totKEn(it);

	cout << "\nTotal Energy: " <<setprecision(20)<< totEn(it);

	totMomX(it)=0;
	totMomY(it)=0;

	for(n=0;n<N;n++)
	{
		totMomX(it)+=massmat(n)*vxp(n);
		totMomY(it)+=massmat(n)*vyp(n);
	}

	cout << "\nTotal X Momentum: " <<setprecision(20)<< totMomX(it);
	cout << "\nTotal Y Momentum: " << setprecision(20)<<totMomY(it);

	cout << "\nTemp=" <<setprecision(20)<< Temp/11700;


	Temp1(it)=Temp;

	LD =sqrt(k*Temp/ndensity);
			if(.5*LD<dx)//Compare to delta x, if too small, display error
			{
			    cout<<"\nError, Mesh not fine enough or not hot enough\nLD="<<LD<<"\ndx="<<dx;
			}
			double PP=sqrt(pi*me*L*L/(N*e*e));//Compute plasma period
			if(.25*PP<DT)//If plasma period is small display error
			{
				cout<<"/nError, large time step\nPP="<<PP<<"\nDT="<<DT;
			}

			//CFL Condition
			double C=DT/dx;
			if(C>1)//If time step is larger than dx
			{
			    cout <<"\nMay be unstable, C>1.\nC="<<C;
			}


//			cout << "\nxp="<<xp(1)/L;
//			cout << "\nyp="<<yp(1)/L;
//			cout << "\nfraz1="<<fraz(1,0);
//			cout << "\nfraz2="<<fraz(1,1);
//			cout << "\nfraz3="<<fraz(1+N,0);
//			cout << "\nfraz4="<<fraz(1+N,1);
//			cout << "\ng"<<g(1,0)<<","<<g(1,1)<<","<<g(1+N,0)<<","<<g(N+1,1);
//
//			cout <<"\nEx="<< Ex(index2(g(1,0),g(1,1),NG))<<",";
//			cout << Ex(index2(g(1,0),g(N+1,1),NG))<<",";
//			cout << Ex(index2(g(N+1,0),g(1,1),NG))<<",";
//			cout << Ex(index2(g(N+1,0),g(N+1,1),NG));
//
//
//			cout <<"\nEy="<< Ey(index2(g(1,0),g(1,1),NG))<<",";
//			cout << Ey(index2(g(1,0),g(N+1,1),NG))<<",";
//			cout << Ey(index2(g(N+1,0),g(1,1),NG))<<",";
//			cout << Ey(index2(g(N+1,0),g(N+1,1),NG));
//
//			cout << "\nEx1="<< Ex1;
//			cout << "\nEy1="<< Ey1;
//			cout << "\nrho1="<< rho1;
//			cout << "\nPhi1="<< Phi1;








	if(it%10000==9999)
	{

		for(n=0;n<10000;n++)
		{
			myfile1 <<Volt1(it-9999+n)<<"\n";

			myfile2 << Volt2(it-9999+n)<<"\n";

			myfile3 << Charge1(it-9999+n)<<"\n";

			myfile4 << Charge2(it-9999+n)<<"\n";

			myfile7 << setprecision(20)<<Temp1(it-9999+n)<<"\n";

			myfile9 << setprecision(20)<<totUEn(it-9999+n)<<"\n";

			myfile10 << setprecision(20)<<totKEn(it-9999+n)<<"\n";

			myfile11 << setprecision(20)<<totEn(it-9999+n)<<"\n";

			myfile23 << TotalMomentumAbsorbed1X(it-9999+n)<<"\n";

			myfile24 << TotalMomentumAbsorbed1Y(it-9999+n)<<"\n";

			myfile25 << TotalMomentumAbsorbed2X(it-9999+n)<<"\n";

			myfile26 << TotalMomentumAbsorbed2Y(it-9999+n)<<"\n";

			myfile27 << ElectricForceSphere1X(it-9999+n)<<"\n";

			myfile28 << ElectricForceSphere1Y(it-9999+n)<<"\n";

			myfile29 << ElectricForceSphere2X(it-9999+n)<<"\n";

			myfile30 << ElectricForceSphere2Y(it-9999+n)<<"\n";

			myfile31 << setprecision(20)<<totMomX(it-9999+n)<<"\n";

			myfile32 << setprecision(20)<<totMomY(it-9999+n)<<"\n";
		}

		myfile5.open(path+to_string(number)+"rho1allit.txt");
		myfile5 << rho1allit;
		myfile5.close();

		myfile6.open(path+to_string(number)+"Phi1allit.txt");
		myfile6 << Phi1allit;
		myfile6.close();

		myfile8.open(path+to_string(number)+"velocitysquared.txt");
		myfile8 << vxp.array().square()+vyp.array().square();
		myfile8.close();

		myfile12.open(path+to_string(number)+"Phi1.txt");
		myfile12 << Phi1;
		myfile12.close();

		myfile13.open(path+to_string(number)+"rho1.txt");
		myfile13 << rho1;
		myfile13.close();

		myfile14.open(path+to_string(number)+"emat.txt");
		myfile14 << emat;
		myfile14.close();

		myfile15.open(path+to_string(number)+"xp.txt");
		myfile15 << xp;
		myfile15.close();

		myfile16.open(path+to_string(number)+"yp.txt");
		myfile16 << yp;
		myfile16.close();

		myfile17.open(path+to_string(number)+"vxp.txt");
		myfile17 << vxp;
		myfile17.close();

		myfile18.open(path+to_string(number)+"vyp.txt");
		myfile18 << vyp;
		myfile18.close();



	}

	}


	myfile1.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
	myfile7.close();
	myfile9.close();
	myfile10.close();
	myfile11.close();
	myfile23.close();
	myfile24.close();
	myfile25.close();
	myfile26.close();
	myfile27.close();
	myfile28.close();
	myfile29.close();
	myfile30.close();

	}


	cout<< "\nTime="<<(clock ()-starttime)/CLOCKS_PER_SEC;



return 0;
}

//
//	 ForceXSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1)'*Ex;
//	 ForceYSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1)'*Ey;
//	 ForceXSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2)'*Ex;
//	 ForceYSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2)'*Ey;
//
//	 Force1=ForceYSphere1(it)
//	 Force2=ForceYSphere2(it)
//
//
//
//
//}

