#include<iostream>
#include<stdio.h>
#include<math.h>
#include <cmath>
#include<string.h>
#include<iomanip>
#include<fstream>
#include<string>
#include <stdlib.h>
#include <time.h>
#include <algorithm>



#define NI 22 //nodes in x direction
#define NJ 42 //nodes in y direction

using namespace std;


//float CELL_ZI[NI+1][NJ+1],CELL_ZETA[NI+1][NJ+1],

/////////////////////
/* geometric data */
////////////////////

float L,H;
int i,j;
double X[NJ][NI], Y[NJ][NI];
double  DELX[NJ][NI], DELY[NJ][NI];
double  Y_FACE_CENTER_X[NJ][NI], X_FACE_CENTER_Y[NJ][NI];
double  X_NODE[NJ+1][NI+1], Y_NODE[NJ+1][NI+1], DXF[4][NJ+1][NI+1];//DXF node distance
double  S[4][NJ][NI]; //surface area
double  VOL[NJ][NI];  //volume 
double delx,dely, ST_FACTOR_X, ST_FACTOR_Y; 
double Rx, Ry, rx, ry, Alphax, Alphay, delxs, delys; //strech mesh parameters
//delxs smallest delx
//delys smallest dely

//////////////////////////
/* heat conduction data */
//////////////////////////

double AW[NJ+1][NI+1], AS[NJ+1][NI+1], AE[NJ+1][NI+1], AN[NJ+1][NI+1], AP[NJ+1][NI+1];
double SP[NJ+1][NI+1], k[NI+1][NJ+1], kw[NI+1][NJ+1], ks[NI+1][NJ+1], ke[NI+1][NJ+1], kn[NI+1][NJ+1], kx[NI+1][NJ+1], ky[NI+1][NJ+1] ; 
double T_OLD[NJ+1][NI+1],T[NJ+1][NI+1], qx[NJ+1][NI+1],qy[NJ+1][NI+1],Q_VOL_GEN ; 
 //S = -1.5
///////////////////////////
/* convergence and error */
///////////////////////////

double RESIDUE[NJ+1][NI+1];
double MAX_ERROR, ABS_ERROR, RMS_ERROR, RMSRESIDUE;
int ITER = 0;

void SET_GEOMETRY_UNIFORM();
void SET_GEOMETRY_NONUNIFORM();
void APPLY_IC();
void APPLY_BC();
void UPDATE();
void CALC_CONDUCT();
double CALC_ABS_ERROR();
double CALC_RMS_ERROR();

//void flux();
//void temporal();
int FILE_WRITE1();
int FILE_WRITE2();

int main ()
{
	L = 1; H = 0.5;
	
	ST_FACTOR_Y = 1.1;
	ST_FACTOR_X = 1.08;
	
	MAX_ERROR = 1e-06;
	Q_VOL_GEN = -1.5;
	
	SET_GEOMETRY_UNIFORM();
//  SET_GEOMETRY_NONUNIFORM();
    APPLY_BC();
    APPLY_IC();
    
    UPDATE();
    ofstream file3("ITER_VS_RESIDUE_42_22.dat");
ITER:

   ++ITER;
       
    CALC_CONDUCT();

    APPLY_BC();
    
    RMSRESIDUE = CALC_RMS_ERROR();
    
    file3<<ITER<<"\t"<<RMSRESIDUE<<"\n";
    cout<<ITER<<" "<<RMSRESIDUE<<" "<<endl;

	if(RMSRESIDUE > MAX_ERROR){
		UPDATE();
		goto ITER;
	}
    
    file3.close();
    FILE_WRITE1();
     

}


void SET_GEOMETRY_UNIFORM()
{
	                     // in compuational domain generates 1x1 squre domain.
	delx = (L/(NI-2));    // X is horizontal direction in physical domain
	dely = (H/(NJ-2));  // Y is vertical direction in physical domain
	//delv = delx * dely ;
//	cout<<delv<<endl ;
//	cout<<delx<<" "<<dely<<endl;

for(j=1;j<=NJ-1;j++)
	{
		for(i=1;i<=NI-1;i++)
		{
            X[j][i]= (i-1)*delx;
			Y[j][i]= (j-1)*dely;	
	    //  cout<<X_FACE_CENTER_Y[j][i]<<endl;
	    } 
	}	
			
for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
		
	    	X_NODE[j][i] = (X[j][i]+X[j][i-1])/2;
			Y_NODE[j][i] = (Y[j][i]+Y[j-1][i])/2;       
		} 
	}	

for(j=2;j<=NJ-1;j++) // left and right boundary
	  {
	  	X_NODE[j][1] = 0;
	  	Y_NODE[j][1] = Y_NODE[j][2];
	  	X_NODE[j][NI] = L;
	  	Y_NODE[j][NI] =Y_NODE[j][NI-1];
	  }
for(i=2;i<=NI-1;i++) //top and bottom boundary
	  {
	  	X_NODE[1][i] = X_NODE[2][i];
	  	Y_NODE[1][i] = 0;
	  	X_NODE[NJ][i] = X_NODE[NJ-1][i];
	  	Y_NODE[NJ][i] = H;
	  }
	  
	  X_NODE[1][1] = 0; Y_NODE[1][1] = 0;  //four corner of mesh
	  X_NODE[1][NI] = L; Y_NODE[1][NI]=0;
	  X_NODE[NJ][1] = 0; Y_NODE[NJ][1] = H; 
	  X_NODE[NJ][NI] = L; Y_NODE[NJ][NI] = H;

//face center
for(j=2;j<=NJ-1;j++)
	{
		for(i=1;i<=NI-1;i++)
		{
			Y_FACE_CENTER_X[j][i]= (i-1)*delx;
		//	cout<<Y_FACE_CENTER_X[j][i]<<endl;
	    } 
	}	
	
for(j=1;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
            X_FACE_CENTER_Y[j][i]= (j-1)*dely;	
	    //  cout<<X_FACE_CENTER_Y[j][i]<<endl;
	    } 
	}	
	  
//NODE DISTANCE

for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
		DXF[0][j][i] = (X_NODE[j][i]-X_NODE[j][i-1]); //WEST
		DXF[1][j][i] = (Y_NODE[j][i]-Y_NODE[j-1][i]); //SOUTH
		DXF[2][j][i] = (X_NODE[j][i+1]-X_NODE[j][i]); //EAST
		DXF[3][j][i] = (Y_NODE[j+1][i]-Y_NODE[j][i]); //NORTH
	  //  cout<<DXF[0][j][i]<<" "<<DXF[1][j][i]<<" "<<DXF[2][j][i]<<" "<<DXF[3][j][i]<<endl;
	  //  cout<<"\n"<<endl;
		}    
    }
		
//SURFACE AREA

for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
			S[0][j][i] = dely; //WEST FACE
			S[1][j][i] = delx;	//SOUTH FACE 
	        S[2][j][i] = dely; //EAST FACE
	        S[3][j][i] = delx;    //NORTH FACE
	        
	        VOL[j][i] = S[2][j][i] * S[3][j][i] ;
	    //    cout<<S[2][j][i]<<" "<<S[3][j][i]<<endl;
	        
        }    
    }
				
	  
}  

void SET_GEOMETRY_NONUNIFORM()
{
	                     // in compuational domain generates 1x1 squre domain.
//	delx = (L/(NI-2));    // X is horizontal direction in physical domain
	dely = (H/(NJ-2));  // Y is vertical direction in physical domain
  
  ry = ST_FACTOR_Y;
  rx = ST_FACTOR_X;
	
  // GP series followed for streching the mesh
//	rx = pow(Rx,(1/NI-1));
    Rx = pow(rx,(NI-3));
	if(Rx > 1){
		Alphax = Rx; 
	}
	else{
		Alphax = 1-(pow(rx,-(NI-2)))+(pow(rx,-1));
	}
	delxs = L*(rx-1)/((Alphax*rx)-1);
  
 /*    
//   ry = pow(Ry,(1/NJ-1));
    Ry = pow(ry,(NJ-3));
    
	if(Ry > 1){
		Alphay = Ry;
	}
	else{
		Alphay = 1-(pow(ry,-(NJ-2)))+(pow(ry,-1));
	}
	delys = L*(ry-1)/((Alphay*ry)-1);	
	
*/	
		
for(j=1;j<=NJ-1;j++)
	{
	if(j == 1){ 	
	  DELX[j][1] = 0;
//	  DELY[j][1] = 0;
        }  

    else if(j == 2){
      DELX[j][1] = 0;
//	  DELY[j][1] = delys; 	
	    } 
    
	else {
	  DELX[j][1] = 0;
//	  DELY[j][1] = DELY[j-1][1] + (pow(ry,(j-2))*delys);	
     	}   
    }

/*
for(i=1;i<=NI-1;i++)
	{
	if(i == 1){ 	
	  DELX[1][i] = 0;
//	  DELY[1][i] = 0;
        }  
   
    else {
	  DELX[1][i] = DELX[1][i-1] + (pow(rx,(i-2))*delxs);
   // DELY[1][i] = 0;	
     	}   
     cout<<DELX[1][i]<<endl;
    }
*/


for(i=NI-1;i>=1;i--)
	{
	if(i == NI-1){ 	
	  DELX[1][i] = L;
//	  DELY[1][i] = 0;
        }  
     
	else {
	  DELX[1][i] = DELX[1][i+1] - (pow(rx,(NI-2-i))*delxs);
   // DELY[1][i] = 0;	
     	}   
    // cout<<DELX[1][i]<<endl;
     }

    
for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
			DELX[j][i] = DELX [1][i];
		//	DELY[j][i] = DELY[j][i-1];	
        }    
    }
	
for(j=1;j<=NJ-1;j++)
	{
		for(i=1;i<=NI-1;i++)
		{
		  X[j][i]= DELX[j][i];
		  Y[j][i]= (j-1)*dely;//DELY[j][i];	
	    }    
    }

for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{	
		   //Node Co-ordinate for TEMPERATURE CV
	      X_NODE[j][i] = (X[j][i]+X[j][i-1])/2;
		  Y_NODE[j][i] = (Y[j][i]+Y[j-1][i])/2;
	    	
	    } 
	 }			
				 
	  
for(j=2;j<=NJ-1;j++) // left and right boundary
	  {
	  	X_NODE[j][1] = 0;
	  	Y_NODE[j][1] = Y_NODE[j][2];
	  	X_NODE[j][NI] = L;
	  	Y_NODE[j][NI] =Y_NODE[j][NI-1];
	  }
for(i=2;i<=NI-1;i++) //top and bottom boundary
	  {
	  	X_NODE[1][i] = X_NODE[2][i];
	  	Y_NODE[1][i] = 0;
	  	X_NODE[NJ][i] = X_NODE[NJ-1][i];
	  	Y_NODE[NJ][i] = H;
	  }
	  
	  X_NODE[1][1] = 0; Y_NODE[1][1] = 0;  //four corner of mesh
	  X_NODE[1][NI] = L; Y_NODE[1][NI]=0;
	  X_NODE[NJ][1] = 0; Y_NODE[NJ][1] = H; 
	  X_NODE[NJ][NI] = L; Y_NODE[NJ][NI] = H;		

//face center
for(j=2;j<=NJ-1;j++)
	{
		for(i=1;i<=NI-1;i++)
		{
			Y_FACE_CENTER_X[j][i]= X[j][i];
		//	cout<<Y_FACE_CENTER_X[j][i]<<endl;
	    } 
	}	
	
for(j=1;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
            X_FACE_CENTER_Y[j][i]= Y[j][i];	
	    } 
	}	

//NODE DISTANCE

for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
		DXF[0][j][i] = (X_NODE[j][i]-X_NODE[j][i-1]); //WEST
		DXF[1][j][i] = (Y_NODE[j][i]-Y_NODE[j-1][i]); //SOUTH
		DXF[2][j][i] = (X_NODE[j][i+1]-X_NODE[j][i]); //EAST
		DXF[3][j][i] = (Y_NODE[j+1][i]-Y_NODE[j][i]); //NORTH
		}    
    }
		
//SURFACE AREA

for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{
			S[0][j][i] = (Y[j][i]-Y[j-1][i]); //WEST FACE
			S[1][j][i] = (DELX[j][i]-DELX[j][i-1]);	//SOUTH FACE 
	        S[2][j][i] = (Y[j][i]-Y[j-1][i]); //EAST FACE
	        S[3][j][i] = (DELX[j][i]-DELX[j][i-1]);    //NORTH FACE
       
	        VOL[j][i] = S[2][j][i] * S[3][j][i] ;
	    }    
    }
				
			  
}  

void APPLY_BC(){
	
   	 
	 for(j=1;j<=NJ;j++) // left and right boundary condition
	{
		T[j][1] = 5; //left
		T[j][NI] = -(5/H) * Y_NODE[j][NI] + 15*cos((2*3.14/H)*Y_NODE[j][NI]);//right 300
	}
	
	for(i=1;i<=NI;i++) // left and right boundary condition
	{
		T[1][i] = 15 ; //bottom 
	    T[NJ][i]= 10 ; //top
	}

}

void APPLY_IC(){
	
	for(j=2;j<=NJ-1;j++) // left and right boundary condition
	{
		for(i=2;i<=NI-1;i++){
			
			T[j][i] = 5;
		}
	}
}

void CALC_CONDUCT(){

	for(j=2;j<=NJ-1;j++)
	{
		for(i=2;i<=NI-1;i++)
		{

		  kw[j][i] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[j][i-1]);
          ks[j][i] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[j][i] + Y_FACE_CENTER_X[j][i-1]));
          ke[j][i] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[j][i]);
          kn[j][i] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[j][i] + Y_FACE_CENTER_X[j][i-1]));
       
          SP[j][i] = Q_VOL_GEN*VOL[j][i]/(T_OLD[j][i]);
		    
		  AW[j][i] = kw[j][i] * S[0][j][i] / (DXF[0][j][i]);
		  AS[j][i] = ks[j][i] * S[1][j][i] / (DXF[1][j][i]);
		  AE[j][i] = ke[j][i] *S[2][j][i] / (DXF[2][j][i]);
		  AN[j][i] = kn[j][i] * S[3][j][i] / (DXF[3][j][i]);
		  AP[j][i] = AE[j][i] + AW[j][i] + AN[j][i] + AS[j][i] - SP[j][i];

         // cout<<kw[j][i]<<" "<<ks[j][i]<<" "<<ke[j][i]<<" "<<kn[j][i]<<endl;
       // cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;
		}	
	}
		
	for(j=2;j<=NJ-1;j++)
	{
	   for(i=2;i<=NI-1;i++)
	   {
		T[j][i] = (AE[j][i] * T[j][i+1] + AW[j][i] * T[j][i-1] + AN[j][i] * T[j+1][i] + AS[j][i] * T[j-1][i] ) / AP[j][i];			
	   }
	}
		
}

void UPDATE(){
	for (j=2;j<=NJ-1;j++)
	{
	   for (i=2;i<=NI-1;i++)
	  {	
	   T_OLD[j][i] = T[j][i];
      }
    }
}

double CALC_ABS_ERROR()
{	
/*	for (j=1;j<=NJ;j++)
	{
		for (i=1;i<=NI;i++)
		{	
		  RESIDUE[j][i] =(T[j][i] - T_OLD[j][i]);			
		}
	}
	
//	RTOT = sqrt(RTOT/(NJ*NI));
//	return RTOT;
*/
}


double CALC_RMS_ERROR()
{	
	double RTOT = 0;
	for (j=2;j<=NJ-1;j++)
	{
		for (i=2;i<=NI-1;i++)
		{	
		  RESIDUE[j][i] =(T[j][i] - T_OLD[j][i])*(T[j][i] - T_OLD[j][i]);
          RTOT=RTOT+RESIDUE[j][i];			
		}
	}
	
	RTOT = sqrt(RTOT/((NJ-2)*(NI-2)));
	return RTOT;
}


int FILE_WRITE1(){  // how to write the data in the file ...
	
	//flux in x-direction
    for (j=2;j<=NJ-1;j++)	
	{
	   for (i=1;i<=NI-1;i++)
		{
		 kx[j][i] = 5 * (1 + (100/L)*Y_FACE_CENTER_X[j][i]) ;	
    	 qx[j][i] = -kx[j][i]*(T[j][i+1] - T[j][i])/(X_NODE[j][i+1] - X_NODE[j][i]); 
	    }
     }
	
	//flux in y-direction
	for (j=1;j<=NJ-1;j++)	
	{
	   for (i=2;i<=NI-1;i++)
		{
		 ky[j][i] = 5 * (1 + (100/L)*0.5*(Y_FACE_CENTER_X[j][i] + Y_FACE_CENTER_X[j][i-1])) ;	
		 qy[j][i] = -ky[j][i]*(T[j+1][i] - T[j][i])/(Y_NODE[j+1][i] - Y_NODE[j][i]);	
	    }
    }
	/*
	//grid lines
	ofstream file4("grid_data_42_22.dat");
    file4<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<endl;
    file4<<"ZONE I="<<NI-1<<", J="<<NJ-1<<", F=POINT"<<endl;

	for (j=1;j<=NJ-1;j++)	
		{
			for (i=1;i<=NI-1;i++)
			{
	 			file4<<X[j][i]<<" "<<Y[j][i]<<endl;
			}
		}
			file4.close();
        */
	
    //temperature values at CV and heat flux (face values)
	ofstream writer("structural_grid_42_22.dat");
	if (!writer){
		cout<<"error"<<endl;
	}
	else{
		writer<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<", "<<'"'<<"T"<<'"'<<", "<<'"'<<"qx"<<'"'<<", "<<'"'<<"qy"<<'"'<<endl;
		writer<<"ZONE I="<<NI<<", J="<<NJ<<", F=POINT"<<endl;
     for (j=1;j<=NJ;j++)	
		{
			for (i=1;i<=NI;i++)
			{
	 			writer<<X_NODE[j][i]<<" "<<Y_NODE[j][i]<<" "<<T[j][i]<<" "<<qx[j][i]<<" "<<qy[j][i]<<" "<<endl;
			}
		}
			writer.close();
	}
}

/*
int FILE_WRITE2(int n){  
	int aa;
	char ch[80];
 
	char str_interface[80];

    aa = sprintf (ch, "%d",n);
			
		strcpy (str_interface,"2d_heat_conduction_with_heat_gen");
		strcat (str_interface,ch);
		strcat (str_interface,".dat");
		
  
         ofstream file(str_interface);
		file<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<endl;
		file<<"ZONE I="<<NI<<", J="<<NJ<<", F=POINT"<<endl;
     for(j=1;j<=NJ;j++) 	
		{
			for(i=1;i<=NI;i++)
			{
				file<<X_NODE[j][i]<<" "<<Y_NODE[j][i]<<" "<<endl;
			}
		}
		//	file.close();
	}
*/
