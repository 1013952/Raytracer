#ifndef _MATRIXHEADERS
#define _MATRIXHEADERS

#include<iostream>

struct matrix3
{
	double els[3][3];
	void eromult(double a, int i){
		if(a!=0 && 0<=i && i<3)
			for(int j=0; j<3; j++)
				els[i][j] *= a;
	}
	void eroswap(int i, int j){
		if(0<=i && 0<=j && i<3 && j<3 && i!=j)
			for(int l=0; l<3; l++){
				double d = els[i][l];
				els[i][l] = els[j][l];
				els[j][l] = d;
			}
	}
	void eroadd(int i, double a, int j){
		if(0<=i && 0<=j && i<3 && j<3 && i!=j)
			for(int l=0; l<3; l++)
				els[j][l] += a*els[i][l];
	}
	void rref(){
		int j=0; int i =0;
		while(i < 3){
			//checks for top-left nonzero element
			while(els[i][j]==0 && i<3){
				j++;
				if(j==3){
					j = 0; i++;
				}
			}//els[i][j] is nonzero
			eromult((1/els[i][j]), i);
			for(int l=0; l<3; l++){
				if(l!=i)
				eroadd(i, -els[l][j], l);
			}
			i++;
		}
	}
	void printmat(){
		for(int i=0; i<3; i++){
		for(int j=0; j<3; j++)
			std::cout<<els[i][j]<<" ";
		std::cout<<std::endl;
	}
	}
};

struct mat34
{
	double els[3][4];
	void eromult(double a, int i){
		if(a!=0 && 0<=i && i<3)
			for(int j=0; j<4; j++)
				els[i][j]*=a;
	}
	void eroswap(int i, int l){
		if(0<=i && 0<=l && i<3 && l<3){
			for(int j=0; j<4; j++){
				double d = els[i][j];
				els[i][j] = els[l][j];
				els[l][j] = els[i][j];
			}
		}
	}
	void eroadd(int i, double a, int l){
		if(0<=i && 0<=l && i<3 && l<3 && i!=l){
			for(int j=0; j<4; j++)
				els[l][j]+= a*els[i][j];
		}
	}
	void rref(){
		int j=0; int i =0;
		while(i < 4){
			//checks for top-left nonzero element
			while(els[i][j]==0 && i<3){
				j++;
				if(j==3){
					j = 0; i++;
				}
			}//els[i][j] is nonzero
			eromult((1/els[i][j]), i);
			for(int l=0; l<3; l++){
				if(l!=i)
				eroadd(i, -els[l][j], l);
			}
			i++;
		}
	}
	void printmat(){
		for(int i=0; i<3; i++){
		for(int j=0; j<4; j++)
			std::cout<<els[i][j]<<" ";
		std::cout<<std::endl;
	}
	}
};


#endif