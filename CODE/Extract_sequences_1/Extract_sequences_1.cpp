// Extracts protein sequences from a given PDB file!

#include "stdafx.h"
#include "windows.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <direct.h>
#include <cmath>
#include <conio.h>
#include <vector>

using namespace std;


void transforms(string a, int &br, int x, int y){
    char numb[256]="\0";    
    char i;

    for (i=x;i<x+y;i++)
      numb[i-x]=a[i];
    
    br=atoi(numb);
}




class element{
   public:
      element* NextLinkE(){ return this->link;};
      void ChainLinkE(element* connect){ this->link=connect;};
	  int FindE(char name);
      char prot_name;
      
   private:
      element* link;   
      
};



int element::FindE(char name){
   int res = 1;
   element* prot;
   prot = this;

   while(prot != NULL){

     if (prot->prot_name == name)
		 res = 0;
     
     prot = prot->NextLinkE();
   } 
   
   return res;
}




int main(int argc, char* argv[]){
    string line, pdb_code;    
    char readln[256];
    char chain [256] = "\0";
    int p_residue;  // previous residue added usin Ca atom
    element* psegment = NULL;
    element* FirstE   = NULL;
    element* temp;  
    char file_name[30]; 
    int i, pdb_numb, cur_res;
    int type = 1;


    char* dir1 = "";
    char* dir2 = "sequences";



    ifstream fin(argv[1]);
    p_residue=0;

    printf("--------------------------------------------------\n");

    while (fin.getline(readln,256)){
        line = readln;
        if (line.substr(0,4)=="ATOM"){
            if (line.substr(12,4)==" CA "){
                if (psegment->FindE(readln[21])){
                    psegment = new element;
                    psegment->ChainLinkE(FirstE);
                    psegment->prot_name = readln[21];
                    
                    FirstE = psegment;
                }						 
            }
        }
    }
    fin.close();


    

    temp = psegment;
    while(temp!=NULL){   
        ifstream fin(argv[1]);
        cout << temp->prot_name << endl;

        i = 0;
        file_name[i] = '\0';
        
        do{
            file_name[i] = argv[1][i];
            i++;
        }while (argv[1][i]!='.');
        
        file_name[i] = '_';
        file_name[++i] = temp->prot_name;
        file_name[++i] = '.';
        file_name[++i] = '\0';
        strcat(file_name, "sequence");

        _chdir(dir2);
        ofstream fout(file_name);

        cur_res=-100;

        while (fin.getline(readln,256)){
            line = readln;
            
            
            if ((line.substr(0, 4)=="ATOM")&&(readln[21]==temp->prot_name)&&(line.substr(26,1)==" ")){
                transforms(line, pdb_numb, 22, 4);
                pdb_code = line.substr(17, 3);

                if (pdb_numb!=cur_res){				
                    cur_res=pdb_numb;
                    
                    if (line.substr(17, 3)=="ALA")  pdb_code="A";
                    if (line.substr(17, 3)=="ARG")  pdb_code="R";
                    if (line.substr(17, 3)=="ASP")  pdb_code="D";
                    if (line.substr(17, 3)=="ASN")  pdb_code="N";
                    if (line.substr(17, 3)=="CYS")  pdb_code="C";
                    if (line.substr(17, 3)=="GLU")  pdb_code="E";
                    if (line.substr(17, 3)=="GLN")  pdb_code="Q";
                    if (line.substr(17, 3)=="GLY")  pdb_code="G";
                    if (line.substr(17, 3)=="HIS")  pdb_code="H";
                    if (line.substr(17, 3)=="ILE")  pdb_code="I";
                    if (line.substr(17, 3)=="LEU")  pdb_code="L";
                    if (line.substr(17, 3)=="LYS")  pdb_code="K";
                    if (line.substr(17, 3)=="MET")  pdb_code="M";
                    if (line.substr(17, 3)=="PHE")  pdb_code="F";
                    if (line.substr(17, 3)=="PRO")  pdb_code="P";
                    if (line.substr(17, 3)=="SER")  pdb_code="S";
                    if (line.substr(17, 3)=="THR")  pdb_code="T";
                    if (line.substr(17, 3)=="TRP")  pdb_code="W";
                    if (line.substr(17, 3)=="TYR")  pdb_code="Y";
                    if (line.substr(17, 3)=="VAL")  pdb_code="V";
                    
                    if (type)
                        fout<<pdb_code;
                    else
                        fout<<pdb_numb<<" "<<pdb_code<<endl;
                }
            }
            else 
                if ((line.substr(0,6)=="HETATM")&&(line.substr(21,1)==chain)){
                    transforms(line,pdb_numb,22,4);
                    pdb_code=line.substr(17,3);

                    if (pdb_numb!=cur_res){				
                        cur_res=pdb_numb;
                        //cout<<pdb_numb<<" "<<i++<<endl;
                        if (line.substr(17,3)=="MSE")  pdb_code="M";
                        
                        if (type)
                            fout<<pdb_code;
                        else
                            fout<<pdb_numb<<" "<<pdb_code<<endl;
                    }
                }
        }
        
        fout<<endl;
        fout.close();
        fin.close();

        _chdir("..");

        temp = temp->NextLinkE();
    } 


	return 0;
}
