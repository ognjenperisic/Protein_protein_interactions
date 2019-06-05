// extract_from_dssp_1.cpp

#include "stdafx.h"
#include "windows.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <direct.h>
#include <cmath>
#include <math.h>
#include <conio.h>

using namespace std;

const float min_d = 7.0;

void transforms(string a, float &br, int x, int y){
    char numb[256] = "";    
    char i;

    for (i = x; i<x+y; i++)
      numb[i-x] = a[i];
    	
    br = atof(numb);
}


void transforms(string a, int &br, int x, int y){
    char numb[256] = "";    
    char i;

    for (i = x; i<x+y; i++)
      numb[i-x] = a[i];
    
    br = atoi(numb);
}







int main(int arcg, char* argv[]){
   char readln[256];   
   string line;
   char file_name[30] = "";   
   float waterm;
   int res_numb;
   int pdb_res_numb;
   int res = 0;
   
   char* dir1 = "";
   char* dir2 = "";

   int i = 0;
   int old_pdb;
   int temp_pdb;
   FILE* out;
   BOOL analyze = FALSE;
   
   ifstream fin(argv[1]);
   analyze = FALSE;
   printf("Filename %s\n", argv[1]);
   
   i = 0;
   do{
       file_name[i] = argv[1][i];
       i++;
   }while (argv[1][i] != '.');
   
   file_name[i] = '\0';
   strcat(file_name, ".sur");
   printf("%s\n", file_name);
   
   out = fopen(file_name,"w");
   
   
   old_pdb = -100;
   res = 0;

   while (fin.getline(readln, 256)){
       line = readln;
       if (line.substr(0,15)=="  #  RESIDUE AA"){
           if (!analyze)
               analyze = TRUE;
       }
       else
           if ((line.substr(0,15)!="  #  RESIDUE AA")&(analyze)&(line.substr(13,1)!="!")){							
               transforms(line, waterm, 35, 3);
               transforms(line, res_numb, 0, 5);
               transforms(line, pdb_res_numb, 7, 3);
               if (pdb_res_numb!=old_pdb){
                   res = res + 1; // calculated res number;
                   //printf("%5d %6.2f %6.2f\n",res_numb,waterm,waterm/10);
                   //fprintf(out,"%5d %6.2f\n",res_numb,waterm);
                   if (waterm >= 35)
                       //fprintf(out,"%4d  %d  %4d %5.0f\n",res_numb,1,pdb_res_numb,waterm);
                       //fprintf(out,"%4d  %d  %4d %5.0f\n", res, 1, pdb_res_numb, waterm);
                       fprintf(out,"%4d %4d %5d %5.0f\n", res, pdb_res_numb, 1, waterm);
                   else
                       //fprintf(out,"%4d  %d  %4d %5.0f\n",res_numb,0,pdb_res_numb,waterm);
                       //fprintf(out,"%4d  %d  %4d %5.0f\n", res, 0, pdb_res_numb, waterm);
                       fprintf(out,"%4d %4d %5d %5.0f\n", res, pdb_res_numb, 0, waterm);

                   old_pdb = pdb_res_numb;
               }               
           }
           
   }
   fin.close();
   fclose(out);
      
   return 0;
}





