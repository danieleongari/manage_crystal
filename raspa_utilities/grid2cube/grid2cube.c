/* 	Daniele Ongari - 16/12/2015 (cota) - 20/11/2016 (raspa)
 
	This program reads the .grid file in $RASPA_DIR/share/grids
	extracting all the information contained (it's a binary file)

        The .grid file is computed, written and read in src/grids.c

        usage: grid2cube file1.grid file2.grid file3.grid
        output: file1.cube file2.cube file3.cube
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#define MAXL 1000                        //grid.c
#define MAXP 10                          //grid.c    
#define MAX_NUMBER_OF_PSEUDO_ATOMS 100   //molecule.h
#define REAL double                      //constants.h

#define ANGS2BOHR 1.88973


typedef struct int_vector3 {int x; int y; int z;}    INT_VECTOR3;  //vector.h
typedef struct vector {double x; double y; double z;} VECTOR;  //vector.h

float *****VDWGrid;

int main(int argc, char *argv[] )
{
    int nfiles=1;
    for(nfiles=1; nfiles<argc; nfiles++){ 

    //------------------------ Read input ------------------------------------------------------

    FILE *FilePtr1=fopen(argv[nfiles],"r"); //input file provided in the command line
   
    //(1) Bin size specified
    double SpacingVDWGrid; 
    fread(&SpacingVDWGrid,1,sizeof(SpacingVDWGrid),FilePtr1);

    //(2) Number of gridpoint for each dimension
    static INT_VECTOR3 NumberOfVDWGridPoints;   	
    fread(&NumberOfVDWGridPoints,1,sizeof(INT_VECTOR3),FilePtr1); 

    //NumberOfVDWGridPoints.x=NumberOfVDWGridPoints.x+1;  //DANIELE: done while VDWGrid
    //NumberOfVDWGridPoints.y=NumberOfVDWGridPoints.y+1;
    //NumberOfVDWGridPoints.z=NumberOfVDWGridPoints.z+1;

    //(3) Size of the unit cell
    static VECTOR SizeGrid;
    fread(&SizeGrid,1,sizeof(VECTOR),FilePtr1);

    //(4) How much the cell should be shifted to make a supercell, because of the orthogonalization.
    static VECTOR ShiftGrid;
    fread(&ShiftGrid,1,sizeof(VECTOR),FilePtr1);

    //(5) Exact binsize computed
    static VECTOR DeltaVDWGrid;
    fread(&DeltaVDWGrid,1,sizeof(VECTOR),FilePtr1);

    //(6) UnitCellSize=SizeGrid 
    static VECTOR unit_cell_size; 
    fread(&unit_cell_size,1,sizeof(VECTOR),FilePtr1);

    //(7) Number of unit cells (es: 1 1 1)
    static INT_VECTOR3 number_of_unit_cells;
    fread(&number_of_unit_cells,1,sizeof(INT_VECTOR3),FilePtr1);

    //(8) After the VDW energy (m=0), the first (m=1,2,3), second (m=4,5,6) and third (m=7) derivatives are stored in the .grid file in units of (J/mol)
    int i, j, k, m;
    int l=0;
    int NumberOfGrids=1;
    int NumberOfPseudoAtoms=1;
    int *GridTypeList;
    //float *****VDWGrid;
    GridTypeList=(int*)calloc(NumberOfGrids,sizeof(int)); //this is assigned for grids of different atoms. Unuseful here
    GridTypeList[l]=0;  
    VDWGrid=(float*****)calloc(NumberOfPseudoAtoms,sizeof(float****));

    VDWGrid[GridTypeList[l]]=(float****)calloc(NumberOfVDWGridPoints.x+1,sizeof(float***));
    for(i=0;i<=NumberOfVDWGridPoints.x;i++)
    {
      VDWGrid[GridTypeList[l]][i]=(float***)calloc(NumberOfVDWGridPoints.y+1,sizeof(float**));
      for(j=0;j<=NumberOfVDWGridPoints.y;j++)
      {
        VDWGrid[GridTypeList[l]][i][j]=(float**)calloc(NumberOfVDWGridPoints.z+1,sizeof(float*));
        for(k=0;k<=NumberOfVDWGridPoints.z;k++)
          VDWGrid[GridTypeList[l]][i][j][k]=(float*)calloc(8,sizeof(float));
      }
    }

  //for(m=0;m<8;m++) //no need for derivatives
    m=0; 
      for(i=0;i<=NumberOfVDWGridPoints.x;i++)
        for(j=0;j<=NumberOfVDWGridPoints.y;j++)
          for(k=0;k<=NumberOfVDWGridPoints.z;k++)
            fread(&VDWGrid[GridTypeList[l]][i][j][k][m],1,sizeof(float),FilePtr1);
    
    fclose(FilePtr1); 
    
    //------------------------ Write output ------------------------------------------------------


    char *c = strdup(argv[nfiles]);
    int  len = strlen(c);
    c[len-4]='c';c[len-3]='u';c[len-2]='b';c[len-1]='e';   //I don't know a more elegant way to do it!
    FILE *FilePtr2=fopen(c,"w"); //output file provided in the command line


    fprintf(FilePtr2,"# GENERATED FROM: %s  (remember: kJ/mol units used)\n", argv[nfiles]);
    fprintf(FilePtr2,"# Cube file with only VdW grid builted from raspa with: SpacingVDWGrid %f NumberOfVDWGridPoints %d %d %d SizeGrid %f %f %f ShiftGrid %f %f %f DeltaVDWGrid %f %f %f unit_cell_size %f %f %f number_of_unit_cells %d %d %d \n",
                                                                                 SpacingVDWGrid,           
                                                      NumberOfVDWGridPoints.x,NumberOfVDWGridPoints.y,NumberOfVDWGridPoints.z,        
                                                                                                                SizeGrid.x,SizeGrid.y,SizeGrid.z,
                                                                                                                                ShiftGrid.x,ShiftGrid.y,ShiftGrid.z,
                                                                                                                                             DeltaVDWGrid.x,DeltaVDWGrid.y,DeltaVDWGrid.z,                                                
                                                                                                                                                                unit_cell_size.x,unit_cell_size.y,unit_cell_size.z,                                      
                                                                                                                                                                           number_of_unit_cells.x,number_of_unit_cells.y,number_of_unit_cells.z);

    fprintf(FilePtr2,"%6i %12.6f %12.6f %12.6f\n",1,ShiftGrid.x*ANGS2BOHR,ShiftGrid.y*ANGS2BOHR,ShiftGrid.z*ANGS2BOHR); // natoms, origin x,y,z

    //Number of grid points in x direction: 
    fprintf(FilePtr2,"%6i %12.6f %12.6f %12.6f \n",(NumberOfVDWGridPoints.x),
                                                   SizeGrid.x/(double)(NumberOfVDWGridPoints.x)*ANGS2BOHR,  //binwidth in bohr 
                                                   0*ANGS2BOHR,              
                                                   0*ANGS2BOHR);
    //Number of grid points in y direction:
    fprintf(FilePtr2,"%6i %12.6f %12.6f %12.6f \n",(NumberOfVDWGridPoints.y),
                                                   0*ANGS2BOHR, 
                                                   SizeGrid.y/(double)(NumberOfVDWGridPoints.y)*ANGS2BOHR, 
                                                   0*ANGS2BOHR);
    // Number of grid points in z direction:
    fprintf(FilePtr2,"%6i %12.6f %12.6f %12.6f \n",(NumberOfVDWGridPoints.z), 
                                                   0*ANGS2BOHR, 
                                                   0*ANGS2BOHR,
                                                   SizeGrid.z/(double)(NumberOfVDWGridPoints.z)*ANGS2BOHR);  

    //One fake atom needed in the cube file to be read by VMD
    fprintf(FilePtr2,"%6i %12.6f %12.6f %12.6f",40,0.0,0.0,0.0); //no /n: it will be added by the next loop

    int x,y,z,column;
    double value;
    double ENERGY_TO_KJ_PER_MOL=0.01;

    for(x=0; x<NumberOfVDWGridPoints.x; x++)
     {
     for(y=0; y<NumberOfVDWGridPoints.y; y++)
      {
      column=1;
      for(z=0; z<NumberOfVDWGridPoints.z; z++)
	{
       	value=VDWGrid[0][x][y][z][0]*ENERGY_TO_KJ_PER_MOL; //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> (kJ/mol)
        if (column==1){
          fprintf(FilePtr2,"\n%.6E",value);
          column++;
          }
        else if (column==6){
          fprintf(FilePtr2,"  %.6E",value);
          column=1;
          }
        else{
          fprintf(FilePtr2,"  %.6E",value);
          column++;
          } 
              		
         }
      }
   }

    fclose(FilePtr2);

    } //processed all nfiles
     
     
    return 0;
}




