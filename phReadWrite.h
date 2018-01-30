#ifndef __phReadWrite_H_
#define __phReadWrite_H_

#include <strings.h>
#include <time.h>
#ifdef __cplusplus
extern "C" {
#endif

  // to read parameters from a phasta file (filename)
  // parameters correspond to nshg & nvar, i.e., size of field-array
  // these parameters are used as reference values 
  // (sometimes needed before reading the field-array)
  void readParametersFromFile(char *filename,
			      char *fieldName,
			      int &nshg, 
			      int &numVars);
  
  // to read array from a phasta file (filename)
  // memory is allocated HERE for 'valueArray'
  // `fieldName' tells which block to read like solution, error etc.
  void readArrayFromFile( char *filename,
			  char *fieldName,
			  double *&valueArray);

  // to write array to a phasta file (filename)
  // NOTE: array should be transposed!!!
  // `fieldName' tells in which block to write like solution, error etc.
  // `outputFormat' tells in which format to write, i.e., binary/ascii
  // `mode' : "write", "append" etc.
  void writeArrayToFile( char *filename,
			 char *fieldName,
			 char *outputFormat,
			 char *mode,
			 int nshg,
			 int numVars,
			 int stepNumber,
			 double *valueArray);
  
#ifdef __cplusplus
}
#endif


#endif
