#ifndef __CTQMC_WORLD__
#define __CTQMC_WORLD__

#include <iostream> 
#include <fstream>

#include "mpi.h"

#include <complex>
#include <vector>
#include "config.h"



#define MatrixType_SEND_TAG 60
#define RealMatrixType_SEND_TAG 61
#define VectorType_SEND_TAG 70
#define RealVectorType_SEND_TAG 51
#define INT_SEND_TAG 50
#define RealType_SEND_TAG 51
#define ScalarType_SEND_TAG 52
#define STRING_SEND_TAG 80

#define I_INT_SEND_TAG 150

#define RIVER_STREAM_SIZE 100

class CTQMC_WORLD {

/*
  A class for all ct-qmc world processes.
  This initiates communicators and makes all send/receive of common objects like matrices, vectors, numbers
 */

ofstream rivers;  // Output stream
MPI_Comm CTQMC_COMM; // Communicator for CT-QMC Processes
int WORLD_SIZE; // Amount of processes we use
MPI_Status CT_STATUS; // Status of receiving message
char **river_names; // The output streams for each process

public :
CTQMC_WORLD(){}; // Empty constructor 

void INIT(int argc, char **argv, int N);
void INIT(int argc, char **argv);
/*-----------------------------------------------------------------------------------*/
/**
  Constructor of world - Done by all processes. It initiates a communicator, the size of the world
  Also initializes all output streams
  @N - the amount of processes we want to initialize. If N is not defined, the N is taken from command prompt
 */
/*-----------------------------------------------------------------------------------*/

int getRank();
/*-----------------------------------------------------------------------------------*/
/**
  Returns rank of a process in CTQMC_WORLD
 */
/*-----------------------------------------------------------------------------------*/

int getSize();
/*-----------------------------------------------------------------------------------*/
/**
  Returns amount of processes in CTQMC_WORLD
 */
/*-----------------------------------------------------------------------------------*/



void IniProcess();
/*-----------------------------------------------------------------------------------*/
/**
  Routines which should be done with start of each process, such as opening its stream
  */
/*-----------------------------------------------------------------------------------*/

ofstream& getStream();
/*-----------------------------------------------------------------------------------*/
/**
  Returns stream of current process
 */
/*-----------------------------------------------------------------------------------*/

// Catching matrix operations 

MatrixType &catch_MatrixType(const int SENDER, int &sizex, int &sizey); 
MatrixType &catch_MatrixType(const int SENDER);

/*-----------------------------------------------------------------------------------*/
/** 
  Catches matrix from process SENDER and writes its dimensions to sizex and sizey. 
  Returns a caught MatrixType
  @SENDER - the id of a process which sends a MatrixType
  @sizex - the x-size (number of rows) of a MatrixType
  @sizey - the y-size (number of cols) of a MatrixType
*/
/*-----------------------------------------------------------------------------------*/

void catch_MatrixType(MatrixType& result, const int SENDER, int &sizex, int &sizey); 
void catch_MatrixType(MatrixType& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a matrix (sizex * sizey) from process "SENDER"
  @result - caught matrix
  @SENDER - the id of a process which sends a MatrixType
  @sizex - the x-size (number of rows) of a MatrixType
  @sizey - the y-size (number of cols) of a MatrixType
 */
/*-----------------------------------------------------------------------------------*/
RealMatrixType &catch_RealMatrixType(const int SENDER, int &sizex, int &sizey); 
RealMatrixType &catch_RealMatrixType(const int SENDER);

/*-----------------------------------------------------------------------------------*/
/** 
  Catches RealMatrixType from process SENDER and writes its dimensions to sizex and sizey. 
  Returns a caught RealMatrixType
  @SENDER - the id of a process which sends a RealMatrixType
  @sizex - the x-size (number of rows) of a RealMatrixType
  @sizey - the y-size (number of cols) of a RealMatrixType
*/
/*-----------------------------------------------------------------------------------*/

void catch_RealMatrixType(RealMatrixType& result, const int SENDER, int &sizex, int &sizey); 
void catch_RealMatrixType(RealMatrixType& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a RealMatrixType (sizex * sizey) from process "SENDER"
  @result - caught RealMatrixType
  @SENDER - the id of a process which sends a RealMatrixType
  @sizex - the x-size (number of rows) of a RealMatrixType
  @sizey - the y-size (number of cols) of a RealMatrixType
 */
/*-----------------------------------------------------------------------------------*/

VectorType &catch_VectorType(const int SENDER, int &size); 
VectorType &catch_VectorType(const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a VectorType (size) from process "SENDER"
  @SENDER - the id of a process which sends a VectorType
  @size - the size (number of cols) of a VectorType
 */
/*-----------------------------------------------------------------------------------*/

void catch_VectorType(VectorType& result, const int SENDER, int &size); 
void catch_VectorType(VectorType& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a VectorType (size) from process "SENDER"
  @result - caught VectorType
  @SENDER - the id of a process which sends a VectorType
  @size - the size (number of cols) of a VectorType
 */
/*-----------------------------------------------------------------------------------*/

RealVectorType &catch_RealVectorType(const int SENDER, int &size); 
RealVectorType &catch_RealVectorType(const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a RealVectorType (size) from process "SENDER"
  @SENDER - the id of a process which sends a RealVectorType
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/

void catch_RealVectorType(RealVectorType& result, const int SENDER, int &size); 
void catch_RealVectorType(RealVectorType& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a RealVectorType (size) from process "SENDER"
  @result - caught RealVectorType
  @SENDER - the id of a process which sends a RealVectorType
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/

IntVectorType &catch_IntVectorType(const int SENDER, int &size); 
IntVectorType &catch_IntVectorType(const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive an IntVectorType (size) from process "SENDER"
  @SENDER - the id of a process which sends an IntVectorType
  @size - the size (number of cols) of an IntVectorType
 */
/*-----------------------------------------------------------------------------------*/

void catch_IntVectorType(IntVectorType& result, const int SENDER, int &size); 
void catch_IntVectorType(IntVectorType& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a IntVectorType (size) from process "SENDER"
  @result - caught IntVectorType
  @SENDER - the id of a process which sends a IntVectorType
  @size - the size (number of cols) of a IntVectorType
 */
/*-----------------------------------------------------------------------------------*/

int catch_int(const int SENDER); 
void catch_int(int& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a single integer from process "SENDER"
  @result - caught integer
  @SENDER - the id of a process which sends an integer
 */
/*-----------------------------------------------------------------------------------*/

RealType catch_float(const int SENDER); 
void catch_float(RealType& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive an single RealType from process "SENDER"
  @result - caught RealType
  @SENDER - the id of a process which sends an RealType
 */
/*-----------------------------------------------------------------------------------*/

ScalarType catch_complex(const int SENDER); 
void catch_complex(ScalarType& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a single ScalarType from process "SENDER"
  @result - caught ScalarType
  @SENDER - the id of a process which sends a ScalarType
 */
/*-----------------------------------------------------------------------------------*/

string& catch_string(const int SENDER, int &size); 
string& catch_string(const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a string (size) from process "SENDER"
  @SENDER - the id of a process which sends a string
  @size - the size (number of rows) of a string
 */
/*-----------------------------------------------------------------------------------*/

void catch_string(string& result, const int SENDER, int &size); 
void catch_string(string& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a string (size) from process "SENDER"
  @result - caught string
  @SENDER - the id of a process which sends a string
  @size - the size (number of rows) of a string
 */
/*-----------------------------------------------------------------------------------*/

/*template <> MatrixType &recv(const int SENDER);
template <> RealMatrixType &recv(const int SENDER);
template <> VectorType &recv(const int SENDER);
template <> RealVectorType &recv(const int SENDER);
template <> int &recv(const int SENDER);
template <> RealType &recv(const int SENDER);
template <> ScalarType &recv(const int SENDER);
template <> string &recv(const int SENDER);*/

void recv (MatrixType &result, const int SENDER);
void recv (RealMatrixType &result, const int SENDER); 
void recv (VectorType &result, const int SENDER);
void recv (RealVectorType &result, const int SENDER);
void recv (int& result, const int SENDER);
void recv (RealType& result, const int SENDER);
void recv (ScalarType& result, const int SENDER);
void recv (string& result, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/**
  Covers for all receiving functions
  @SENDER - the id of a sending process
  @x - the variable to be received
 */
/*-----------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------*/
void send_MatrixType (const int RECEIVER, MatrixType &x, int sizex, int sizey);
void send_MatrixType (const int RECEIVER, MatrixType &x); 
/*-----------------------------------------------------------------------------------*/
/**
  Send MatrixType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a MatrixType
  @x - the MatrixType to be sent
  @sizex - the x-size (number of rows) of a MatrixType
  @sizey - the y-size (number of cols) of a MatrixType
  
 */
/*-----------------------------------------------------------------------------------*/

void send_RealMatrixType (const int RECEIVER, RealMatrixType &x, int sizex, int sizey);
void send_RealMatrixType (const int RECEIVER, RealMatrixType &x); 
/*-----------------------------------------------------------------------------------*/
/**
  Send RealMatrixType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a RealMatrixType
  @x - the RealMatrixType to be sent
  @sizex - the x-size (number of rows) of a RealMatrixType
  @sizey - the y-size (number of cols) of a RealMatrixType
  
 */
/*-----------------------------------------------------------------------------------*/

void send_VectorType (const int RECEIVER, VectorType &x, int size);
void send_VectorType (const int RECEIVER, VectorType &x); 
/*-----------------------------------------------------------------------------------*/
/**
  Send VectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a VectorType
  @x - the VectorType to be sent
  @size - the size (number of cols) of a VectorType
 */
/*-----------------------------------------------------------------------------------*/


void send_RealVectorType (const int RECEIVER, RealVectorType &x, int size);
void send_RealVectorType (const int RECEIVER, RealVectorType &x); 
/*-----------------------------------------------------------------------------------*/
/**
  Send RealVectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a RealVectorType
  @x - the RealVectorType to be sent
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/
void send_IntVectorType (const int RECEIVER, IntVectorType &x, int size);
void send_IntVectorType (const int RECEIVER, IntVectorType &x); 
/*-----------------------------------------------------------------------------------*/
/**
  Send RealVectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a IntVectorType
  @x - the RealVectorType to be sent
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/

void send_int(const int RECEIVER, int x);
/*-----------------------------------------------------------------------------------*/
/**
  Send integer x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of an integer
  @x - the integer to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send_float(const int RECEIVER, RealType x);
/*-----------------------------------------------------------------------------------*/
/**
  Send RealType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of an RealType
  @x - the RealType to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send_complex(const int RECEIVER, ScalarType x);
/*-----------------------------------------------------------------------------------*/
/**
  Send integer x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a ScalarType
  @x - the ScalarType to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send_string(const int RECEIVER, string &x, int size);
void send_string(const int RECEIVER, string &x);

/*-----------------------------------------------------------------------------------*/
/**
  Send string x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a string
  @x - the string to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send (const int RECEIVER, MatrixType &x);
void send (const int RECEIVER, RealMatrixType &x); 
void send (const int RECEIVER, VectorType &x);
void send (const int RECEIVER, RealVectorType &x);
void send (const int RECEIVER, IntVectorType &x);
void send (const int RECEIVER, int x);
void send (const int RECEIVER, RealType x);
void send (const int RECEIVER, ScalarType x);
void send (const int RECEIVER, string &x);
/*-----------------------------------------------------------------------------------*/
/**
  The cover for all sending procedures - sends x of mentioned types to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a MatrixType
  @x - the variable to be sent
 */
/*-----------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------*/
/*-------------------------NON-BLOCKING ROUTINES-------------------------------------*/
/*-----------------------------------------------------------------------------------*/


MPI_Request& isend_int(const int RECEIVER, int x);
MPI_Request& isend(const int RECEIVER, int x);
/*-----------------------------------------------------------------------------------*/
/**
  Immediate send of integer x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient
  @x - the integer to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& irecv_int(int &x, const int SENDER);
MPI_Request& irecv(int &x, const int SENDER);
/*-----------------------------------------------------------------------------------*/
/**
  Initialize and perform non-blocking receive of integer x from process "SENDER"
  @SENDER - the id of a process - recipient
  @x - the integer to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/



//void isend_raw(const int RECEIVER, vector<ScalarType> &x, int size);
//void isend_raw(const int RECEIVER, vector<RealType> &x, int size);
/*-----------------------------------------------------------------------------------*/
/**
  Immediate send of raw data x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a string
  @x - the string to be sent
  Returns 2 requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

//void irecv_raw(const int SENDER, vector<ScalarType> &x, int& size);
//void irecv_raw(const int SENDER, vector<RealType> &x, int& size);
/*-----------------------------------------------------------------------------------*/
/**
  Non-blocking receive of raw data x from process "SENDER"
  @SENDER - the id of a process - recipient of a string
  @x - the string to be sent
  Returns 2 requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

void wait(MPI_Request &REQ);
int test(MPI_Request &REQ);
/*-----------------------------------------------------------------------------------*/
/**
  Blocking and non-blocking test procedures for non-blocking sending/receiving raw data
  @REQ - the request array of 2 elements to be checked
 */
/*-----------------------------------------------------------------------------------*/

~CTQMC_WORLD();
/*-----------------------------------------------------------------------------------*/
/**
  Destructor - Finishes all MPI Processes
  */
/*-----------------------------------------------------------------------------------*/

};

/*-----------------------------------------------------------------------------------*/


#endif // endif :: __CTQMC_WORLD__


