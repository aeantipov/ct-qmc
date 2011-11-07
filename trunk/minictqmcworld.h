/* This file is a part of 'ct-qmc 1.6' project
**
** minictqmcworld.h
** Declaration of the CTQMC_WORLD class. Doesn't support matrix definitions
**
** Author    		     : Andrey Antipov, antipov@shg.ru
*/

#ifndef __CTQMC_WORLD__
#define __CTQMC_WORLD__

#include <iostream> 
#include <fstream>

#include "mpi.h"

#include <complex>
#include <vector>
#include "config.h"
#include <string.h>

/* Instruction for MPI to sum complex values */
void CTQMC_MPI_COMPLEX_ADD ( ComplexType* invec, ComplexType *inoutvec, int *len, MPI_Datatype *dtype);

/* TAGS FOR CORRESPONDING TYPES */

class CTQMC_TAG {

/*
  A class for tags for sending data.
  Usually pre-defined tags are used. But when some complex data or signals need to be sent - custom tags may be useful
 */
int tag_;
public : 
  CTQMC_TAG (int x) : tag_(x) {};
  CTQMC_TAG () {};
  CTQMC_TAG& operator = (const int x) {(*this).tag_ = x; return *this;};
  CTQMC_TAG& operator = (const CTQMC_TAG &RHS) {(*this).tag_ = RHS.tag_; return *this;};
  operator int () { return tag_;}
};


static CTQMC_TAG MatrixType_SEND_TAG = 60;
static CTQMC_TAG RealMatrixType_SEND_TAG = 61;
static CTQMC_TAG VectorType_SEND_TAG = 70;
static CTQMC_TAG RealVectorType_SEND_TAG = 71;
static CTQMC_TAG IntVectorType_SEND_TAG = 72;
static CTQMC_TAG INT_SEND_TAG = 50;
static CTQMC_TAG n_type_SEND_TAG = 51;
static CTQMC_TAG ComplexType_SEND_TAG = 52;
static CTQMC_TAG STRING_SEND_TAG = 81;

static CTQMC_TAG I_INT_SEND_TAG = 150;
static CTQMC_TAG I_FLOAT_SEND_TAG = 151;
static CTQMC_TAG I_ComplexType_SEND_TAG = 152;
static CTQMC_TAG I_VectorType_SEND_TAG = 170;
static CTQMC_TAG I_RealVectorType_SEND_TAG = 171;
static CTQMC_TAG I_IntVectorType_SEND_TAG = 172;

static CTQMC_TAG RAW_COMPLEX_1D_SEND_TAG = 200;
static CTQMC_TAG RAW_COMPLEX_2D_SEND_TAG = 201;
static CTQMC_TAG RAW_COMPLEX_3D_SEND_TAG = 202;
static CTQMC_TAG RAW_COMPLEX_4D_SEND_TAG = 203;

#define RIVER_STREAM_SIZE 100

class CTQMC_WORLD {

/*
  A class for all ct-qmc world processes.
  This initiates communicators and makes all send/receive of common objects like matrices, vectors, numbers
 */

ofstream rivers;  // Output stream
int WORLD_SIZE; // Amount of processes we use
MPI_Status CT_STATUS; // Status of receiving message
char **river_names; // The output streams for each process

public :
MPI_Datatype CTQMCComplex;
MPI_Comm CTQMC_COMM; // Communicator for CT-QMC Processes
MPI_Op ComplexAdd; // A summation MPI Operation
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

VectorType &recv_VectorType(const int SENDER, int &size, CTQMC_TAG tag = VectorType_SEND_TAG); 
VectorType &recv_VectorType(const int SENDER, CTQMC_TAG tag = VectorType_SEND_TAG);
void recv_VectorType(const int SENDER, VectorType& result, int &size, CTQMC_TAG tag = VectorType_SEND_TAG); 
void recv_VectorType(const int SENDER, VectorType& result, CTQMC_TAG tag = VectorType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a VectorType (size) from process "SENDER"
  @result - caught VectorType
  @SENDER - the id of a process which sends a VectorType
  @size - the size (number of cols) of a VectorType
 */
/*-----------------------------------------------------------------------------------*/

RealVectorType &recv_RealVectorType(const int SENDER, int &size, CTQMC_TAG tag = RealVectorType_SEND_TAG); 
RealVectorType &recv_RealVectorType(const int SENDER, CTQMC_TAG tag = RealVectorType_SEND_TAG);
void recv_RealVectorType(const int SENDER, RealVectorType& result, int &size, CTQMC_TAG tag = RealVectorType_SEND_TAG); 
void recv_RealVectorType(const int SENDER, RealVectorType& result, CTQMC_TAG tag = RealVectorType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a RealVectorType (size) from process "SENDER"
  @result - caught RealVectorType
  @SENDER - the id of a process which sends a RealVectorType
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/

IntVectorType &recv_IntVectorType(const int SENDER, int &size, CTQMC_TAG tag = IntVectorType_SEND_TAG); 
IntVectorType &recv_IntVectorType(const int SENDER, CTQMC_TAG tag = IntVectorType_SEND_TAG);
void recv_IntVectorType(const int SENDER, IntVectorType& result, int &size, CTQMC_TAG tag = IntVectorType_SEND_TAG); 
void recv_IntVectorType(const int SENDER, IntVectorType& result, CTQMC_TAG tag = IntVectorType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a IntVectorType (size) from process "SENDER"
  @result - caught IntVectorType
  @SENDER - the id of a process which sends a IntVectorType
  @size - the size (number of cols) of a IntVectorType
 */
/*-----------------------------------------------------------------------------------*/

int recv_int(const int SENDER, CTQMC_TAG tag = INT_SEND_TAG); 
void recv_int(const int SENDER, int& result, CTQMC_TAG tag = INT_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a single integer from process "SENDER"
  @result - caught integer
  @SENDER - the id of a process which sends an integer
 */
/*-----------------------------------------------------------------------------------*/

n_type recv_float(const int SENDER, CTQMC_TAG tag = n_type_SEND_TAG); 
void recv_float(const int SENDER, n_type& result, CTQMC_TAG tag = n_type_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive an single n_type from process "SENDER"
  @result - caught n_type
  @SENDER - the id of a process which sends an n_type
 */
/*-----------------------------------------------------------------------------------*/

ComplexType recv_complex(const int SENDER, CTQMC_TAG tag = ComplexType_SEND_TAG); 
void recv_complex(const int SENDER, ComplexType& result, CTQMC_TAG tag = ComplexType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a single ComplexType from process "SENDER"
  @result - caught ComplexType
  @SENDER - the id of a process which sends a ComplexType
 */
/*-----------------------------------------------------------------------------------*/

string& recv_string(const int SENDER, int &size, CTQMC_TAG tag = STRING_SEND_TAG); 
string& recv_string(const int SENDER, CTQMC_TAG tag = STRING_SEND_TAG);
void recv_string(const int SENDER, string& result,  int &size, CTQMC_TAG tag = STRING_SEND_TAG); 
void recv_string(const int SENDER, string& result, CTQMC_TAG tag = STRING_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a string (size) from process "SENDER"
  @result - caught string
  @SENDER - the id of a process which sends a string
  @size - the size (number of rows) of a string
 */
/*-----------------------------------------------------------------------------------*/

void recv (const int SENDER, VectorType &result, CTQMC_TAG tag = VectorType_SEND_TAG);
void recv (const int SENDER, RealVectorType &result, CTQMC_TAG tag = RealVectorType_SEND_TAG);
void recv (const int SENDER, IntVectorType &result, CTQMC_TAG tag = IntVectorType_SEND_TAG);
void recv (const int SENDER, int& result, CTQMC_TAG tag = INT_SEND_TAG);
void recv (const int SENDER, n_type& result, CTQMC_TAG tag = n_type_SEND_TAG);
void recv (const int SENDER, ComplexType& result, CTQMC_TAG tag = ComplexType_SEND_TAG);
void recv (const int SENDER, string& result, CTQMC_TAG tag = STRING_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Covers for all receiving functions
  @SENDER - the id of a sending process
  @x - the variable to be received
 */
/*-----------------------------------------------------------------------------------*/

void send_VectorType (const int RECEIVER, VectorType &x, int size, CTQMC_TAG tag = VectorType_SEND_TAG);
void send_VectorType (const int RECEIVER, VectorType &x, CTQMC_TAG tag = VectorType_SEND_TAG); 
/*-----------------------------------------------------------------------------------*/
/**
  Send VectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a VectorType
  @x - the VectorType to be sent
  @size - the size (number of cols) of a VectorType
 */
/*-----------------------------------------------------------------------------------*/


void send_RealVectorType (const int RECEIVER, RealVectorType &x, int size, CTQMC_TAG tag = RealVectorType_SEND_TAG);
void send_RealVectorType (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag = RealVectorType_SEND_TAG); 
/*-----------------------------------------------------------------------------------*/
/**
  Send RealVectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a RealVectorType
  @x - the RealVectorType to be sent
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/
void send_IntVectorType (const int RECEIVER, IntVectorType &x, int size, CTQMC_TAG tag = IntVectorType_SEND_TAG);
void send_IntVectorType (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag = IntVectorType_SEND_TAG); 
/*-----------------------------------------------------------------------------------*/
/**
  Send RealVectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a IntVectorType
  @x - the RealVectorType to be sent
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/

void send_int(const int RECEIVER, int x, CTQMC_TAG tag = INT_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Send integer x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of an integer
  @x - the integer to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send_float(const int RECEIVER, n_type x, CTQMC_TAG tag = n_type_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Send n_type x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of an n_type
  @x - the n_type to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send_complex(const int RECEIVER, ComplexType x, CTQMC_TAG tag = ComplexType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Send integer x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a ComplexType
  @x - the ComplexType to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send_string(const int RECEIVER, string &x, int size, CTQMC_TAG tag = STRING_SEND_TAG);
void send_string(const int RECEIVER, string &x, CTQMC_TAG tag = STRING_SEND_TAG);

/*-----------------------------------------------------------------------------------*/
/**
  Send string x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a string
  @x - the string to be sent
 */
/*-----------------------------------------------------------------------------------*/

void send (const int RECEIVER, VectorType &x, CTQMC_TAG tag = VectorType_SEND_TAG);
void send (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag = RealVectorType_SEND_TAG);
void send (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag = IntVectorType_SEND_TAG);
void send (const int RECEIVER, int x, CTQMC_TAG tag = INT_SEND_TAG);
void send (const int RECEIVER, n_type x, CTQMC_TAG tag = n_type_SEND_TAG);
void send (const int RECEIVER, ComplexType x, CTQMC_TAG tag = ComplexType_SEND_TAG);
void send (const int RECEIVER, string &x, CTQMC_TAG tag = STRING_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  The cover for all sending procedures - sends x of mentioned types to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a MatrixType
  @x - the variable to be sent
 */
/*-----------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------*/
/*----------------------NON-BLOCKING SENDING ROUTINES--------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& isend_VectorType (const int RECEIVER, VectorType &x, int size, CTQMC_TAG tag = I_VectorType_SEND_TAG);
MPI_Request& isend_VectorType (const int RECEIVER, VectorType &x, CTQMC_TAG tag = I_VectorType_SEND_TAG); 
MPI_Request& isend (const int RECEIVER, VectorType &x, CTQMC_TAG tag = I_VectorType_SEND_TAG); 
/*-----------------------------------------------------------------------------------*/
/**
  Non-blocking sending VectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a VectorType
  @x - the VectorType to be sent
  @size - the size (number of cols) of a VectorType
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& isend_RealVectorType (const int RECEIVER, RealVectorType &x, int size, CTQMC_TAG tag = I_RealVectorType_SEND_TAG);
MPI_Request& isend_RealVectorType (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag = I_RealVectorType_SEND_TAG); 
MPI_Request& isend (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag = I_RealVectorType_SEND_TAG); 
/*-----------------------------------------------------------------------------------*/
/**
  Non-blocking sending RealVectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a RealVectorType
  @x - the RealVectorType to be sent
  @size - the size (number of cols) of a RealVectorType
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/
MPI_Request& isend_IntVectorType (const int RECEIVER, IntVectorType &x, int size, CTQMC_TAG tag = I_IntVectorType_SEND_TAG);
MPI_Request& isend_IntVectorType (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag = I_IntVectorType_SEND_TAG); 
MPI_Request& isend (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag = I_IntVectorType_SEND_TAG); 
/*-----------------------------------------------------------------------------------*/
/**
  Non-blocking sending RealVectorType x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a IntVectorType
  @x - the RealVectorType to be sent
  @size - the size (number of cols) of a RealVectorType
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/


MPI_Request& isend_int(const int RECEIVER, int x, CTQMC_TAG tag = I_INT_SEND_TAG);
MPI_Request& isend(const int RECEIVER, int x, CTQMC_TAG tag = I_INT_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Non-blocking sending of integer x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient
  @x - the integer to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& isend_float(const int RECEIVER, n_type x, CTQMC_TAG tag = I_FLOAT_SEND_TAG);
MPI_Request& isend(const int RECEIVER, n_type x, CTQMC_TAG tag = I_FLOAT_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Non-blocking sending of n_type x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of an n_type
  @x - the n_type to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& isend_complex(const int RECEIVER, ComplexType x, CTQMC_TAG tag = I_ComplexType_SEND_TAG);
MPI_Request& isend(const int RECEIVER, ComplexType x, CTQMC_TAG tag = I_ComplexType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Non-blocking sending of integer x to process "RECEIVER"
  @RECEIVER - the id of a process - recipient of a ComplexType
  @x - the ComplexType to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------*/
/*----------------------NON-BLOCKING RECEIVING ROUTINES------------------------------*/
/*-----------------------------------------------------------------------------------*/


MPI_Request& irecv_VectorType(const int SENDER, VectorType& result, int size, CTQMC_TAG tag = I_VectorType_SEND_TAG); 
MPI_Request& irecv_VectorType(const int SENDER, VectorType& result, CTQMC_TAG tag = I_VectorType_SEND_TAG);
MPI_Request& irecv(const int SENDER, VectorType& result, CTQMC_TAG tag = I_VectorType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a VectorType (size) from process "SENDER"
  @result - caught VectorType
  @SENDER - the id of a process which sends a VectorType
  @size - the size (number of cols) of a VectorType
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& irecv_RealVectorType(const int SENDER, RealVectorType& result, int size, CTQMC_TAG tag = I_RealVectorType_SEND_TAG); 
MPI_Request& irecv_RealVectorType(const int SENDER, RealVectorType& result, CTQMC_TAG tag = I_RealVectorType_SEND_TAG);
MPI_Request& irecv(const int SENDER, RealVectorType& result, CTQMC_TAG tag = I_RealVectorType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a RealVectorType (size) from process "SENDER"
  @result - caught RealVectorType
  @SENDER - the id of a process which sends a RealVectorType
  @size - the size (number of cols) of a RealVectorType
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& irecv_IntVectorType(const int SENDER, IntVectorType& result, int size, CTQMC_TAG tag = I_IntVectorType_SEND_TAG); 
MPI_Request& irecv_IntVectorType(const int SENDER, IntVectorType& result, CTQMC_TAG tag = I_IntVectorType_SEND_TAG);
MPI_Request& irecv(const int SENDER, IntVectorType& result, CTQMC_TAG tag = I_IntVectorType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Receive a IntVectorType (size) from process "SENDER"
  @result - caught IntVectorType
  @SENDER - the id of a process which sends a IntVectorType
  @size - the size (number of cols) of a IntVectorType
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& irecv_float(const int SENDER, n_type &x, CTQMC_TAG tag = I_FLOAT_SEND_TAG);
MPI_Request& irecv(const int SENDER, n_type &x, CTQMC_TAG tag = I_FLOAT_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Initialize and perform non-blocking receive of integer x from process "SENDER"
  @SENDER - the id of a process - recipient
  @x - the integer to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& irecv_complex(const int SENDER, ComplexType &x, CTQMC_TAG tag = I_ComplexType_SEND_TAG);
MPI_Request& irecv(const int SENDER, ComplexType &x, CTQMC_TAG tag = I_ComplexType_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Initialize and perform non-blocking receive of integer x from process "SENDER"
  @SENDER - the id of a process - recipient
  @x - the integer to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/

MPI_Request& irecv_int(const int SENDER, int &x, CTQMC_TAG tag = I_INT_SEND_TAG);
MPI_Request& irecv(const int SENDER, int &x, CTQMC_TAG tag = I_INT_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/**
  Initialize and perform non-blocking receive of integer x from process "SENDER"
  @SENDER - the id of a process - recipient
  @x - the integer to be sent
  Returns requests to be checked afterwards
 */
/*-----------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------*/
/*------------------------------RAW DATA ROUTINES------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void send_complexRaw1d (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG);
void send_complexRaw2d (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG);
void send_complexRaw3d (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG);
void send_complexRaw4d (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG);

void send (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG);
void send (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG);
void send (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG);
void send (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG);

MPI_Request* isend_complexRaw1d (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG);
MPI_Request* isend_complexRaw2d (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG);
MPI_Request* isend_complexRaw3d (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG);
MPI_Request* isend_complexRaw4d (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG);

MPI_Request* isend (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG);
MPI_Request* isend (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG);
MPI_Request* isend (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG);
MPI_Request* isend (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG);
/*-----------------------------------------------------------------------------------*/
/** 
  Non-blocking send a ComplexType RawData (size) to process "RECEIVER". 
  @x - raw data to be sent
  @RECEIVER - the id of a process which sends an object
  @size$x - the size of a dimension of an object
 */
/*-----------------------------------------------------------------------------------*/

void recv_complexRaw1d(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG); 
void recv_complexRaw2d(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG); 
void recv_complexRaw3d(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG); 
void recv_complexRaw4d(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG); 

void recv(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG); 
void recv(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG); 
void recv(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG); 
void recv(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG); 

MPI_Request* irecv_complexRaw1d(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG); 
MPI_Request* irecv_complexRaw2d(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG); 
MPI_Request* irecv_complexRaw3d(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG); 
MPI_Request* irecv_complexRaw4d(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG);

MPI_Request* irecv(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag = RAW_COMPLEX_1D_SEND_TAG); 
MPI_Request* irecv(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag = RAW_COMPLEX_2D_SEND_TAG); 
MPI_Request* irecv(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag = RAW_COMPLEX_3D_SEND_TAG); 
MPI_Request* irecv(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag = RAW_COMPLEX_4D_SEND_TAG); 
/*-----------------------------------------------------------------------------------*/
/** 
  Non-blocking receive a ComplexType RawData (size) from process "SENDER". 
  WARNING! DOESN'T ALLOCATE MEMORY!
  @result - received raw data
  @SENDER - the id of a process which sends an object
  @sizex - the size of a dimension of an object
 */
/*-----------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------*/
/*------------------------------COMMON ROUTINES--------------------------------------*/
/*-----------------------------------------------------------------------------------*/

void wait(MPI_Request &REQ);
void waitall(MPI_Request *REQ, int size);
int test(MPI_Request &REQ);
int testall(MPI_Request *REQ, int size);
int sync();
void summReduce(void *sendbuffer, void *receivebuffer, int size, int root_process);
void genericReduce(void *sendbuffer, void *receivebuffer, int size, MPI_Datatype datatype, MPI_Op Operation, int root_process);

/*-----------------------------------------------------------------------------------*/
/**
  Blocking and non-blocking test procedures for non-blocking sending/receiving raw data
  @REQ - requests to be checked
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


