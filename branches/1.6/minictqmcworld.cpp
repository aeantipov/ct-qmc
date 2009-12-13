/* This file is a part of 'ct-qmc 1.6' project
**
** minictqmcworld.cpp
** Implementations of the CTQMC_WORLD class methods.
**
** Author 			: Andrey Antipov, antipov@shg.ru
*/



#include "minictqmcworld.h"
#include <cstdlib>

/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::INIT(int argc, char **argv, int N)
{ int p;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  if (N>p) 
    {
      // Need to throw exception
      cout << "SHIT" << endl;
      exit(0);
    }

  MPI_Group MPI_GROUP_WORLD;
  MPI_Group CTQMC_PROCS; // THe group name for CT-QMC Processes

  int* process_ranks = new int[N];
  /* Make a list of the processes in the new communicator */
  for (int proc = 0; proc < N; proc++) process_ranks [proc] = proc;
  /* Get the group underlying MPI_COMM_WORLD */
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  /* Create the new group */
  MPI_Group_incl(MPI_GROUP_WORLD, N, process_ranks,&CTQMC_PROCS);
  /* Create the new communicator */
  MPI_Comm_create(MPI_COMM_WORLD, CTQMC_PROCS, &CTQMC_COMM);

  WORLD_SIZE=N;

  river_names=new char *[N];
  for (int i=0;i<N;i++) 
    { 
      river_names[i]=new char[RIVER_STREAM_SIZE]; 
      sprintf(river_names[i],"River%d.dat",i);
    }
}

/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::INIT(int argc, char **argv)
{ int N;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &N);
  MPI_Group MPI_GROUP_WORLD;
  MPI_Group CTQMC_PROCS; // THe group name for CT-QMC Processes

  int* process_ranks = new int[N];
  /* Make a list of the processes in the new communicator */
  for (int proc = 0; proc < N; proc++) process_ranks [proc] = proc;
  /* Get the group underlying MPI_COMM_WORLD */
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  /* Create the new group */
  MPI_Group_incl(MPI_GROUP_WORLD, N, process_ranks,&CTQMC_PROCS);
  /* Create the new communicator */
  MPI_Comm_create(MPI_COMM_WORLD, CTQMC_PROCS, &CTQMC_COMM);

  WORLD_SIZE=N;

  river_names=new char *[N];
  for (int i=0;i<N;i++) 
    { 
      river_names[i]=new char[RIVER_STREAM_SIZE]; 
      sprintf(river_names[i],"River%d.dat",i);
    }
}

/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::IniProcess()
{
  int q=(*this).getRank();
  rivers.open(river_names[q]);
  rivers << "This is output of process " << q << endl;
}

/*-----------------------------------------------------------------------------------*/

int CTQMC_WORLD::getRank()
{ int my_rank;
  MPI_Comm_rank(CTQMC_COMM, &my_rank);
  return my_rank;
}

/*-----------------------------------------------------------------------------------*/

int CTQMC_WORLD::getSize()
{ 
  return WORLD_SIZE;
}


/*-----------------------------------------------------------------------------------*/

ofstream & CTQMC_WORLD::getStream()
{
 return rivers;
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING VectorType--------------------------------*/
/*-----------------------------------------------------------------------------------*/
VectorType& CTQMC_WORLD::recv_VectorType(const int SENDER, int &size, CTQMC_TAG tag)
{
 static VectorType result;
 (*this).recv_VectorType(SENDER, result, size,tag);
 return result;
}

VectorType& CTQMC_WORLD::recv_VectorType(const int SENDER, CTQMC_TAG tag)
{
 static VectorType result;
 (*this).recv_VectorType(SENDER, result,tag);
 return result;
}

void CTQMC_WORLD::recv_VectorType(const int SENDER, VectorType& result, int &size, CTQMC_TAG tag ) 
{ 
  MPI_Recv(&size, 1, MPI::INT, SENDER, tag, CTQMC_COMM, &CT_STATUS);
  result.resize(size);
  MPI_Recv(&result[0], size, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &CT_STATUS);
}

void CTQMC_WORLD::recv_VectorType(const int SENDER, VectorType& result, CTQMC_TAG tag)
{ 
 int size;
 (*this).recv_VectorType(SENDER,result,size,tag);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING RealVectorType----------------------------*/
/*-----------------------------------------------------------------------------------*/
RealVectorType& CTQMC_WORLD::recv_RealVectorType(const int SENDER, int &size, CTQMC_TAG tag) 
{
 static RealVectorType result;
 (*this).recv_RealVectorType(SENDER, result, size, tag);
 return result;
}
RealVectorType& CTQMC_WORLD::recv_RealVectorType(const int SENDER, CTQMC_TAG tag)
{
 static RealVectorType result;
 (*this).recv_RealVectorType(SENDER, result, tag);
 return result;
}
void CTQMC_WORLD::recv_RealVectorType(const int SENDER, RealVectorType& result, int &size, CTQMC_TAG tag)
{
  MPI_Recv(&size, 1, MPI::INT, SENDER, tag, CTQMC_COMM, &CT_STATUS);
  result.resize(size);
  MPI_Recv(&result[0], size, MPI::DOUBLE, SENDER, tag, CTQMC_COMM, &CT_STATUS);
}
void CTQMC_WORLD::recv_RealVectorType(const int SENDER, RealVectorType& result, CTQMC_TAG tag)
{
 int size;
 (*this).recv_RealVectorType(SENDER, result,size, tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING IntVectorType-----------------------------*/
/*-----------------------------------------------------------------------------------*/
IntVectorType& CTQMC_WORLD::recv_IntVectorType(const int SENDER, int &size, CTQMC_TAG tag) 
{
 static IntVectorType result;
 (*this).recv_IntVectorType(SENDER, result, size, tag);
 return result;
}
IntVectorType& CTQMC_WORLD::recv_IntVectorType(const int SENDER, CTQMC_TAG tag)
{
 static IntVectorType result;
 (*this).recv_IntVectorType(SENDER, result, tag);
 return result;
}
void CTQMC_WORLD::recv_IntVectorType(const int SENDER, IntVectorType& result, int &size, CTQMC_TAG tag)
{
  MPI_Recv(&size, 1, MPI::INT, SENDER, tag, CTQMC_COMM, &CT_STATUS);
  result.resize(size);
  MPI_Recv(&result[0], size, MPI::INT, SENDER, tag, CTQMC_COMM, &CT_STATUS);
}

void CTQMC_WORLD::recv_IntVectorType(const int SENDER, IntVectorType& result, CTQMC_TAG tag)
{
 int size;
 (*this).recv_IntVectorType(SENDER, result,size, tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING INTEGER-----------------------------------*/
/*-----------------------------------------------------------------------------------*/
int CTQMC_WORLD::recv_int(const int SENDER, CTQMC_TAG tag) 
{
 int result;
 MPI_Recv(&result, 1, MPI::INT, SENDER, tag, CTQMC_COMM, &CT_STATUS);
 return result;
}
void CTQMC_WORLD::recv_int(const int SENDER, int& result, CTQMC_TAG tag)
{
  MPI_Recv(&result, 1, MPI::INT, SENDER, tag, CTQMC_COMM, &CT_STATUS);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING n_type----------------------------------*/
/*-----------------------------------------------------------------------------------*/
n_type CTQMC_WORLD::recv_float(const int SENDER, CTQMC_TAG tag)
{
  n_type result;
  MPI_Recv(&result, 1, MPI::DOUBLE, SENDER, tag, CTQMC_COMM, &CT_STATUS);
  return result;
}
void CTQMC_WORLD::recv_float(const int SENDER, n_type& result, CTQMC_TAG tag)
{
  MPI_Recv(&result, 1, MPI::DOUBLE, SENDER, tag, CTQMC_COMM, &CT_STATUS);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING ComplexType--------------------------------*/
/*-----------------------------------------------------------------------------------*/
ComplexType CTQMC_WORLD::recv_complex(const int SENDER, CTQMC_TAG tag)
{
  ComplexType result;
  MPI_Recv(&result, 1, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &CT_STATUS);
  return result;
}
void CTQMC_WORLD::recv_complex(const int SENDER, ComplexType& result, CTQMC_TAG tag)
{
  MPI_Recv(&result, 1, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &CT_STATUS);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING STD::STRING-------------------------------*/
/*-----------------------------------------------------------------------------------*/
string& CTQMC_WORLD::recv_string(const int SENDER, int &size, CTQMC_TAG tag) 
{
  static string str;
  (*this).recv_string(SENDER,str,size,tag);
  return str;
}
string& CTQMC_WORLD::recv_string(const int SENDER, CTQMC_TAG tag)
{
  static string str;
  (*this).recv_string(SENDER, str,tag);
  return str;
}
void CTQMC_WORLD::recv_string(const int SENDER, string& result, int &size, CTQMC_TAG tag)
{
  MPI_Recv(&size, 1, MPI::INT, SENDER, tag, CTQMC_COMM, &CT_STATUS);
  char *in = new char [size];
  MPI_Recv(in, size, MPI::CHAR, SENDER, tag, CTQMC_COMM, &CT_STATUS);
  result=in;
  delete in;

}
void CTQMC_WORLD::recv_string(const int SENDER, string& result, CTQMC_TAG tag)
{ 
 int size;
 (*this).recv_string(SENDER,result,size,tag);
}

/*-------------------------------COVERS FOR ALL RECV ROUTINES -----------------------*/

void CTQMC_WORLD::recv(const int SENDER, VectorType& result, CTQMC_TAG tag) {(*this).recv_VectorType(SENDER, result, tag);}
void CTQMC_WORLD::recv(const int SENDER, RealVectorType& result, CTQMC_TAG tag) {(*this).recv_RealVectorType(SENDER, result, tag);}
void CTQMC_WORLD::recv(const int SENDER, IntVectorType& result, CTQMC_TAG tag) {(*this).recv_IntVectorType(SENDER, result, tag);}
void CTQMC_WORLD::recv(const int SENDER, int& result, CTQMC_TAG tag) {(*this).recv_int(SENDER, result, tag);}
void CTQMC_WORLD::recv(const int SENDER, n_type& result, CTQMC_TAG tag) {(*this).recv_float(SENDER, result, tag);}
void CTQMC_WORLD::recv(const int SENDER, ComplexType& result, CTQMC_TAG tag) {(*this).recv_complex(SENDER, result, tag);}
void CTQMC_WORLD::recv(const int SENDER, string& result, CTQMC_TAG tag) {(*this).recv_string(SENDER, result, tag);}


/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING VectorType----------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_VectorType (const int RECEIVER, VectorType &x, int size, CTQMC_TAG tag)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM);
  MPI_Send(&x[0], size, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM);
}

void CTQMC_WORLD::send_VectorType (const int RECEIVER, VectorType &x, CTQMC_TAG tag)
{
 int size=x.size();
 (*this).send_VectorType(RECEIVER,x,size,tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING RealVectorType------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_RealVectorType (const int RECEIVER, RealVectorType &x, int size, CTQMC_TAG tag)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM);
  MPI_Send(&x[0], size, MPI::DOUBLE, RECEIVER, tag, CTQMC_COMM);
}
void CTQMC_WORLD::send_RealVectorType (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag)
{
 int size=x.size();
 (*this).send_RealVectorType(RECEIVER,x,size,tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING INT_VECTOR----------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_IntVectorType (const int RECEIVER, IntVectorType &x, int size, CTQMC_TAG tag)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM);
  MPI_Send(&x[0], size, MPI::INT, RECEIVER, tag, CTQMC_COMM);
}
void CTQMC_WORLD::send_IntVectorType (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag)
{
 int size=x.size();
 (*this).send_IntVectorType(RECEIVER,x,size,tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING INT-----------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_int(const int RECEIVER, int x, CTQMC_TAG tag)
{
  MPI_Send(&x, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING n_type--------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_float(const int RECEIVER, n_type x, CTQMC_TAG tag)
{
  MPI_Send(&x, 1 , MPI::DOUBLE, RECEIVER, tag, CTQMC_COMM);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING ComplexType---------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_complex(const int RECEIVER, ComplexType x, CTQMC_TAG tag)
{
  MPI_Send(&x, 1 , MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING STRING--------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_string(const int RECEIVER, string &x, int size, CTQMC_TAG tag)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM);
  char *buf = new char[size];
  strcpy(buf,x.c_str());
  MPI_Send(&buf[0], size, MPI::CHAR, RECEIVER, tag, CTQMC_COMM);
  delete buf;
}

void CTQMC_WORLD::send_string(const int RECEIVER, string &x, CTQMC_TAG tag)
{
 (*this).send_string(RECEIVER,x,x.size()+1,tag);
}

/*-------------------------------COVERS FOR ALL SEND ROUTINES -----------------------*/

void CTQMC_WORLD::send (const int RECEIVER, VectorType &x, CTQMC_TAG tag) {(*this).send_VectorType(RECEIVER,x,tag); }
void CTQMC_WORLD::send (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag) {(*this).send_RealVectorType(RECEIVER,x,tag); }
void CTQMC_WORLD::send (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag) {(*this).send_IntVectorType(RECEIVER,x,tag); }
void CTQMC_WORLD::send (const int RECEIVER, int x, CTQMC_TAG tag) {(*this).send_int(RECEIVER,x,tag); }
void CTQMC_WORLD::send (const int RECEIVER, n_type x, CTQMC_TAG tag) {(*this).send_float(RECEIVER,x,tag); }
void CTQMC_WORLD::send (const int RECEIVER, ComplexType x, CTQMC_TAG tag) {(*this).send_complex(RECEIVER,x,tag); }
void CTQMC_WORLD::send (const int RECEIVER, string &x, CTQMC_TAG tag) {(*this).send_string(RECEIVER,x,tag); }

/*-----------------------------------------------------------------------------------*/
/*-------------------------------NON-BLOCKING ROUTINES-------------------------------*/
/*-----------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------*/
/*-------------------------------ISENDING VectorType---------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::isend_VectorType (const int RECEIVER, VectorType &x, int size, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  MPI_Isend(&size, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM, &REQ);
  MPI_Wait(&REQ,&CT_STATUS); 
  MPI_Isend(&x[0], size, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend_VectorType (const int RECEIVER, VectorType &x, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  int size=x.size();
  REQ=(*this).isend_VectorType(RECEIVER,x,size,tag);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend (const int RECEIVER, VectorType &x, CTQMC_TAG tag)
{
  return (*this).isend_VectorType(RECEIVER,x,tag);
}


/*-----------------------------------------------------------------------------------*/
/*-------------------------------ISENDING RealVectorType-----------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::isend_RealVectorType (const int RECEIVER, RealVectorType &x, int size, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  MPI_Isend(&size, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM, &REQ);
  MPI_Wait(&REQ,&CT_STATUS); 
  MPI_Isend(&x[0], size, MPI::DOUBLE, RECEIVER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend_RealVectorType (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  int size=x.size();
  REQ=(*this).isend_RealVectorType(RECEIVER,x,size,tag);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend (const int RECEIVER, RealVectorType &x, CTQMC_TAG tag)
{
  return (*this).isend_RealVectorType(RECEIVER,x,tag);
}



/*-----------------------------------------------------------------------------------*/
/*-------------------------------ISENDING INT_VECTOR---------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::isend_IntVectorType (const int RECEIVER, IntVectorType &x, int size, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  MPI_Isend(&size, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM, &REQ);
  MPI_Wait(&REQ,&CT_STATUS); 
  MPI_Isend(&x[0], size, MPI::INT, RECEIVER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend_IntVectorType (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  int size=x.size();
  REQ=(*this).isend_IntVectorType(RECEIVER,x,size,tag);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend (const int RECEIVER, IntVectorType &x, CTQMC_TAG tag)
{
  return (*this).isend_IntVectorType(RECEIVER,x,tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------ISENDING INT----------------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::isend_int(const int RECEIVER, int x, CTQMC_TAG tag)
{ 
  static MPI_Request REQ;
  MPI_Isend(&x, 1 , MPI::INT, RECEIVER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend(const int RECEIVER, int x, CTQMC_TAG tag)
{
 return (*this).isend_int(RECEIVER,x,tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------ISENDING n_type-------------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::isend_float(const int RECEIVER, n_type x, CTQMC_TAG tag)
{ 
  static MPI_Request REQ;
  MPI_Isend(&x, 1 , MPI::DOUBLE, RECEIVER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend(const int RECEIVER, n_type x, CTQMC_TAG tag)
{
 return (*this).isend_float(RECEIVER,x,tag);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------ISENDING ComplexType--------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::isend_complex(const int RECEIVER, ComplexType x, CTQMC_TAG tag)
{ 
  static MPI_Request REQ;
  MPI_Isend(&x, 1 , MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend(const int RECEIVER, ComplexType x, CTQMC_TAG tag)
{
 return (*this).isend_complex(RECEIVER,x,tag);
}

/*-----------------------------------------------------------------------------------*/
/*--------------------------NON-BLOCKING RECEIVING ROUTINES--------------------------*/
/*-----------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------*/
/*-------------------------------IRECEIVING VectorType-------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::irecv_VectorType(const int SENDER, VectorType& result, int &size, CTQMC_TAG tag) 
{ 
  static MPI_Request REQ;
  MPI_Irecv(&size, 1 , MPI::INT, SENDER, tag, CTQMC_COMM, &REQ);
  MPI_Wait(&REQ,&CT_STATUS); 
  result.resize(size);
  MPI_Irecv(&result[0], size , MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::irecv_VectorType(const int SENDER, VectorType& result, CTQMC_TAG tag)
{ 
 int size;
 return (*this).irecv_VectorType(SENDER,result,size,tag);
}

MPI_Request& CTQMC_WORLD::irecv(const int SENDER, VectorType& result, CTQMC_TAG tag) {return (*this).irecv_VectorType(SENDER, result,tag); }

/*-----------------------------------------------------------------------------------*/
/*----------------------------IRECEIVING RealVectorType------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::irecv_RealVectorType(const int SENDER, RealVectorType& result, int &size, CTQMC_TAG tag) 
{ 
  static MPI_Request REQ;
  MPI_Irecv(&size, 1 , MPI::INT, SENDER, tag, CTQMC_COMM, &REQ);
  MPI_Wait(&REQ,&CT_STATUS);
  result.resize(size);
  MPI_Irecv(&result[0], size , MPI::DOUBLE, SENDER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::irecv_RealVectorType(const int SENDER, RealVectorType& result, CTQMC_TAG tag)
{ 
 int size;
 return (*this).irecv_RealVectorType(SENDER,result,size,tag);
}

MPI_Request& CTQMC_WORLD::irecv(const int SENDER, RealVectorType& result, CTQMC_TAG tag) {return (*this).irecv_RealVectorType(SENDER, result,tag); }

/*-----------------------------------------------------------------------------------*/
/*-----------------------------IRECEIVING IntVectorType------------------------------*/
/*-----------------------------------------------------------------------------------*/

MPI_Request& CTQMC_WORLD::irecv_IntVectorType(const int SENDER, IntVectorType& result, int &size, CTQMC_TAG tag) 
{ 
  static MPI_Request REQ;
  MPI_Irecv(&size, 1 , MPI::INT, SENDER, tag, CTQMC_COMM, &REQ);
  MPI_Wait(&REQ,&CT_STATUS); 
  result.resize(size);
  MPI_Irecv(&result[0], size , MPI::INT, SENDER, tag, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::irecv_IntVectorType(const int SENDER, IntVectorType& result, CTQMC_TAG tag)
{ 
 int size;
 return (*this).irecv_IntVectorType(SENDER,result,size,tag);
}

MPI_Request& CTQMC_WORLD::irecv(const int SENDER, IntVectorType& result, CTQMC_TAG tag) {return (*this).irecv_IntVectorType(SENDER, result); }


/*-----------------------------------------------------------------------------------*/
MPI_Request& CTQMC_WORLD::irecv_complex(const int SENDER, ComplexType &x, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  MPI_Irecv(&x, 1 , MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &REQ);
  return REQ;
}


MPI_Request& CTQMC_WORLD::irecv(const int SENDER, ComplexType &x, CTQMC_TAG tag) {return (*this).irecv_complex(SENDER, x,tag); }


/*-----------------------------------------------------------------------------------*/
MPI_Request& CTQMC_WORLD::irecv_float(const int SENDER, n_type &x, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  MPI_Irecv(&x, 1 , MPI::DOUBLE, SENDER, tag, CTQMC_COMM, &REQ);
  return REQ;
}


MPI_Request& CTQMC_WORLD::irecv(const int SENDER, n_type &x, CTQMC_TAG tag) {return (*this).irecv_float(SENDER, x,tag); }


/*-----------------------------------------------------------------------------------*/
MPI_Request& CTQMC_WORLD::irecv_int(const int SENDER, int &x, CTQMC_TAG tag)
{
  static MPI_Request REQ;
  MPI_Irecv(&x, 1 , MPI::INT, SENDER, tag, CTQMC_COMM, &REQ);
  return REQ;
}


MPI_Request& CTQMC_WORLD::irecv(const int SENDER, int &x, CTQMC_TAG tag) {return (*this).irecv_int(SENDER, x,tag); }

/*-----------------------------------------------------------------------------------*/
/*------------------------------RAW DATA ROUTINES------------------------------------*/
/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::send_complexRaw1d (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag)
{
  MPI_Send(x, size, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM);
}

void CTQMC_WORLD::send_complexRaw2d (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag)
{
  for (int i=0;i<size1;i++) 
   {
      MPI_Send(x[i], size2, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM);
   }
}

void CTQMC_WORLD::send_complexRaw3d (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag)
{
  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
   {
      MPI_Send(x[i][j], size3, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM);
   }
}

void CTQMC_WORLD::send_complexRaw4d (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag)
{
  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
      for (int k=0;k<size3;k++)
   {
      MPI_Send(x[i][j][k], size4, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM);
   }
}

void CTQMC_WORLD::send (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag) 
{(*this).send_complexRaw1d(RECEIVER,x,size,tag);};
void CTQMC_WORLD::send (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag)
{(*this).send_complexRaw2d(RECEIVER,x,size1,size2,tag);};
void CTQMC_WORLD::send (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag)
{(*this).send_complexRaw3d(RECEIVER,x,size1,size2,size3,tag);};
void CTQMC_WORLD::send (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag)
{(*this).send_complexRaw4d(RECEIVER,x,size1,size2,size3,size4,tag);};

/*-----------------------------------------------------------------------------------*/

MPI_Request* CTQMC_WORLD::isend_complexRaw1d (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag)
{
  MPI_Request * REQ = new MPI_Request;
  MPI_Isend(x, size, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM,REQ);
  return REQ;
}

MPI_Request* CTQMC_WORLD::isend_complexRaw2d (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag)
{
  MPI_Request *REQ = new MPI_Request [size1];
  for (int i=0;i<size1;i++) 
   {
      MPI_Isend(x[i], size2, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM, &REQ[i]);
   }
  return REQ;
}

MPI_Request* CTQMC_WORLD::isend_complexRaw3d (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag)
{
  MPI_Request *REQ = new MPI_Request [size1*size2];
  int count = 0;
  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
   {
      MPI_Isend(x[i][j], size3, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM, &REQ[count]);
      count++;
   }
  return REQ;
}

MPI_Request* CTQMC_WORLD::isend_complexRaw4d (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag)
{
  MPI_Request *REQ = new MPI_Request [size1*size2*size3];
  int count=0;
  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
      for (int k=0;k<size3;k++)
   {
      MPI_Isend(x[i][j][k], size4, MPI::DOUBLE_COMPLEX, RECEIVER, tag, CTQMC_COMM, &REQ[count]);
      count++;
   }
  return REQ;
}

MPI_Request* CTQMC_WORLD::isend (const int RECEIVER, ComplexType *x, int size, CTQMC_TAG tag) 
{return (*this).isend_complexRaw1d(RECEIVER,x,size,tag);};
MPI_Request* CTQMC_WORLD::isend (const int RECEIVER, ComplexType **x, int size1, int size2, CTQMC_TAG tag) 
{return (*this).isend_complexRaw2d(RECEIVER,x,size1,size2,tag);};
MPI_Request* CTQMC_WORLD::isend (const int RECEIVER, ComplexType ***x, int size1, int size2, int size3, CTQMC_TAG tag) 
{return (*this).isend_complexRaw3d(RECEIVER,x,size1,size2,size3,tag);};
MPI_Request* CTQMC_WORLD::isend (const int RECEIVER, ComplexType ****x, int size1, int size2, int size3, int size4, CTQMC_TAG tag) 
{return (*this).isend_complexRaw4d(RECEIVER,x,size1,size2,size3,size4,tag);};


/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::recv_complexRaw1d(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag)
{
  MPI_Recv(result, size, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &CT_STATUS);
}

void CTQMC_WORLD::recv_complexRaw2d(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag)
{

  for (int i=0;i<size1;i++)
    {
       MPI_Recv(result[i], size2, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &CT_STATUS);
    }

}

void CTQMC_WORLD::recv_complexRaw3d(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag)
{
  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
    {
       MPI_Recv(result[i][j], size3, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &CT_STATUS);
    }

}

void CTQMC_WORLD::recv_complexRaw4d(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag)
{
  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
      for (int k=0;k<size3;k++)
    {
       MPI_Recv(result[i][j][k], size4, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &CT_STATUS);
    }

}

void CTQMC_WORLD::recv(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag)
{(*this).recv_complexRaw1d(SENDER,result,size,tag);}
void CTQMC_WORLD::recv(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag)
{(*this).recv_complexRaw2d(SENDER,result,size1,size2,tag);}
void CTQMC_WORLD::recv(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag)
{(*this).recv_complexRaw3d(SENDER,result,size1,size2,size3,tag);}
void CTQMC_WORLD::recv(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag)
{(*this).recv_complexRaw4d(SENDER,result,size1,size2,size3,size4,tag);}
/*-----------------------------------------------------------------------------------*/

MPI_Request* CTQMC_WORLD::irecv_complexRaw1d(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag)
{
  MPI_Request *REQ = new MPI_Request;
  MPI_Irecv(result, size, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, REQ);
  return REQ;
}
MPI_Request* CTQMC_WORLD::irecv_complexRaw2d(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag)
{
  MPI_Request *REQ = new MPI_Request [size1];
  int count=0;
  for (int i=0;i<size1;i++)
    {
       MPI_Irecv(result[i], size2, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &REQ[count]);
       count++;
    }
  return REQ;
}


MPI_Request* CTQMC_WORLD::irecv_complexRaw3d(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag)
{
  MPI_Request *REQ = new MPI_Request [size1*size2];
  int count=0;

  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
    {
       MPI_Irecv(result[i][j], size3, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &REQ[count]);
       count++;
    }
  return REQ;
}

MPI_Request* CTQMC_WORLD::irecv_complexRaw4d(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag)
{
  MPI_Request *REQ = new MPI_Request [size1*size2*size3];
  int count=0;
  for (int i=0;i<size1;i++)
    for (int j=0;j<size2;j++)
      for (int k=0;k<size3;k++)
    {
       MPI_Irecv(result[i][j][k], size4, MPI::DOUBLE_COMPLEX, SENDER, tag, CTQMC_COMM, &REQ[count]);
       count++;
    }
  return REQ;
}

MPI_Request* CTQMC_WORLD::irecv(const int SENDER, ComplexType* result, int size, CTQMC_TAG tag)
{return (*this).irecv_complexRaw1d(SENDER,result,size,tag);}
MPI_Request* CTQMC_WORLD::irecv(const int SENDER, ComplexType** result, int size1, int size2, CTQMC_TAG tag)
{return (*this).irecv_complexRaw2d(SENDER,result,size1,size2,tag);}
MPI_Request* CTQMC_WORLD::irecv(const int SENDER, ComplexType*** result, int size1, int size2, int size3, CTQMC_TAG tag)
{return (*this).irecv_complexRaw3d(SENDER,result,size1,size2,size3,tag);}
MPI_Request* CTQMC_WORLD::irecv(const int SENDER, ComplexType**** result, int size1, int size2, int size3, int size4, CTQMC_TAG tag)
{return (*this).irecv_complexRaw4d(SENDER,result,size1,size2,size3,size4,tag);}

/*-----------------------------------------------------------------------------------*/
/*--------------------------------COMMON PROCEDURES----------------------------------*/
/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::wait (MPI_Request &REQ)
{
  MPI_Wait(&REQ,&CT_STATUS); 
}

int CTQMC_WORLD::test (MPI_Request &REQ)
{
  int flag;
  MPI_Test(&REQ,&flag,&CT_STATUS);
  return flag;
}

int CTQMC_WORLD::testall (MPI_Request *REQ, int size)
{
  int flag;
  MPI_Status *temp = new MPI_Status[size];
  MPI_Testall(size,REQ,&flag,temp);
  delete[] temp;
  return flag;
}

int CTQMC_WORLD::sync ()
{
  return MPI_Barrier(CTQMC_COMM); 
}


/*-----------------------------------------------------------------------------------*/
CTQMC_WORLD::~CTQMC_WORLD()
{
  for (int i=0;i<WORLD_SIZE;i++)
    {
//      rivers[i].close();
//      delete river_names[i]; 
    }
  MPI_Finalize();
}


/*--------------------------END OF CLASS CTQMC_WORLD---------------------------------*/



