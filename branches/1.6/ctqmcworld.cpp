/**
   IMPLEMENTATIONS OF CTQMCWORLD.H ROUTINES
 */


#include "ctqmcworld.h"

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
/*-------------------------------RECEIVING MatrixType-----------------------------------*/
/*-----------------------------------------------------------------------------------*/

MatrixType& CTQMC_WORLD::catch_MatrixType(const int SENDER, int &sizex, int &sizey) 
// Receive a matrix (sizex * sizey) from process "from"
{
  static MatrixType x;
  // Catch dimensions
  MPI_Recv(&sizex, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  MPI_Recv(&sizey, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Catch data as a buffer
  //  RealType in[sizex*sizey];
  ScalarType *in = new ScalarType [sizex*sizey];
  MPI_Recv(in, sizex*sizey, MPI::DOUBLE_COMPLEX, SENDER, MatrixType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  x=Eigen::Map<MatrixType> (in,sizex,sizey).transpose();
  delete in;
  return x;
}

/*-----------------------------------------------------------------------------------*/

MatrixType& CTQMC_WORLD::catch_MatrixType(const int SENDER) 
{int sizex,sizey;
 return catch_MatrixType(SENDER,sizex,sizey);
}

/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::catch_MatrixType(MatrixType& result, const int SENDER, int &sizex, int &sizey) // Receive a matrix (sizex * sizey) from process "from"
{
  // Catch dimensions
  MPI_Recv(&sizex, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  MPI_Recv(&sizey, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Catch data as a buffer
  //  RealType in[sizex*sizey];
  ScalarType *in = new ScalarType [sizex*sizey];
  MPI_Recv(in, sizex*sizey, MPI::DOUBLE_COMPLEX, SENDER, MatrixType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  result=Eigen::Map<MatrixType> (in,sizex,sizey).transpose();
  delete in;
}

/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::catch_MatrixType(MatrixType& result, const int SENDER)
{int sizex,sizey;
 catch_MatrixType(result,SENDER,sizex,sizey);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING RealMatrixType-----------------------------------*/
/*-----------------------------------------------------------------------------------*/

RealMatrixType& CTQMC_WORLD::catch_RealMatrixType(const int SENDER, int &sizex, int &sizey) 
// Receive a matrix (sizex * sizey) from process "from"
{
  static RealMatrixType x;
  // Catch dimensions
  MPI_Recv(&sizex, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  MPI_Recv(&sizey, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Catch data as a buffer
  //  RealType in[sizex*sizey];
  RealType *in = new RealType [sizex*sizey];
  MPI_Recv(in, sizex*sizey, MPI::DOUBLE, SENDER, RealMatrixType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  x=Eigen::Map<RealMatrixType> (in,sizex,sizey).transpose();
  delete in;
  return x;
}

/*-----------------------------------------------------------------------------------*/

RealMatrixType& CTQMC_WORLD::catch_RealMatrixType(const int SENDER) 
{int sizex,sizey;
 return catch_RealMatrixType(SENDER,sizex,sizey);
}

/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::catch_RealMatrixType(RealMatrixType& result, const int SENDER, int &sizex, int &sizey) // Receive a matrix (sizex * sizey) from process "from"
{
  // Catch dimensions
  MPI_Recv(&sizex, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  MPI_Recv(&sizey, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Catch data as a buffer
  //  RealType in[sizex*sizey];
  RealType *in = new RealType [sizex*sizey];
  MPI_Recv(in, sizex*sizey, MPI::DOUBLE, SENDER, RealMatrixType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  result=Eigen::Map<RealMatrixType> (in,sizex,sizey).transpose();
  delete in;
}

/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::catch_RealMatrixType(RealMatrixType& result, const int SENDER)
{int sizex,sizey;
 catch_RealMatrixType(result,SENDER,sizex,sizey);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING VectorType-----------------------------------*/
/*-----------------------------------------------------------------------------------*/
VectorType& CTQMC_WORLD::catch_VectorType(const int SENDER, int &size)
{
 static VectorType result;
 (*this).catch_VectorType(result, SENDER, size);
 return result;
}

VectorType& CTQMC_WORLD::catch_VectorType(const int SENDER)
{
 static VectorType result;
 (*this).catch_VectorType(result, SENDER);
 return result;
}

void CTQMC_WORLD::catch_VectorType(VectorType& result, const int SENDER, int &size) 
{ 
  MPI_Recv(&size, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  ScalarType *in = new ScalarType [size];
  MPI_Recv(in, size, MPI::DOUBLE_COMPLEX, SENDER, VectorType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  result=Eigen::Map<VectorType> (in,size);
  delete in;

}

void CTQMC_WORLD::catch_VectorType(VectorType& result, const int SENDER)
{ 
 int size;
 (*this).catch_VectorType(result,SENDER,size);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING RealVectorType-----------------------------------*/
/*-----------------------------------------------------------------------------------*/
RealVectorType& CTQMC_WORLD::catch_RealVectorType(const int SENDER, int &size) 
{
 static RealVectorType result;
 (*this).catch_RealVectorType(result, SENDER, size);
 return result;
}
RealVectorType& CTQMC_WORLD::catch_RealVectorType(const int SENDER)
{
 static RealVectorType result;
 (*this).catch_RealVectorType(result, SENDER);
 return result;
}
void CTQMC_WORLD::catch_RealVectorType(RealVectorType& result, const int SENDER, int &size)
{
  MPI_Recv(&size, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  RealType *in = new RealType [size];
  MPI_Recv(in, size, MPI::DOUBLE, SENDER, RealVectorType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  result=Eigen::Map<RealVectorType> (in,size);
  delete in;

}
void CTQMC_WORLD::catch_RealVectorType(RealVectorType& result, const int SENDER)
{
 int size;
 (*this).catch_RealVectorType(result,SENDER,size);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING IntVectorType-----------------------------------*/
/*-----------------------------------------------------------------------------------*/
IntVectorType& CTQMC_WORLD::catch_IntVectorType(const int SENDER, int &size) 
{
 static IntVectorType result;
 (*this).catch_IntVectorType(result, SENDER, size);
 return result;
}
IntVectorType& CTQMC_WORLD::catch_IntVectorType(const int SENDER)
{
 static IntVectorType result;
 (*this).catch_IntVectorType(result, SENDER);
 return result;
}
void CTQMC_WORLD::catch_IntVectorType(IntVectorType& result, const int SENDER, int &size)
{
  MPI_Recv(&size, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  int *in = new int [size];
  MPI_Recv(in, size, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  result=Eigen::Map<IntVectorType> (in,size);
  delete in;

}
void CTQMC_WORLD::catch_IntVectorType(IntVectorType& result, const int SENDER)
{
 int size;
 (*this).catch_IntVectorType(result,SENDER,size);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING INTEGER-----------------------------------*/
/*-----------------------------------------------------------------------------------*/
int CTQMC_WORLD::catch_int(const int SENDER) 
{
 int result;
 MPI_Recv(&result, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
 return result;
}
void CTQMC_WORLD::catch_int(int& result, const int SENDER)
{
  MPI_Recv(&result, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING RealType------------------------------------*/
/*-----------------------------------------------------------------------------------*/
RealType CTQMC_WORLD::catch_float(const int SENDER)
{
  RealType result;
  MPI_Recv(&result, 1, MPI::DOUBLE, SENDER, RealType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  return result;
}
void CTQMC_WORLD::catch_float(RealType& result, const int SENDER)
{
  MPI_Recv(&result, 1, MPI::DOUBLE, SENDER, RealType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING ScalarType------------------------------------*/
/*-----------------------------------------------------------------------------------*/
ScalarType CTQMC_WORLD::catch_complex(const int SENDER)
{
  ScalarType result;
  MPI_Recv(&result, 1, MPI::DOUBLE_COMPLEX, SENDER, ScalarType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  return result;
}
void CTQMC_WORLD::catch_complex(ScalarType& result, const int SENDER)
{
  MPI_Recv(&result, 1, MPI::DOUBLE_COMPLEX, SENDER, ScalarType_SEND_TAG, CTQMC_COMM, &CT_STATUS);
}
/*-----------------------------------------------------------------------------------*/
/*-------------------------------RECEIVING STD::STRING-------------------------------*/
/*-----------------------------------------------------------------------------------*/
string& CTQMC_WORLD::catch_string(const int SENDER, int &size) 
{
  static string str;
  (*this).catch_string(str,SENDER,size);
  return str;
}
string& CTQMC_WORLD::catch_string(const int SENDER)
{
  static string str;
  (*this).catch_string(str,SENDER);
  return str;
}
void CTQMC_WORLD::catch_string(string& result, const int SENDER, int &size)
{
  MPI_Recv(&size, 1, MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  char *in = new char [size];
  MPI_Recv(in, size, MPI::CHAR, SENDER, STRING_SEND_TAG, CTQMC_COMM, &CT_STATUS);
  // Map to MatrixType_type
  result=in;
  delete in;

}
void CTQMC_WORLD::catch_string(string& result, const int SENDER)
{ 
 int size;
 (*this).catch_string(result,SENDER,size);
}

/*-------------------------------COVERS FOR ALL SEND ROUTINES -----------------------*/

void CTQMC_WORLD::recv(MatrixType& result, const int SENDER) {(*this).catch_MatrixType(result,SENDER);} 
void CTQMC_WORLD::recv(RealMatrixType& result, const int SENDER) {(*this).catch_RealMatrixType(result,SENDER);}
void CTQMC_WORLD::recv(VectorType& result, const int SENDER) {(*this).catch_VectorType(result,SENDER);}
void CTQMC_WORLD::recv(RealVectorType& result, const int SENDER) {(*this).catch_RealVectorType(result,SENDER);}
void CTQMC_WORLD::recv(int& result, const int SENDER) {(*this).catch_int(result,SENDER);}
void CTQMC_WORLD::recv(RealType& result, const int SENDER) {(*this).catch_float(result,SENDER);}
void CTQMC_WORLD::recv(ScalarType& result, const int SENDER) {(*this).catch_complex(result,SENDER);}
void CTQMC_WORLD::recv(string& result, const int SENDER) {(*this).catch_string(result,SENDER);}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING MatrixType-------------------------------------*/
/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::send_MatrixType (const int RECEIVER, MatrixType &x, int sizex, int sizey) 
{
  MPI_Send(&sizex, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  MPI_Send(&sizey, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  vector<ScalarType> buf (sizex*sizey);
  for (int i=0;i<sizey;i++) for (int j=0;j<sizex;j++) buf[i*sizex+j]=x(i,j);
  MPI_Send(&buf[0], sizex*sizey, MPI::DOUBLE_COMPLEX, RECEIVER, MatrixType_SEND_TAG, CTQMC_COMM);
  buf.resize(0);
}

/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::send_MatrixType (const int RECEIVER, MatrixType &x) 
{
  int sizex=x.rows();
  int sizey=x.cols();
  send_MatrixType(RECEIVER,x,sizex,sizey);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING RealMatrixType-------------------------------------*/
/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::send_RealMatrixType (const int RECEIVER, RealMatrixType &x, int sizex, int sizey) 
{
  MPI_Send(&sizex, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  MPI_Send(&sizey, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  vector<RealType> buf (sizex*sizey);
  for (int i=0;i<sizey;i++) for (int j=0;j<sizex;j++) buf[i*sizex+j]=x(i,j);
  MPI_Send(&buf[0], sizex*sizey, MPI::DOUBLE, RECEIVER, RealMatrixType_SEND_TAG, CTQMC_COMM);
  buf.resize(0);
}

/*-----------------------------------------------------------------------------------*/

void CTQMC_WORLD::send_RealMatrixType (const int RECEIVER, RealMatrixType &x) 
{
  int sizex=x.rows();
  int sizey=x.cols();
  send_RealMatrixType(RECEIVER,x,sizex,sizey);
}

/*-----------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING VectorType-------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_VectorType (const int RECEIVER, VectorType &x, int size)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  vector<ScalarType> buf (size);
  for (int i=0;i<size;i++) buf[i]=x(i);
  MPI_Send(&buf[0], size, MPI::DOUBLE_COMPLEX, RECEIVER, VectorType_SEND_TAG, CTQMC_COMM);
  buf.resize(0);

}

void CTQMC_WORLD::send_VectorType (const int RECEIVER, VectorType &x)
{
 int size=x.size();
 (*this).send_VectorType(RECEIVER,x,size);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING RealVectorType-------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_RealVectorType (const int RECEIVER, RealVectorType &x, int size)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  vector<RealType> buf (size);
  for (int i=0;i<size;i++) buf[i]=x(i);
  MPI_Send(&buf[0], size, MPI::DOUBLE, RECEIVER, RealVectorType_SEND_TAG, CTQMC_COMM);
  buf.resize(0);
}
void CTQMC_WORLD::send_RealVectorType (const int RECEIVER, RealVectorType &x)
{
 int size=x.size();
 (*this).send_RealVectorType(RECEIVER,x,size);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING INT_VECTOR----------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_IntVectorType (const int RECEIVER, IntVectorType &x, int size)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  vector<int> buf (size);
  for (int i=0;i<size;i++) buf[i]=x(i);
  MPI_Send(&buf[0], size, MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  buf.resize(0);
}
void CTQMC_WORLD::send_IntVectorType (const int RECEIVER, IntVectorType &x)
{
 int size=x.size();
 (*this).send_IntVectorType(RECEIVER,x,size);
}


/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING INT-----------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_int(const int RECEIVER, int x)
{
  MPI_Send(&x, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING RealType--------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_float(const int RECEIVER, RealType x)
{
  MPI_Send(&x, 1 , MPI::DOUBLE, RECEIVER, RealType_SEND_TAG, CTQMC_COMM);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING ScalarType--------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_complex(const int RECEIVER, ScalarType x)
{
  MPI_Send(&x, 1 , MPI::DOUBLE_COMPLEX, RECEIVER, ScalarType_SEND_TAG, CTQMC_COMM);
}

/*-----------------------------------------------------------------------------------*/
/*-------------------------------SENDING STRING--------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void CTQMC_WORLD::send_string(const int RECEIVER, string &x, int size)
{
  MPI_Send(&size, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM);
  char *buf = new char[size];
  strcpy(buf,x.c_str());
  MPI_Send(&buf[0], size, MPI::CHAR, RECEIVER, STRING_SEND_TAG, CTQMC_COMM);
  delete buf;
}

void CTQMC_WORLD::send_string(const int RECEIVER, string &x)
{
 (*this).send_string(RECEIVER,x,x.size()+1);
}

/*-------------------------------COVERS FOR ALL SEND ROUTINES -----------------------*/

void CTQMC_WORLD::send (const int RECEIVER, MatrixType &x) {(*this).send_MatrixType(RECEIVER,x); }
void CTQMC_WORLD::send (const int RECEIVER, RealMatrixType &x) {(*this).send_RealMatrixType(RECEIVER,x); } 
void CTQMC_WORLD::send (const int RECEIVER, VectorType &x) {(*this).send_VectorType(RECEIVER,x); }
void CTQMC_WORLD::send (const int RECEIVER, RealVectorType &x) {(*this).send_RealVectorType(RECEIVER,x); }
void CTQMC_WORLD::send (const int RECEIVER, IntVectorType &x) {(*this).send_IntVectorType(RECEIVER,x); }
void CTQMC_WORLD::send (const int RECEIVER, int x) {(*this).send_int(RECEIVER,x); }
void CTQMC_WORLD::send (const int RECEIVER, RealType x) {(*this).send_float(RECEIVER,x); }
void CTQMC_WORLD::send (const int RECEIVER, ScalarType x) {(*this).send_complex(RECEIVER,x); }
void CTQMC_WORLD::send (const int RECEIVER, string &x) {(*this).send_string(RECEIVER,x); }

/*-----------------------------------------------------------------------------------*/
/*-------------------------------NON-BLOCKING ROUTINES-------------------------------*/
/*-----------------------------------------------------------------------------------*/
MPI_Request& CTQMC_WORLD::isend_int(const int RECEIVER, int x)
{ 
  static MPI_Request REQ;
  MPI_Isend(&x, 1 , MPI::INT, RECEIVER, I_INT_SEND_TAG, CTQMC_COMM, &REQ);
  return REQ;
}

MPI_Request& CTQMC_WORLD::isend(const int RECEIVER, int x)
{
 return (*this).isend_int(RECEIVER,x);
}

/*-----------------------------------------------------------------------------------*/
MPI_Request& CTQMC_WORLD::irecv_int(int &x, const int SENDER)
{
  static MPI_Request REQ;
  MPI_Irecv(&x, 1 , MPI::INT, SENDER, I_INT_SEND_TAG, CTQMC_COMM, &REQ);
  return REQ;
}


MPI_Request& CTQMC_WORLD::irecv(int &x, const int SENDER) {return (*this).irecv_int(x,SENDER); }

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


/*-----------------------------------------------------------------------------------*/
/*-------------------------------RAW DATA ROUTINES-----------------------------------*/
/*-----------------------------------------------------------------------------------*/
/*void CTQMC_WORLD::isend_raw(const int RECEIVER, vector <ScalarType> &x, int size)
{
  static MPI_Request REQ[2];
  MPI_Isend(&size, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM, &REQ[0]);
  MPI_Isend(&x[0], size, MPI::DOUBLE_COMPLEX, RECEIVER, VectorType_SEND_TAG, CTQMC_COMM, &REQ[1]);
}

void CTQMC_WORLD::isend_raw(const int RECEIVER, vector <RealType> &x, int size)
{
  static MPI_Request REQ[2];
  MPI_Isend(&size, 1 , MPI::INT, RECEIVER, INT_SEND_TAG, CTQMC_COMM, &REQ[0]);
  MPI_Isend(&x[0], size, MPI::DOUBLE, RECEIVER, RealVectorType_SEND_TAG, CTQMC_COMM, &REQ[1]);
}*/
/*-----------------------------------------------------------------------------------*/
/*void CTQMC_WORLD::irecv_raw(const int SENDER, vector <ScalarType> &x, int &size)
{
  static MPI_Request REQ[2];
  MPI_Irecv(&size, 1 , MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &REQ[0]);
  MPI_Irecv(&x[0], size, MPI::DOUBLE_COMPLEX, SENDER, VectorType_SEND_TAG, CTQMC_COMM, &REQ[1]);
}

void CTQMC_WORLD::irecv_raw(const int SENDER, vector <RealType> &x, int &size)
{
  static MPI_Request REQ[2];
  MPI_Irecv(&size, 1 , MPI::INT, SENDER, INT_SEND_TAG, CTQMC_COMM, &REQ[0]);
  MPI_Irecv(&x[0], size, MPI::DOUBLE, SENDER, RealVectorType_SEND_TAG, CTQMC_COMM, &REQ[1]);
}*/

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



