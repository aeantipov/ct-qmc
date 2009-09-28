typedef double n_type; //accuracy
int const n_zone=2, n_part=1;
#include "headers.cpp"


void main_(int argc, char **argv)
{
   char * text; int L;
   {
   	ifstream inp("_input.dat",ios::in|ios::binary|ios::ate);
      if (!inp.is_open()) {ofstream ou("_input.dat",ios::binary); ou<<"Added by \"change\"\n$"<<argv[argc-2]<<"$\n"<<argv[argc-1]<<"\n\n"; return;}
   	L=inp.tellg();inp.seekg (0, ios::beg);
    	text=new char[L+256];inp.read (text, L);
   ;}

   int length_name=0; for (;argv[argc-2][length_name]!=0;length_name++);


   	bool found=false; int i1, i2;
   	for (int i=0; i<L-1; i++)
      if (text[i]=='$')
      {int j;
      	for (j=i+1; j<L; j++)
         if (text[j]=='$')
         {
         	found=(length_name==j-i-1);
         	for (int k=i+1; k<j; k++) if (text[k]!=argv[argc-2][k-i-1]) found=false;
         	break;
         ;}
         if (found)
         {
         	i1=j;  for(;text[i1]!='\n';i1++);
            i2=i1+1; for(;text[i2]!='\n';i2++);
         	break;
         ;}
      ;}

//      cout<<i1<<"  "<<i2;


      ofstream ou("_input.dat",ios::binary);
//      if (!found) {ou<<text<<"\n$"<<argv[argc-2]<<"$\n"<<argv[argc-1]<<"\nAdded by \"change\"\n"; return;}
      if (!found) {for (int j=0; j<L; j++) ou<<text[j]; ou<<"\n$"<<argv[argc-2]<<"$\n"<<argv[argc-1]<<"\nAdded by \"change\"\n"; return;}
      {for (int j=0; j<=i1; j++) ou<<text[j];}
      ou<<argv[argc-1];
      {for (int j=i2; j<L; j++) ou<<text[j];}

	return;
;}


int main (int argc, char *argv[])
{
	for (;argc>1;argc-=2) main_(argc, argv);
   return 1;
;}
