typedef double n_type; //accuracy
int const n_zone=2, n_part=1;
#include "headers.cpp"


void main_(int argc, char **argv)
{
	cout<<argv[argc]<<"=";
	n_type x=value(argv[argc]);
	cout<<x<<"   ";
;}


int main (int argc, char *argv[])
{
	for (int a=1;a<argc;a++) main_(a, argv);
   	return 1;
;}
