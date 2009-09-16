#include<iostream>
#include<fstream>
#include<complex>
#include <time.h>

using namespace std;
#define complex complex<n_type>


//===================file input ======================

char input_values [255]="_input.dat", default_input [255]="/home/antipov/ct-qmc-1.5/input.dat";


int MAX_NUMBER_OF_VARS=2000, MAX_NUMBER_OF_COMMENT_CHARS=10000;

n_type value_(char* name, int index=0)
{

   char x [255];
	ifstream inf(default_input);
   int f, count=0;
   do
   {
   	int g=0; do {inf>>x;g++;} while (x[0]!='$' && g<MAX_NUMBER_OF_COMMENT_CHARS);

      {for (int i=0; x[i]!=0 && i<1000; i++) {x[i]=x[i+1]; if (x[i]=='$') {x[i]=0; x[i+1]=0;};};}
      f=1; int i=-1;
      do
      	{i++; if (x[i]!=name[i]) f=0;}
      while (x[i]!=0 && name[i]!=0 && i<255 && f==1);
      count++;
      if (g==MAX_NUMBER_OF_COMMENT_CHARS) count=MAX_NUMBER_OF_VARS; //EOF
   ;}
   while (f==0 && count<MAX_NUMBER_OF_VARS);

   if (count==MAX_NUMBER_OF_VARS)
   {
   	cout<<"NO INPUT VARIABLE FOR THE VARIABLE \""<<name<<"\"\n";
   	return 0;
   ;}
   n_type r;  for(int i=0; i<index+1; i++) inf>>r; return r;
;}

n_type value(char* name, int index=0)
{

	static int numb=0; numb++; if (numb>MAX_NUMBER_OF_VARS) {cout<<"VALUE is called too many times. Check the code!\n"; numb=0;}

   char x [255];
	ifstream inf(input_values);
   int f, count=0;
   do
   {
   	int g=0; do {inf>>x;g++;} while (x[0]!='$' && g<MAX_NUMBER_OF_COMMENT_CHARS);

      {for (int i=0; x[i]!=0 && i<1000; i++) {x[i]=x[i+1]; if (x[i]=='$') {x[i]=0; x[i+1]=0;};};}
      f=1; int i=-1;
      do
      	{i++; if (x[i]!=name[i]) f=0;}
      while (x[i]!=0 && name[i]!=0 && i<255 && f==1);
      count++;
      if (g==MAX_NUMBER_OF_COMMENT_CHARS) count=MAX_NUMBER_OF_VARS; //EOF
   ;}
   while (f==0 && count<MAX_NUMBER_OF_VARS);

   if (count==MAX_NUMBER_OF_VARS)   	return value_(name,index);
   n_type r;  for(int i=0; i<index+1; i++) inf>>r; return r;

;}



int int_value_(char* name, int index=0)
{

   char x [255];
	ifstream inf(default_input);
   int f, count=0;
   do
   {
   	int g=0; do {inf>>x;g++;} while (x[0]!='$' && g<MAX_NUMBER_OF_COMMENT_CHARS);
      {for (int i=0; x[i]!=0 && i<1000; i++) {x[i]=x[i+1]; if (x[i]=='$') {x[i]=0; x[i+1]=0;};};}
      f=1; int i=-1;
      do
      	{i++; if (x[i]!=name[i]) f=0;}
      while (x[i]!=0 && name[i]!=0 && i<255 && f==1);
      count++;
      if (g==MAX_NUMBER_OF_COMMENT_CHARS) count=MAX_NUMBER_OF_VARS; //EOF
   ;}
   while (f==0 && count<MAX_NUMBER_OF_VARS);

   if (count==MAX_NUMBER_OF_VARS)
   {cout<<"NO INPUT VARIABLE FOR THE VARIABLE \""<<name<<"\"\n"; return 0;}
   int r;  for(int i=0; i<index+1; i++) inf>>r; return r;

;}


int int_value(char* name, int index=0)
{

	static int numb=0; numb++; if (numb>MAX_NUMBER_OF_VARS) {cout<<"INT_VALUE is called too many times. Check the code!\n"; numb=0;}

   char x [255];
	ifstream inf(input_values);
   int f, count=0;
   do
   {
   	int g=0; do {inf>>x;g++;} while (x[0]!='$' && g<MAX_NUMBER_OF_COMMENT_CHARS);
      {for (int i=0; x[i]!=0 && i<1000; i++) {x[i]=x[i+1]; if (x[i]=='$') {x[i]=0; x[i+1]=0;};};}
      f=1; int i=-1;
      do
      	{i++; if (x[i]!=name[i]) f=0;}
      while (x[i]!=0 && name[i]!=0 && i<255 && f==1);
      count++;
      if (g==MAX_NUMBER_OF_COMMENT_CHARS) count=MAX_NUMBER_OF_VARS; //EOF
   ;}
   while (f==0 && count<MAX_NUMBER_OF_VARS);

   if (count==MAX_NUMBER_OF_VARS) return int_value_(name,index);
   int r=3; for(int i=0; i<index+1; i++) inf>>r; return r;

;}



void change_parameter(char *argv, n_type val)
{
   char * text; int L;
   {
   	ifstream inp("_input.dat",ios::in|ios::binary|ios::ate);
      if (!inp.is_open()) {ofstream ou("_input.dat",ios::binary); ou<<"Added by \"change\"\n$"<<argv<<"$\n"<<val<<"\n\n"; return;}
   	L=inp.tellg();inp.seekg (0, ios::beg);
    	text=new char[L+256];inp.read (text, L);
   ;}

   int length_name=0; for (;argv[length_name]!=0;length_name++);


   	bool found=false; int i1, i2;
   	for (int i=0; i<L-1; i++)
      if (text[i]=='$')
      {
      	int j; for (j=i+1; j<L; j++)
         if (text[j]=='$')
         {
         	found=(length_name==j-i-1);
         	for (int k=i+1; k<j; k++) if (text[k]!=argv[k-i-1]) found=false;
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
      if (!found) {ou<<text<<"\n$"<<argv<<"$\n"<<val<<"\nAdded by \"change\"\n"; return;}
      {for (int j=0; j<=i1; j++) ou<<text[j];}
      ou<<val;
      {for (int j=i2; j<L; j++) ou<<text[j];}

	return;
;}

