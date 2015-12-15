////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////
////////////                       GrROOT
////////////
////////////          Purpose:
////////////                   To assist in the analysis of data from
////////////                 the gretina/S800 experimental setup.
////////////                          
////////////          Current Maintainers:
////////////                 Kathrin Wimmer  (wimmer@phys.s.u-tokyo.ac.jp)
////////////                 Eric Lunderberg (lunderberg@nscl.msu.edu)
////////////
////////////          Distribution:
////////////                   Please do not redistribute this software directly.
////////////                   If someone new wants a copy of this software,
////////////                 email one of the maintainers for the download link.
////////////                   This allows us to keep track of who has the software,
////////////                 and properly distribute updates and bug fixes.
////////////                 
////////////          Suggestions:
////////////                   We view the development of the software as a collaborative
////////////                 effort, and as such, welcome and appreciate any suggestions
////////////                 for bug fixes and improvements.
////////////
////////////          Disclaimer:
////////////                 This software is provided as-is, with no warranty.
////////////                 No current or future support is guaranteed for this software.
////////////
////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <stdlib.h>

using namespace std;

void printlist(int* list, int size)
{
  for (int i=0;i<size;i++) cout << list[i];
  cout <<endl;
}
void swap(int *list,int i1, int i2) // swap two elements of list
{
  int tmp = list[i2];
  list[i2] = list[i1];
  list[i1] = tmp;
}

void sort(int * list, int size) // numerically sort the list
{
  for (int i=0;i<size;i++)
    if (list[i] < list[0]) swap(list,0,i);
  if (size>1) sort(&list[1],size-1);
}
void get_nth_perm(int *list, int size, int n)
{
  if (size <=1 || n==0) return;
  sort(list,size);

  int i,X=1;
  for(i=1;i<size;i++) X*=i;

  for (i=n/X;i>0;i--) swap(list, i,i-1);

  get_nth_perm(&list[1],size-1,n%X);
}
bool get_next_perm(int *list, int size)
{
  int i,j;
  if (size < 2) return false; // false means done
  for(i=size-1; i>0 && list[i] < list[i-1]; i--);
  if (i==0) return false; // false means done (already on last element)
  sort(&list[i],size-i);
  i--;
  for(j=i+1;j<size && list[i]>list[j]; j++);
  swap(list,i,j);
  return true; // true means there are more permutations to be had
}
