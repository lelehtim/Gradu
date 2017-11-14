#include<stdio.h>
int main() {
  FILE *f;
  f=fopen("test","w");
  fclose(f);
}
