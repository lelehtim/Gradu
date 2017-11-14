#include <stdio.h>
void main() {
  int i=0;
  int a=(int*) malloc(sizeof(int)) ;
  for (int j=0; j<64; j++) {
    for (int k=0; k<64; k++) {
      for (int l=0; l<64; l++) {
	i++;
      }
    }
  }
  printf("%d\n",i);
}
