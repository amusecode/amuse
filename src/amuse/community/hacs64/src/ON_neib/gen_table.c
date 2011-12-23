/* generate 7-bit table for morton key */
#include <stdio.h>
int main(){
	int i;
	for(i=0; i<128; i++){
		int mask = 64;
		printf("0");
		while(mask){
			printf("%d", (i&mask) ? 1 : 0);
			mask /= 2;
		}
		printf(",\n");
	}
	return 0;
}
