#include <stdio.h>
#include <stdlib.h>

int main(void)
{
	FILE *file;
	char text[] = "Testing";
	
	file = fopen("test.txt","wb");
	fwrite(text, sizeof(text), 1, file);
	return 1;
}		