#include "histSaver.h"
#include "HISTFITTER.h"
#include "EigenVectorCalc.h"
int main(int argc, char const *argv[])
{
	float a[3][3] = {{1,2,3},{4,6,7},{3,2,4}};
	float **aa = (float**)malloc(3*sizeof(float*));
	float *aaa[3];
	int b = 3;
	float *c;
	for (int i = 0; i < 3; ++i)
	{
		aaa[i] = (float*)malloc(3*sizeof(float));
		for (int j = 0; j < 3; ++j)
		{
			aaa[i][j] = a[i][j];
		}
		aa[i] = aaa[i];
	}
	EigenVectorCalc(aa,b,c);
	return 0;
}
