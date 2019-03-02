#include "histSaver.h"
#include "HISTFITTER.h"
#include "EigenVectorCalc.h"
int main(int argc, char const *argv[])
{
	float a[3][3] = {{1,2,3},{4,6,7},{3,2,4}};
	float **aa = (float**)malloc(3*sizeof(float*));
	float *aaa[3];
	float **dd = (float**)malloc(3*sizeof(float*));
	int b = 3;
	float *c = (float*)malloc(3*sizeof(float));
	for (int i = 0; i < 3; ++i)
	{
		dd[i] = (float*)malloc(3*sizeof(float));
		aaa[i] = (float*)malloc(3*sizeof(float));
		for (int j = 0; j < 3; ++j)
		{
			aaa[i][j] = a[i][j];
		}
		aa[i] = aaa[i];
	}
	printf("start evaluate\n");
	EigenVectorCalc(aa,b,c,dd);
	for (int i = 0; i < 3; ++i)
	{
		printf("Matrix: %f %f %f\n", a[0][i] , a[1][i], a[2][i]);
	}
	for (int i = 0; i < 3; ++i)
	{
		printf("Eigen Value: %f, EigenVector: (%f %f %f)\n", c[i], dd[i][0] , dd[i][1], dd[i][2]);
	}
	return 0;
}
