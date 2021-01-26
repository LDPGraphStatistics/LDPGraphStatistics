#ifndef MEMOPE
#define MEMOPE
#include <stdlib.h>
#include <malloc.h>

template < typename Type >
void malloc1D(Type **source, int i)
{
	*source = (Type*)calloc(i, sizeof(Type));
	if (source == NULL)
	{
		printf("malloc error.");
		exit(-1);
	}
}

template < typename Type >
void malloc2D(Type ***source, int i, int j)
{
	Type ** ret = (Type**)calloc(i, sizeof(Type*));
	for (int counter = 0; counter < i; counter++)
		malloc1D(&(ret[counter]), j);
	*source = ret;
}

template < typename Type >
void malloc3D(Type ****source, int i, int j, int k)
{
	Type *** ret = (Type***)calloc(i, sizeof(Type**));
	for (int counter = 0; counter < i; counter++)
		malloc2D(&(ret[counter]), j, k);

	*source = ret;
}

template <typename Type>
void free1D(Type *source)
{
	free(source);

	source = NULL;
}

template < typename Type >
void free2D(Type **source, int i)
{
	for (int counter = 0; counter < i; counter++)
		free(source[counter]);
	free(source);

	source = NULL;
}

template < typename Type >
void free3D(Type **source, int i, int j)
{
	for (int counter = 0; counter < i; counter++)
		free2D(source[counter], j);
	free(source);

	source = NULL;
}

#endif
