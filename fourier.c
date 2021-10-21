#include <math.h>
#include "fourier.h"


/* function that analises that if its the column or line*/
void fft_vector_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height, int forward, int line)
{

    int element;
    int dimension_1, dimension_2;

    /* check if its the line of the matrix*/
    if (line)
    {
        dimension_1 = width;
        dimension_2 = height;
    }
    /* if not, its the column of the matrix*/
    else
    {
        dimension_1 = height;
        dimension_2 = width;
    }
    
    double complex new_vector[dimension_1]; /* create a new vector */
    double complex new_vector_answer[dimension_1]; /* create the vector that will receive the answer */

    /* iterates dimension_2*/
    for (int vector = 0; vector < dimension_2; vector++)
    {   
        /* iterates the vector */
        for (element = 0; element < dimension_1; element++)
        {   
            /* generates the new vector */
            if (line)
            {
                new_vector[element] = matrix[vector][element];
            }
            else
            {
                new_vector[element] = matrix[element][vector];
            }
        }

        /* aplies the forward or inverse transformation to the line */
        if (forward) 
        {
            fft_forward(new_vector, new_vector_answer, dimension_1);
        }
        else
        {
            fft_inverse(new_vector, new_vector_answer, dimension_1);
        }

        /* reinsert the values in the original matrix*/
        for (element = 0; element < dimension_1; element++)
        {
            if (line)
            {
                matrix[vector][element] = new_vector_answer[element];
            }
            else
            {
                matrix[element][vector] = new_vector_answer[element];
            }
        }
    }
}

/* function that aplies the bidimensional transform */
void fft_main_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height, int forward)
{
    /* line */
    fft_vector_2d(matrix, width, height, forward, 1);
    /* column */
    fft_vector_2d(matrix, width, height, forward, 0);
}

/* separate the vector in even or odd lists*/
void split_list(double complex main_list[], double complex even_list[], double complex odd_list[], int n)
{
    int counter = 0;
    for (int index = 0; index < n; index += 2)
    {
        even_list[counter] = main_list[index];
        odd_list[counter] = main_list[index + 1];
        counter++;
    }
}

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign)
{
    for (int k = 0; k < n; k++)
    {
        t[k] = 0;

        for (int j = 0; j < n; j++)
        {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n)
{
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n)
{
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++)
    {
        s[k] /= n;
    }
}


void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign)
{ 
    /* create the even and odd vectors*/
    double complex even_list_s[n / 2];
    double complex odd_list_s[n / 2];

    /* create the even and odd vectors*/
    double complex even_list_t[n / 2];
    double complex odd_list_t[n / 2];

    /* if n equals 1, k is between 0 and n-1, that is, k = 0*/
    if (n == 1)
    {
        t[0] = s[0];
        return;
    }

    /* split the list into the vectors*/
    split_list(s, even_list_s, odd_list_s, n);

    /*Fourier transform of even list*/
    fft(even_list_s, even_list_t, n / 2, sign);
    /*Fourier transform of odd list*/
    fft(odd_list_s, odd_list_t, n / 2, sign);

    for (int k = 0; k < n / 2; k++)
    {
        t[k] = even_list_t[k] + odd_list_t[k] * cexp(sign * 2 * PI * k * I / n);
        t[k + n / 2] = even_list_t[k] - odd_list_t[k] * cexp(sign * 2 * PI * k * I / n);
    }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n)
{
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n)
{
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++)
    {
        s[k] /= n;
    }
}

/* function that aplies the bidimensional transform accoring to direction*/
void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height)
{
    fft_main_2d(matrix, width, height, 1);
}
/* function that aplies the inverse bidimensional transform accoring to direction*/
void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height)
{
    fft_main_2d(matrix, width, height, 0);
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip)
{
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip)
            {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height)
{
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height)
{
    filter(input, output, width, height, 1);
}
