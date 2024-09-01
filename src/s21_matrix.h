#ifndef S21_MATRIX_H_
#define S21_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define OK 0
#define ERROR 1
#define CALCULATION_ERROR 2

#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

// вспомогательная функция для поиска детерминанта матрицы:
// принимает номер столбца (равный номеру строки) и матрицу.
// ищет максимум в столбце ниже главной диагонали.
// если максимум обнаружен не на главной диагонали, то происходит замена тех
// строк, где найден максимум и той, номер которой подался на вход. если замена
// произошла возвращает -1; иначе 1.
int swap_rows_from_main_diagonal(int n, matrix_t *A);
// вспомогательная функция для поиска матрицы алгебраических дополнений:
// принимает матрицу, для которой нужно найти матрицы алгебраических дополнений,
// номер вырезаемой строки, номер вырезаемого столбца, порядок матрицы минора,
// матрицу минора. возвращает минор элемента.
double calc_minor(matrix_t *A, int crossed_out_rows, int crossed_out_columns,
                  int order, matrix_t *minor_element);
// вспомогательная функция для проверки детерминанта:
// принимает матрицу, возвращает её детерминант.
double check_determinant(matrix_t *A);
// вспомогательная функция для мечати матрицы:
// принимает на вход указатель на матрицу, которую нужно распечатать.
void printf_matrix(matrix_t *result);

#endif  // S21_MATRIX_H_
