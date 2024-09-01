#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (rows < 1 || columns < 1 || result == NULL) {
    return ERROR;
  }
  result->rows = rows;
  result->columns = columns;
  result->matrix = (double **)calloc(
      rows * columns * sizeof(double) + rows * sizeof(double), sizeof(double));
  double *pointer = (double *)(result->matrix + rows);
  for (int i = 0; i < rows; i++) {
    result->matrix[i] = pointer + columns * i;
  }
  return OK;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->columns > 0 && A->rows > 0 && A->matrix != NULL && A != NULL) {
    free(A->matrix);
    A->rows = 0;
    A->columns = 0;
    A->matrix = NULL;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (A->rows != B->rows || A->columns != B->columns) {
    return FAILURE;
  }
  int status = SUCCESS;
  for (int i = 0; i < A->rows && status == SUCCESS; i++) {
    for (int j = 0; j < A->columns && status == SUCCESS; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
        status = FAILURE;
      }
    }
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || B->rows < 1 ||
      B->columns < 1 || B->matrix == NULL || result == NULL) {
    return ERROR;
  }
  if (A->rows != B->rows || A->columns != B->columns) {
    return CALCULATION_ERROR;
  }
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return OK;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || B->rows < 1 ||
      B->columns < 1 || B->matrix == NULL || result == NULL) {
    return ERROR;
  }
  if (A->rows != B->rows || A->columns != B->columns) {
    return CALCULATION_ERROR;
  }
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
  return OK;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || result == NULL) {
    return ERROR;
  }
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }
  return OK;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || B->rows < 1 ||
      B->columns < 1 || B->matrix == NULL || result == NULL) {
    return ERROR;
  }
  if (A->columns != B->rows) {
    return CALCULATION_ERROR;
  }
  s21_create_matrix(A->rows, B->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < B->columns; j++) {
      for (int m = 0; m < A->columns; m++) {
        result->matrix[i][j] += A->matrix[i][m] * B->matrix[m][j];
      }
    }
  }
  return OK;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || result == NULL) {
    return ERROR;
  }
  s21_create_matrix(A->columns, A->rows, result);
  for (int i = 0; i < A->columns; i++) {
    for (int j = 0; j < A->rows; j++) {
      result->matrix[i][j] = A->matrix[j][i];
    }
  }
  return OK;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || result == NULL) {
    return ERROR;
  }
  if (A->columns != A->rows || A->columns == 1) {
    return CALCULATION_ERROR;
  }
  int order = A->rows - 1;
  matrix_t minor_element = {0, 0, 0};
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] =
          pow(-1, i + j) * calc_minor(A, i, j, order, &minor_element);
    }
  }
  return OK;
}

int s21_determinant(matrix_t *A, double *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || result == NULL) {
    return ERROR;
  }
  if (A->columns != A->rows) {
    return CALCULATION_ERROR;
  }
  *result = 1;
  int n = 0;
  for (int j = 0; j < A->columns; j++) {
    int sign = swap_rows_from_main_diagonal(j, A);
    for (int i = A->rows - 1; i > n; i--) {
      double number = A->matrix[i][j] / A->matrix[j][j];
      for (int m = 0; m < A->columns; m++) {
        A->matrix[i][m] -= number * A->matrix[j][m];
      }
    }
    *result *= sign * A->matrix[j][j];
    n++;
  }
  if (fabs(*result) <= 1e-6) {
    *result = fabs(*result);
  }
  return OK;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (A->rows < 1 || A->columns < 1 || A->matrix == NULL || result == NULL) {
    return ERROR;
  }
  if (A->columns != A->rows || A->columns == 1 ||
      fabs(check_determinant(A)) <= 1e-6) {
    return CALCULATION_ERROR;
  }
  matrix_t matrix_complements = {0, 0, 0};
  matrix_t matrix_complements_transpose = {0, 0, 0};
  s21_calc_complements(A, &matrix_complements);
  s21_transpose(&matrix_complements, &matrix_complements_transpose);
  s21_mult_number(&matrix_complements_transpose, (1 / check_determinant(A)),
                  result);
  s21_remove_matrix(&matrix_complements);
  s21_remove_matrix(&matrix_complements_transpose);
  return OK;
}

double check_determinant(matrix_t *A) {
  matrix_t B = {0, 0, 0};
  s21_create_matrix(A->rows, A->columns, &B);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      B.matrix[i][j] = A->matrix[i][j];
    }
  }
  double det = 0;
  s21_determinant(&B, &det);
  s21_remove_matrix(&B);
  return det;
}

double calc_minor(matrix_t *A, int crossed_out_rows, int crossed_out_columns,
                  int order, matrix_t *minor_element) {
  double minor = 0;
  s21_create_matrix(order, order, minor_element);
  for (int n = 0, i = 0; i < A->rows; i++, n++) {
    if (i != crossed_out_rows) {
      for (int m = 0, j = 0; j < A->columns; j++, m++) {
        if (j != crossed_out_columns) {
          minor_element->matrix[n][m] = A->matrix[i][j];
        } else {
          m--;
        }
      }
    } else {
      n--;
    }
  }
  s21_determinant(minor_element, &minor);
  s21_remove_matrix(minor_element);
  return minor;
}

int swap_rows_from_main_diagonal(int m, matrix_t *A) {
  int flag = 1;
  double max = fabs(A->matrix[m][m]);
  int max_rows = 0;
  for (int i = m; i < A->rows; i++) {
    if (max < fabs(A->matrix[i][m])) {
      max = fabs(A->matrix[i][m]);
      max_rows = i;
      flag = -1;
    }
  }
  if (flag == -1) {
    for (int j = 0; j < A->rows; j++) {
      double tmp = A->matrix[m][j];
      A->matrix[m][j] = A->matrix[max_rows][j];
      A->matrix[max_rows][j] = tmp;
    }
  }
  return flag;
}

void printf_matrix(matrix_t *result) {
  printf("\nrows = %d\ncolumns = %d\n", result->rows, result->columns);
  for (int i = 0; i < result->rows; i++) {
    printf("\n");
    for (int j = 0; j < result->columns; j++) {
      printf("%lf ", result->matrix[i][j]);
    }
  }
}
