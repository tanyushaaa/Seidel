#include "sm.h"

float epsilon = 0.01;

float norm(float* vec1, float* vec2, int n){
	float sum = 0;
	for (int i = 0; i < n; i++){
		sum += fabs(vec1[i] - vec2[i]);
	}
	return sum;
}

std::pair<float**, int> new_mtrx(std::pair<float**, int> m, float* x){
	float** massiv = new float* [m.second];
	for (int i = 0; i < m.second; i++){
		massiv[i] = new float [m.second + 1];
	}

	float divisor = 0;
	for (int i = 0; i < m.second; i++){ 		//��������� ������ ������ �������� ������� ���������
		divisor = m.first[i][i]; 				//�� ������� ����� ���������� ����-��
		for (int j = 0; j < m.second + 1; j++){ //��������� (� �� � 0 ��� ��-��� ������� ���������) ��������� ����� (������ ����)
			massiv[i][j] = (i == j) ? 0 :			//����� ��������� ������� ��������� ������ (��������� ����)
							((j == m.second) ?
							m.first[i][j] / divisor : //��� ����� �� ������������ ������ ��-�,
						  - m.first[i][j] / divisor); //����� ������ �������� �������
		}
		x[i] = massiv[i][m.second];
	}
	return std::make_pair(massiv, m.second);
}

float* iteration(std::pair<float**, int> m, float* x){
	float* x_copy = new float [m.second];
	for (int i = 0; i < m.second; i++){
		x_copy[i] = x[i];
	}

	int n = m.second;
	float* X = static_cast<float*>(calloc(n, sizeof(float)));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n + 1; j++){
			if (j != n)
				X[i] += m.first[i][j] * x[j];
			else
				X[i] += m.first[i][j];
		}
		x[i] = X[i];
	}

	if (norm(x_copy, X, m.second) > epsilon)
		return iteration(m, X);
	return X;
}

float* sm(const Matrix& mtrx){ //����� � � ������ ������� �������� ������� �������� �� �==0 � ���������==0
	float* x = new float [mtrx.GetA().second];
	float* answer = new float [mtrx.GetA().second];
	answer = iteration(new_mtrx(mtrx.GetA(), x), x);
//	for (int i = 0; i < mtrx.GetA().second; i++){
//		std::cout << answer[i] << " ";
//	}
	return answer;
}
