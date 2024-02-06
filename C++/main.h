#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <math.h>
#include <vector>
#include <tuple>
#include <chrono>

using namespace std;

#define PI						3.14159265358979323846
#define JACOBI_ERROR			1.0e-9

long Get_Time();
vector<vector<double>> read_file(const string& file_name, const char& delimiter);