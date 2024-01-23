#pragma once
#include <bits/stdc++.h>
#include <random>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <inttypes.h>
#include <assert.h>
#include <errno.h>
#include <semaphore.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <queue>
#include <cstdio>//
#include <limits.h>
#include <iomanip>
#include <stdexcept>//
#include "BlockManager.h"


using namespace std;
#define BYTE 8
#define PAGESIZE 4096
#define SPACE_CHAR ' '
#define TAB_CHAR ' '
#define SPLIT_CHAR ' '
#define BRANCH 5
#define PI 3.14159

#define GETMASK(index, size) (((1 << (size)) - 1) << (index))
#define READFROM(data, index, size) (((data) & GETMASK((index), (size))) >> (index))
#define WRITETO(data, index, size, value) ((data) = ((data) & (~GETMASK((index), (size)))) | ((value) << (index)))
#define FIELD(data, name, index, size) \
  inline decltype(data) name() { return READFROM(data, index, size); } \
  inline void set_##name(decltype(data) value) { WRITETO(data, index, size, value); }

extern int M;//dimension number
extern int N;//data set size
extern int n;//sample size
extern int dtype;
extern int independent;
extern vector<int> d;//cardinality
extern vector<int> dbit;//the bit number of cardinality
extern int dbit_sum;
extern int page_capacity;//the tuple number in a page
extern int subset_num;
extern string depend_pattern;
extern int bucket_num;
extern int B;
extern vector<int> logical_size0;
extern vector<int> logical_size;
extern int lnum_max;
extern int batch;
extern vector<int> shape;
extern vector<int> empty_tuple;
typedef pair<double, double> point;
typedef vector<point> points;
static const int CHILDS = 3;
static const double INF = 999999999;

extern int bitsperkey;
extern BlockManager* block1;
extern int fcurpageid;
extern char* sdata1;
extern int beginbyte1;///write the filter file
extern vector<vector<int>> filter_offset;
extern vector<string> filters;
extern int last_validchunk;




///------------------------------------------------------------------------------------------------------------------------------------------
void read_Cardinality(const char* datafolder){
    string data_folder(datafolder);
    string scardpath = data_folder + "cardinality.txt";
    const char* cardpath = scardpath.c_str();
    ifstream fin(cardpath);
	string line;
	d.clear();
	dbit.clear();
	dbit_sum=0;
	if (fin)
	{
		while (getline(fin, line))
		{
			istringstream iss(line);
			string temp;
			while (getline(iss, temp, SPACE_CHAR)) {
                d.push_back(stoi(temp));
                double a = log(stof(temp)) / log(2);
                dbit.push_back((int)ceil(a));
                dbit_sum += dbit[dbit.size()-1];
			}
        }
    }
    fin.clear();
    fin.close();
    page_capacity = PAGESIZE * BYTE / dbit_sum;
    cout<<dbit_sum<<endl;
    cout<<page_capacity<<endl;
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
int compare_Twotuples(vector<int> a, vector<int> b){
    int i;
    int same = 0;
    for(i = 0; i < a.size(); i++){
        if(a[i]+1==b[i]) {same++;continue;}///the two tuples are continous
        if(a[i]!=b[i]) return 0;
    }
    if(same == 1) return 2;
    else if(same > 1) return 0;
    return 1;///the two tuples are the same
}
///------------------------------------------------------------------------------------------------------------------------------------------
void splitString(const string& s, vector<string>& tokens, char delim) {
	tokens.clear();
	auto string_find_first_not = [s, delim](size_t pos = 0) -> size_t {
		for (size_t i = pos; i < s.size(); i++) {
			if (s[i] != delim) return i;
		}
		return string::npos;
	};
	size_t lastPos = string_find_first_not(0);
	size_t pos = s.find(delim, lastPos);
	while (lastPos != string::npos) {
		tokens.emplace_back(s.substr(lastPos, pos - lastPos));
		lastPos = string_find_first_not(pos);
		pos = s.find(delim, lastPos);
	}
	return;
}
///------------------------------------------------------------------------------------------------//
void loadQuery(const char* querypath, vector<vector<int>> &qarray) {
	ifstream fin(querypath);
	string line;
	if (fin)
	{
		while (getline(fin, line))
		{
			istringstream iss(line);
			string temp;
			vector<int> sv;
			while (getline(iss, temp, ' ')) {
				sv.push_back(stoi(temp));
				//if (qnum < 1) cout << stoi(temp) << endl;
			}
			qarray.push_back(sv);
			//qnum++;
		}
	}
	//cout<<qarray[0][0]<<" "<<qarray[0][1]<<" "<<qarray[0][2]<<" "<<qarray[0][3]<<endl;
	fin.clear();
	fin.close();
}
///------------------------------------------------------------------------------------------------//
void strmncpy(char* s, int start1, int len, char* t, int start2){
    for(int i = 0; i < len; i++)
        t[start2+i] = s[start1+i];
    return;
}