#include "common.h"
#include "bfilter.h"
#include "rfilter.h"


int main()
{
    string datafolder = argv[1];
    const char* dataFolder = datafolder.c_str();
    string datapath = dataf + "data.txt";
    const char* dataPath = datapath.c_str();
    string querypath = dataf + "query.txt";
    const char* queryPath = querypath.c_str();
    string binarypath1 = dataf + "binary1.txt";
    const char* binaryPath1 = binarypath1.c_str();
    string offsetpath = dataf + "offset.txt";
    const char* offsetPath = offsetpath.c_str();
    string filterpath = dataf + "filter.txt";
    const char* filterPath = filterpath.c_str();
    string resultpath1 = dataf + "result1.txt";
    const char* resultPath1 = resultpath1.c_str();
    string resultpath2 = dataf + "result2.txt";
    const char* resultPath2 = resultpath2.c_str();

    int i, j;
    M = atoi(argv[2]);
	for(i = 0; i < M; i++){ 
		d.push_back(atoi(argv[3+i]));
		logical_size.push_back(atoi(argv[3+M+i]));
	}
    for(i = 0; i < M; i++)
        logical_size0.push_back((int)pow(ceil(sqrt(logical_size[i])), 2));
    dbit_sum = 0;
    for(i = 0; i < d.size(); i++){
        double a = log(d[i]) / log(2);
        dbit.push_back((int)ceil(a));
        dbit_sum += dbit[dbit.size()-1];
    }
    page_capacity = PAGESIZE * BYTE / dbit_sum;
    for(i = 0; i < M; i++){
        empty_tuple.push_back(0);
    }
	//build the filter 
    Rfilter* rf = new Rfilter();
    rf->construct_Rangefilter(dataPath, binaryPath1, filterPath, offsetPath);
    rf->process_Queries(binaryPath1,queryPath,offsetPath,filterPath,resultPath1);

    cout << "Ends!" << endl;
    return 0;
}
