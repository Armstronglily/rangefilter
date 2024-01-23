#pragma once
#include "common.h"
#include "BlockManager.h"
#include "QuickSort.h"
#include "Timer.h"
#include "bloom_filter.hpp"
#include "basic.hpp"
#include "bfilter.h"



class Rfilter
{
public:

    int m;
    vector<int> dimbit;///the bits of cardinality on each dimension
    vector<int> chunkbit;///the bits of chunk length (logical_size) on each dimension
    vector<int> piecesbit;///the bits of chunk number on each dimension
    vector<int> powpieces;
    vector<int> intervallen;///the length of partition interval on each dim
    vector<int> intervalbit;
    int chunknum;///total chunk number
    int cknum = 0;///non-empty chunk number
    int chunk_logsize;///fixed logical size
    vector<int> chunksize;///tuple number (pysical size) of each chunk
    vector<int> page_startid;
    vector<int> page_endid;
    vector<int> filterbyte;
    char* sdata;///write the data file

    int bitmapbytes;///the byte number of bitmap on each chunk
    char* rfbitmap;///the bitmap filter of each chunk
    vector<vector<vector<int>>> oneDrange;///all the possible one-dimensional ranges for each value on each dim
    vector<int> oneDrange_bit;///the bit of one-dimensional range number on each dim
    int oneDrangebits;///the sum of oneDrange_bit
    uint64_t allmranges;///the number of all possible multi-dimensional ranges
    vector<vector<vector<uint16_t>>> mulDrange;
    vector<string> filters;





    Rfilter(){
        int i;
        int j = 0, k = 0;
        m = d.size();
        dimbit.resize(m);
        chunkbit.resize(m);
        piecesbit.resize(m);
        powpieces.resize(m);
        intervallen.resize(m);
        intervalbit.resize(m);
        oneDrange.resize(m);
        oneDrange_bit.resize(m);
        oneDrangebits = 0;
        sdata = new char[PAGESIZE];
        sdata1 = new char[PAGESIZE];
        beginbyte1 = 0;
        for(i = 0; i < m; i++){
            dimbit[i] = (int)ceil(log(d[i]) / log(2));
            chunkbit[i] = (int)ceil(log(logical_size[i]) / log(2));
            k += chunkbit[i];
            piecesbit[i] = dimbit[i] - chunkbit[i];
            j += piecesbit[i];
            powpieces[i] = pow(2,piecesbit[i]);
            intervallen[i] = (int)ceil(sqrt(logical_size[i]));
            intervalbit[i] = (int)ceil(log(intervallen[i]) / log(2));
        }
        chunknum = pow(2, j);
        chunk_logsize = pow(2, k);
        mulDrange.resize(chunk_logsize);
        chunksize.resize(chunknum, 0);
        page_startid.resize(chunknum, -1);
        page_endid.resize(chunknum, -1);
        filterbyte.resize(chunknum);
        filter_offset.resize(chunknum);
        filters.resize(chunknum);
        for(i = 0; i < chunknum; i++) filter_offset[i].resize(5);
        for(i = 0; i < m; i++)
            oneDrange_bit[i] = ceil(log(3.5 * logical_size0[i] - 1.5 * intervallen[i] - 2) / log(2));
            oneDrangebits += oneDrange_bit[i];
        }
        allmranges = pow(2, oneDrangebits);

    }



    void transfer_Txt_ToBinaryfile(const char* datapath, const char* binarypath);
    void transfer_Tuple_ToBinary(vector<int> tup, unsigned char &c, int &beginbyte, int &beginbit);
    void read_Page(BlockManager* block, int pid, vector<vector<int>> &tuples);
    void read_Offsets(const char* offsetpath);

    void construct_Rangefilter(const char* datapath, const char* binarypath, const char* filterpath, const char* offsetpath);
    void compute_Total1Drange();///one-dimensional ranges
    void compute_Rangeset(int chunkid, vector<vector<int>> alltuples, vector<uint64_t> &rangeIDs);///multi-dimensional ranges, build bitmap as filter
    void write_RFbitmap(int chunkid);

    void read_Filters(const char* offsetpath, const char* filterpath);///read the filters from file
    int search_RFbitmap(int chunkid, uint64_t arange);///Identify whether a range is empty on rfbitmap
    void get_MulRanges4query(vector<int> aquery, vector<uint64_t> &mranges);
    void process_Queries(const char* binarypath, const char* querypath, const char* offsetpath, const char* filterpath, const char* resultpath);



};





