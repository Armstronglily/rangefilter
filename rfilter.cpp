#include "rfilter.h"



///--------------------------------------------------------------------------------------------------------------------
void Rfilter::transfer_Txt_ToBinaryfile(const char* datapath, const char* binarypath){
    int i, j, g;
    ///read data, the first time
    string line;
    ifstream fin(datapath);
    if (fin)
    {
        while (getline(fin, line))
        {
            istringstream iss(line);
            string temp;
            vector<int> tup;
            while (getline(iss, temp, TAB_CHAR)) {
                tup.push_back(stoi(temp));
            }
            int chunkid = 0;
            chunkid = tup[0] / logical_size[0];
            for(i = 1; i < m; i++){
                chunkid = chunkid << piecesbit[i];
                chunkid += (tup[i] / logical_size[i]);
            }
            chunksize[chunkid]++;
        }
    }
    fin.clear();
    fin.close();
    int former = 0;
    cknum = 0;
    for(i = 0; i < chunknum; i++){
        if(chunksize[i] == 0) continue;
        page_startid[i] = former;
        former += ceil((double)chunksize[i] / page_capacity);
        page_endid[i] = former - 1;
        cknum++;
        last_validchunk = i;
    }
    ///read data the second time, and write the chunks
    vector<vector<int>> buffer[chunknum];
    int beginbyte = 0, beginbit = 0;
    unsigned char c = 0;
    sdata[0] = '\0';
    vector<int> page_currentid;
    page_currentid = page_startid;
//    memcpy(page_currentid, page_startid, sizeof(page_currentid));
    BlockManager* block = new BlockManager(binarypath, O_CREAT, PAGESIZE);
    ifstream fin1(datapath);
    if (fin1)
    {
        while (getline(fin1, line))
        {
            istringstream iss(line);
            string temp;
            vector<int> tup;
            while (getline(iss, temp, TAB_CHAR)) {
                tup.push_back(stoi(temp));
            }
            int chunkid = 0;
            chunkid = tup[0] / logical_size[0];
            for(i = 1; i < m; i++){
                chunkid = chunkid << piecesbit[i];
                chunkid += tup[i] / logical_size[i];
            }
            buffer[chunkid].push_back(tup);
            if(buffer[chunkid].size() == page_capacity){
                for(j = 0; j < page_capacity; j++){
                    transfer_Tuple_ToBinary(buffer[chunkid][j], c, beginbyte, beginbit);
                }
                block->WriteBlock(sdata, page_currentid[chunkid]++, PAGESIZE);
                sdata[0] = '\0'; beginbyte = 0; beginbit = 0; c = 0;
                buffer[chunkid].clear();
            }
        }
    }
    fin1.clear();
    fin1.close();
    ///write the left
    for(i = 0; i < chunknum; i++){
        if(buffer[i].size() == 0) continue;
        for(g = 0; g < buffer[i].size(); g++){
            transfer_Tuple_ToBinary(buffer[i][g], c, beginbyte, beginbit);
        }
        for(g < buffer[i].size(); g < page_capacity; g++){
            transfer_Tuple_ToBinary(empty_tuple, c, beginbyte, beginbit);
        }
        block->WriteBlock(sdata, page_currentid[i]++, PAGESIZE);
        sdata[0] = '\0'; beginbyte = 0; beginbit = 0; c = 0;
        buffer[i].clear();
    }
    delete block;
    return;
}

///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::transfer_Tuple_ToBinary(vector<int> tup, unsigned char &c, int &beginbyte, int &beginbit){
    int i, j;
    int length;
    for (j = 0; j < M; j++) {//对tuple的每一维取后len[j]位
        length = dbit[j];
        while (beginbit + length >= BYTE) {
            c = (c << (BYTE - beginbit)) + READFROM(tup[j], length - BYTE + beginbit, BYTE - beginbit);
            length -= (BYTE - beginbit);
            sdata[beginbyte] = c;
            beginbyte++;
            c = 0;
            beginbit = 0;
        }
        if (length != 0) {
            c = (c << length) + READFROM(tup[j], 0, length);
            beginbit += length;
        }
    }
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::read_Page(BlockManager* block, int pid, vector<vector<int>> &tuples){
    int k, g;
    unsigned char c = 0;
    sdata[0] = '\0';
	int beginbyte = 0, beginbit = 0;
	int length = 0;
    block->ReadBlock(sdata, pid, PAGESIZE);
    for (k = 0; k < page_capacity; k++) { //每个线性化值写到arr中，再存入struct
		vector<int> arr(M);
		beginbyte = k * dbit_sum / BYTE;
		beginbit = k * dbit_sum % BYTE;
		for (g = 0; g < M; g++) {
			length = dbit[g];
			while (BYTE - beginbit <= length) {
				c = *(sdata + beginbyte);
				arr[g] = arr[g] << (BYTE - beginbit);
				arr[g] += READFROM(c, 0, BYTE - beginbit);//从右侧第0位向左读BYTE - beginbit位;
				length -= BYTE - beginbit;
				beginbit = 0;
				beginbyte++;
			}
			if (length != 0) {
				c = *(sdata + beginbyte);
				arr[g] = arr[g] << length;
				arr[g] += READFROM(c, BYTE - beginbit - length, length);
				beginbit += length;
			}
		}
		if(compare_Twotuples(arr, empty_tuple)==1) break;
		tuples.push_back(arr);
	}
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::read_Offsets(const char* offsetpath){
    int i, j;
    ifstream fin(offsetpath);
    string line;
    if (fin)
    {
        while (getline(fin, line))
        {
            istringstream iss(line);
            string temp;
            vector<int> a;
            while (getline(iss, temp, ' ')) {
                a.push_back(stoi(temp));
            }
            page_startid[a[0]] = a[1];
            page_endid[a[0]] = a[2];
            //filter_offset[a[0]].resize(5);
            filter_offset[a[0]][0] = a[3];///filter type
            filter_offset[a[0]][1] = a[4];///filter beginning page
            filter_offset[a[0]][2] = a[5];///......byte
            filter_offset[a[0]][3] = a[6];///filter ending page
            filter_offset[a[0]][4] = a[7];///......byte
            last_validchunk = a[0];
            if(a[1]==a[2]) filterbyte[a[0]] = a[7] - a[5];
            else filterbyte[a[0]] = PAGESIZE * (a[2]-a[1]) + PAGESIZE - a[5] + 1 + a[7];   
		}
    }
    fin.clear();
    fin.close();
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::construct_Rangefilter(const char* datapath, const char* binarypath, const char* filterpath, const char* offsetpath){
    int i, j, k, g;
    int currentid;
    ///partition array to chunks and generate binary file
    transfer_Txt_ToBinaryfile(datapath, binarypath);
    ofstream fout(offsetpath, ios::out);
    vector<vector<int>> tuples;//the read tuples of a page
    vector<vector<int>> tuplesintotal;//the total tuples of a chunk
    int beginbyte = 0, beginbit = 0;
    unsigned char c = 0;
    sdata[0] = '\0';  //data
    sdata1[0] = '\0'; //filter
    fcurpageid = 0; beginbyte1 = 0;
    BlockManager* block = new BlockManager(binarypath, O_CREAT, PAGESIZE);
    block1 = new BlockManager(filterpath, O_CREAT, PAGESIZE);

    ///compute one-dimensional ranges
    compute_Total1Drange();
    for(i = 0; i < chunknum; i++){
        if(chunksize[i] == 0) continue;
        currentid = page_startid[i];
        for(g = page_startid[i]; g <= page_endid[i]; g++){
            read_Page(block, g, tuples);
            tuplesintotal.insert(tuplesintotal.end(), tuples.begin(), tuples.end());
            tuples.clear();
        }
        ///compute range-set (multidimensional ranges) and construct range filter
        vector<uint64_t> rangeids;
        compute_Rangeset(i, tuplesintotal, rangeids);
        if(allmranges / rangeids.size() >= bitsperkey) { ///build bloom filter based on bitsperkey
            cout<<"chunkid "<<i<<" bloom filter"<<endl;
            Bfilter *bf = new Bfilter(rangeids.size());
            bf->construct_Bloomfilter(i, rangeids);
            delete bf;
        }
        else{
            write_RFbitmap(i);
        }
        tuplesintotal.clear();
        fout<<i<<" "<<page_startid[i]<<" "<<page_endid[i]<<" "<<filter_offset[i][0]<<" "<<filter_offset[i][1]<<" "<<filter_offset[i][2]<<" "<<filter_offset[i][3]<<" "<<filter_offset[i][4]<<endl;
    }
    fout.clear();
    fout.close();
    return;
}

///------------------------------------------------------------------------------------------------------------------------------------------
///compute all the one-dimensional ranges for each value on each dim
void Rfilter::compute_Total1Drange(){
    int i, j, k, g;
    vector<int> intervalnum(m);
    for(i = 0; i < m; i++){
        oneDrange[i].resize(logical_size[i]);
        //intervalnum[i] = (int)((double)logical_size[i] / intervallen[i]);
        int offset = 3 * logical_size0[i] - 2 * intervallen[i] - 2;///the number of the first two kinds ranges
        for (j = 0; j < logical_size[i]; j++){
            ///1. [a,b], a and b are not partition values
            int low, high;
            int interval = j / intervallen[i];
            if(((j + 1) % intervallen[i]) == 0) {///j is partition value
                if(j == intervallen[i] - 1) {///the first partition value
                    low = intervallen[i] - 2;
                    high = low + 2 * (intervallen[i] - 1) - 2;
                }
                else{
                    if(j == logical_size0[i] - 1){///the last partition value
                        high = 2 * (intervallen[i] - 1) * ((j + 1) / intervallen[i]) - 3;
                        low = high - intervallen[i] + 2;
                    }
                    else{///other partition values
                            //cout<<"j = "<<j<<endl;
                        low = 2 * (intervallen[i] - 1) * interval - 2 + intervallen[i] - 1;
                        high = low + 2 * intervallen[i] - 3;
                    }
                }
            }
            else{///j is not partition value
                if(j < intervallen[i]){///fall into the first interval
                    if(j==0){
                        low = 0; high = intervallen[i] - 3;
                    }
                    else{
                        low = j - 1; high = low + intervallen[i] - 1 - 1;
                    }
                }
                else{///other intervals
                    low = 2 * (intervallen[i] - 1) * interval - 2 + j % intervallen[i];
                    high = low + intervallen[i] - 1;
                }
            }
            for(k = low; k <= high; k++) oneDrange[i][j].push_back(k);
            ///2. [a,a]
            oneDrange[i][j].push_back((uint8_t)(2*(logical_size0[i] - intervallen[i] -1) + j));
            ///3. [a,b], a and b are partition values
            int sum = offset, beginp;
            for(k = 0; k < interval + 1; k++){//row
                sum = sum + intervallen[i] - k;
                beginp = sum - (intervallen[i] - interval);
                for(g = interval; g < intervallen[i]; g++){
                    oneDrange[i][j].push_back( (uint8_t)beginp++);
                }
            }
        }//for j

    }//for i
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
///compute multidimensional ranges
void Rfilter::compute_Rangeset(int chunkid, vector<vector<int>> alltuples, vector<uint64_t> &rangeIDs){
    int i, j, g, k;
    int mbyte, mbit;
    bitmapbytes = (int)ceil((double)allmranges / BYTE);
    rfbitmap = new char[bitmapbytes];
    for(i = 0; i < bitmapbytes; i++) rfbitmap[i] = 0;
    ///calculate chunk coordinates
    vector<int> coords(m);
    int cid = chunkid;
    for(j = m-1; j >=0; j--){
        coords[j] = cid & (powpieces[j] - 1);
        cid  = cid >> piecesbit[j];
    }
    cout<<"chunkid "<<chunkid<<", tuple num "<<alltuples.size()<<endl;
    for(j = 0; j < alltuples.size(); j++){
        ///compute the offsets of the tuple on each dim
        vector<int> tupoffset(m);
        int tupleid;
        int rangenum = 1;
        for(k = 0; k < m; k++) {
            tupoffset[k] = alltuples[j][k] - coords[k] * logical_size[k];
            rangenum *= oneDrange[k][tupoffset[k]].size();
        }
        vector<uint64_t> mranges(rangenum);
        vector<uint64_t> lengths(m, 1);
        lengths[0] = 1;
        for(k = 1; k < m; k++) {
            lengths[k] = lengths[k-1] * oneDrange[k-1][tupoffset[k-1]].size();
        }
        ///compute the multi-dimensional ranges
        for(g = 0; g < rangenum; g++) {
            mranges[g] = oneDrange[0][tupoffset[0]][g % oneDrange[0][tupoffset[0]].size()];
        }
        for(k = 1; k < m; k++) {
            for(g = 0; g < rangenum; g++){
				mranges[g] = mranges[g] << oneDrange_bit[k];
                mranges[g] += oneDrange[k][tupoffset[k]][g / lengths[k] % oneDrange[k][tupoffset[k]].size()];
            }
            //cout<<""<<endl;
        }
        for(g = 0; g < rangenum; g++){
            mbyte = mranges[g] / BYTE;
            mbit = mranges[g] % BYTE;
            char a = rfbitmap[mbyte];
            a = a >> (7-mbit);
            int f = a & 1;
            if(f != 1){
                rfbitmap[mbyte] = (rfbitmap[mbyte] | (1 << (7-mbit)));
            }
        }

    }///tuple
    for(i = 0; i < bitmapbytes; i++){
        char ar = rfbitmap[i];//cout<<"c = "<<ar+0;
        for(j = 0; j < 8; j++){
            if(((ar >> (7-j)) & 1)==1) {
                rangeIDs.push_back(i*8+j);
            }
        }
    }
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::write_RFbitmap(int chunkid){
    int i;char v;
    int offset = 0, offset1;///the offset on rfbitmap
    filter_offset[chunkid][0] = 0;///0 represents bitmap as filter
    filter_offset[chunkid][1] = fcurpageid;
    filter_offset[chunkid][2] = beginbyte1;
    cout<<"chunkid = "<<chunkid<<", bitmapbytes = "<<bitmapbytes<<endl;
    if(PAGESIZE - beginbyte1 > bitmapbytes){
        strmncpy(rfbitmap, 0, bitmapbytes, sdata1, beginbyte1);
        beginbyte1 += bitmapbytes;
    }
    else{
        strmncpy(rfbitmap, 0, PAGESIZE - beginbyte1, sdata1, beginbyte1);
        block1->WriteBlock(sdata1, fcurpageid++, PAGESIZE);
        offset += (PAGESIZE - beginbyte1);
        offset1 = offset;
        beginbyte1 = 0;
        sdata1[0] = '\0';
        for(int k = 0; k < (bitmapbytes - offset1) / PAGESIZE; k++){
            strmncpy(rfbitmap, offset, PAGESIZE, sdata1, beginbyte1);
            block1->WriteBlock(sdata1, fcurpageid++, PAGESIZE);
            beginbyte1 = 0;
            sdata1[0] = '\0';
            offset += PAGESIZE;
        }
        if(offset != bitmapbytes){
            strmncpy(rfbitmap, offset, bitmapbytes-offset, sdata1, beginbyte1);
            beginbyte1 = bitmapbytes-offset;
        }
    }
    filter_offset[chunkid][3] = fcurpageid;
    filter_offset[chunkid][4] = beginbyte1;///the next beginning byte
    if(chunkid==last_validchunk) block1->WriteBlock(sdata1, fcurpageid++, PAGESIZE);
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
int Rfilter::search_RFbitmap(int chunkid, uint64_t arange){///Identify whether a range is empty on bitmap
    int i, j, k;
    int byte = arange / BYTE;
    int bit = arange % BYTE;
    char a = filters[chunkid][byte];
    if(((a >> (7-bit)) & 1) == 1) return 1;
    return 0;
}

///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::read_Filters(const char* offsetpath, const char* filterpath){///read the filters from file
    int i, j;
    sdata[0] = '\0';
    int sign = 0;
    read_Offsets(offsetpath);
    BlockManager* block = new BlockManager(filterpath, O_CREAT, PAGESIZE);

    for(i = 0; i < chunknum; i++){//chunknum
            cout<<"read chunk "<<i<<endl;
        if(page_startid[i]==-1) continue;
        int length = 0, offset = 0;
        if(filter_offset[i][1] == filter_offset[i][3])
            length = filter_offset[i][4] - filter_offset[i][2];
        else{
            length = PAGESIZE - filter_offset[i][2];
            for(j = filter_offset[i][1]+1; j <= filter_offset[i][3]; j++){
                length += PAGESIZE;
            }
            length += filter_offset[i][4];
        }
        char* meta = new char[length];
        meta[0] = '\0';
        if(sign == 0) block->ReadBlock(sdata, filter_offset[i][1], PAGESIZE);///need to read the first page
        if(filter_offset[i][1] == filter_offset[i][3]){///begin page is also end page
            strmncpy(sdata, filter_offset[i][2], filter_offset[i][4] - filter_offset[i][2], meta, offset);
            offset = offset + filter_offset[i][4] - filter_offset[i][2];
        }
        else{
            strmncpy(sdata, filter_offset[i][2], PAGESIZE - filter_offset[i][2], meta, offset);
            offset = offset + PAGESIZE - filter_offset[i][2];
            sdata[0] = '\0';
            for(j = filter_offset[i][1]+1; j < filter_offset[i][3]; j++){
                block->ReadBlock(sdata, j, PAGESIZE);
                strmncpy(sdata, 0, PAGESIZE, meta, offset);
                offset += PAGESIZE;
                sdata[0] = '\0';
            }
            block->ReadBlock(sdata, j, PAGESIZE);
            strmncpy(sdata, 0, filter_offset[i][4], meta, offset);
            sign = 1;
        }
        string afilter(meta, meta+length);
        filters[i] = afilter;
    }
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::process_Queries(const char* binarypath, const char* querypath, const char* offsetpath, const char* filterpath, const char* resultpath){
    int i, j, k, g;
    int lowchunk, highchunk;
    vector<int> p1(m);
    vector<int> p2(m);
    vector<int> coor(m);
    vector<int> q(2*m);
    int nemptychunks, borderchunks, overlapchunks;
    int cid, inrange, isborderchunk;
    double filtertime, processtime;
    Timer* timer = new Timer();

    Bfilter* bf = new Bfilter(100);
    ofstream fout(resultpath, ios::out);
    BlockManager* block = new BlockManager(binarypath, O_CREAT, PAGESIZE);
    read_Filters(offsetpath, filterpath);
    vector<vector<int>> query;
    loadQuery(querypath, query);
    for(i = 0; i < query.size(); i++){
        lowchunk = 0, highchunk = 0;
        vector<int> chunklist;
        for(j = 0; j < m; j++){
            p1[j] = query[i][2*j] / logical_size[j];
            lowchunk = (lowchunk << piecesbit[j]) + p1[j];
            p2[j] = query[i][2*j+1] / logical_size[j];
            highchunk = (highchunk << piecesbit[j]) + p2[j];
        }
        ///compute the overlapping chunks, find the border chunks and filter them
        borderchunks = 0, overlapchunks = 0;
        timer->Start();
        for(k = lowchunk; k <= highchunk; k++){
            cid = k;
            inrange = 1;
            isborderchunk = 0;
            for(j = m-1; j >= 0; j--){
                coor[j] = cid & (powpieces[j] - 1);
                cid = cid >> piecesbit[j];
                if(coor[j]<p1[j] || coor[j]>p2[j]) {
                    inrange = 0;
                    break;
                }
                if(coor[j]==p1[j] || coor[j]==p2[j]) isborderchunk = 1;
            }
            if(inrange==0) continue;
            overlapchunks++;
            if(isborderchunk==0) continue;
            borderchunks++;///overlapping border chunks, not considering empty or non-empty-----///
            ///search range on the chunk's filter
            for(j = 0; j < m; j++){
                q[2*j] = max(query[i][2*j], coor[j]*logical_size[j]) - coor[j]*logical_size[j];
                q[2*j+1] = min(query[i][2*j+1], (coor[j] + 1) * logical_size[j] - 1) - coor[j] * logical_size[j];
            }
            vector<uint64_t> mranges;
            get_MulRanges4query(q, mranges);
            ///filtering
            if(filter_offset[k][0]==0){///bitmap
                for(j = 0; j < mranges.size(); j++){
                    if(search_RFbitmap(k, mranges[j])==1){
                        chunklist.push_back(k);
                        break;
                    }
                }
            }
            else{///bloom filter
                for(j = 0; j < mranges.size(); j++){
                    if(bf->search_Bloomfilter(filters[k], mranges[j], filterbyte[k])==1){
                        chunklist.push_back(k);
                        break;
                    }
                }
            }
        }//chunk
        filtertime = timer->GetRunningTime();
        timer->Stop();
        ///read and output non-empty chunks
        vector<vector<int>> tuples;
        int sign1, sign2;
        nemptychunks = 0;///real non-empty chunks
        timer->Start();
        for(k = 0; k < chunklist.size(); k++){
            for(j = page_startid[chunklist[k]]; j <= page_endid[chunklist[k]]; j++){//page
                read_Page(block, j, tuples);
                sign1 = 0;
                for(g = 0; g < tuples.size(); g++){
                    sign2 = 1;
                    for(int v = 0; v < m; v++){
                        if(tuples[g][v] < query[i][2*v] || tuples[g][v] > query[i][2*v+1]) {
                            sign2 = 0; break;
                        }
                    }
                    if(sign2==1) {sign1=1; break;}
                }
                if(sign1==1) {nemptychunks++; break;}
            }
        }
        processtime = timer->GetRunningTime();
        timer->Stop();
        fout<<overlapchunks<<" "<<borderchunks<<" "<<filtertime<<" "<<processtime<<endl;
    }
    fout.clear();
    fout.close();
    return;
}
///------------------------------------------------------------------------------------------------------------------------------------------
void Rfilter::get_MulRanges4query(vector<int> aquery, vector<uint64_t> &mranges){
    int i, j, k, g;
    int interval1, interval2, offset1, offset2, sum, p1, p2, num1, num2;
    vector<vector<int>> oneranges(m);
    ///get one-dimensional ranges
    for(i = 0; i < aquery.size()/2; i++){
        interval1 = aquery[2*i] / intervallen[i];
        interval2 = aquery[2*i+1] / intervallen[i];
        offset1 = aquery[2*i] % intervallen[i];
        offset2 = aquery[2*i+1] % intervallen[i];
        p1 = (aquery[2*i]+1) / intervallen[i];
        p2 = (aquery[2*i+1]+1) / intervallen[i];
        num1 = 2*(logical_size0[i] - intervallen[i] - 1);
        num2 = 3*logical_size0[i] - 2*intervallen[i] - 2;
        if(aquery[2*i] == aquery[2*i+1]){///[a,a]
            oneranges[i].push_back(num1 + aquery[2*i]);
            continue;
        }
        if(aquery[2*i]==0 || offset1 == intervallen[i]-1){///[a,b], a is partition value
            if(offset2 == intervallen[i]-1){///b is a partition value
                sum = 0;
                for(j = 0; j < p1; j++)///row
                    sum = sum + intervallen[i] - j;
                oneranges[i].push_back(num2 + sum + p2 - p1 - 1);
            }
            else{///b is not a partition value
                if(aquery[2*i+1] - aquery[2*i] < intervallen[i]){ ///b is near, distance is smaller than interval length
                    if(aquery[2*i]==0) {
                        oneranges[i].push_back(offset2 - 1);
                    }
                    else {
                        oneranges[i].push_back(2*(intervallen[i]-1) * p1 - 2 + offset2);
                    }
                }
                else{///b is far
                    sum = 0;
                    for(j = 0; j < p1; j++)///row
                        sum = sum + intervallen[i] - j;
                    oneranges[i].push_back(num2 + sum + p2 - p1 - 1);
                    oneranges[i].push_back(2*(intervallen[i]-1) * interval2 - 2 + offset2);
                }
            }
            continue;
        }
        if(offset2 == intervallen[i]-1){///[a,b], b is a partition value
            if(aquery[2*i+1] - aquery[2*i] < intervallen[i]){///a is near
                oneranges[i].push_back(2*(intervallen[i]-1) * interval1 - 2 + intervallen[i]-1 + offset1);
            }
            else{///a is far
                sum = 0;
                for(j = 0; j <= p1; j++)///row
                    sum = sum + intervallen[i] - j;
                oneranges[i].push_back(num2 + sum + interval2 - interval1 - 1);
                oneranges[i].push_back(2*(intervallen[i]-1) * interval1 - 2 + intervallen[i]-1 + offset1);
            }
            continue;
        }
        ///[a,b], a and b are not partition values
        if(interval1==interval2){
            for(j = offset1; j <= offset2; j++){
                oneranges[i].push_back(num1 + aquery[2*i]+j);
            }
        }
        else{
            if(interval1 + 1 == interval2){///adjcent intervals
                oneranges[i].push_back(2*(intervallen[i]-1) * interval1 - 2 + intervallen[i]-1 + offset1);
                oneranges[i].push_back(2*(intervallen[i]-1) * interval2 - 2 + offset2);
            }
            else{
                oneranges[i].push_back(2*(intervallen[i]-1) * interval1 - 2 + intervallen[i]-1 + offset1);
                sum = 0;
                for(j = 0; j <= p1; j++)///row
                    sum = sum + intervallen[i] - j;
                oneranges[i].push_back(num2 + sum + interval2 - interval1 - 2);
                oneranges[i].push_back(2*(intervallen[i]-1) * interval2 - 2 + offset2);
            }
        }
    }
    ///combine multidimensional ranges
    vector<uint64_t> lengths(m, 1);
    lengths[0] = 1;
    for(k = 1; k < m; k++) {
        lengths[k] = lengths[k-1] * oneranges[k-1].size();
    }
    int rnum = lengths[m-1] * oneranges[m-1].size();
    mranges.resize(rnum);
    ///compute the multi-dimensional ranges
    for(g = 0; g < rnum; g++) {
        mranges[g] = oneranges[0][g % oneranges[0].size()];/////////////
    }
    for(k = 1; k < m; k++) {
        for(g = 0; g < rnum; g++){
            mranges[g] = mranges[g] << oneDrange_bit[k];//////////////
            mranges[g] += oneranges[k][g / lengths[k] % oneranges[k].size()];//////////////
        }
    }
    return;
}

