#include <iostream>
#include <vector>
#include <iomanip>
#include "./include/cppjieba/Jieba.hpp"
#include "./include/simhash/Simhasher.hpp"

using namespace std;


bool Loaded_context(string Path, string& context);
bool Output_result(string Path,double res);
unsigned hamdist( string str1,  string str2);
bool cal_HamSim(double & res, int dis);


const char* const DICT_PATH = "include/dict/jieba.dict.utf8";
const char* const HMM_PATH = "./include/dict/hmm_model.utf8";
const char* const USER_DICT_PATH = "./include/dict/user.dict.utf8";
const char* const IDF_PATH = "./include/dict/idf.utf8";
const char* const STOP_WORD_PATH = "./include/dict/stop_words.utf8";
const char* const HMM_MODEL_PATH = "./include/dict/hmm_model.utf8";
const char* const KEYWORDS_1_PATH = "./keyword1.txt";
const char* const KEYWORDS_2_PATH = "./keyword2.txt";


string context_1;
string context_2;
vector<pair<string ,double>> result_1;
vector<pair<string ,double>> result_2;

uint64_t simHashKey_contxt_1 = 0;
uint64_t simHashKey_contxt_2 = 0;

size_t topN = 1000;

int main(int argc, char * argv[]){

    const char* const SRC_1_PATH = "src1.txt";
    const char* const SRC_2_PATH = "src2.txt";
    const char* const RES_PATH = "res.txt";


    Loaded_context(SRC_1_PATH,context_1);
    Loaded_context(SRC_2_PATH,context_2);


    simhash::Simhasher simhasher(DICT_PATH, HMM_MODEL_PATH, IDF_PATH, STOP_WORD_PATH);


    simhasher.extract(context_1, result_1, topN);
    simhasher.extract(context_2, result_2, topN);


    simhasher.make(context_1, topN, simHashKey_contxt_1);
    simhasher.make(context_2, topN, simHashKey_contxt_2);


    string s1;
    string s2;
    simhash::Simhasher::toBinaryString(simHashKey_contxt_1,s1);
    cout <<"Simhash_context_1: " << s1 << endl;
    simhash::Simhasher::toBinaryString(simHashKey_contxt_2,s2);
    cout <<"Simhash_context_2: "<< s2 << endl;


    int dis = simhash::Simhasher::isEqual(simHashKey_contxt_1, simHashKey_contxt_2);
    cout <<"hamming: "<< dis << endl;
    double res;
    cal_HamSim(res,dis);
    cout<<fixed<<setprecision(2);
    cout <<"similarity: " << res << endl;

    Output_result(RES_PATH, res);

    return 0;

}

bool Loaded_context(string Path, string& Context){
    ifstream context_ifs(Path);
    if(context_ifs.fail()){
        cout << "error: SRC_1_PATH null!" << endl;
        context_ifs.close();
        return false;
    }else{
        while(!context_ifs.eof()){
            string temp;
            context_ifs >> temp;
            Context.append(temp);
        }
        context_ifs.close();
        cout << "File path: " << Path << " loaded " << endl;
        return 1;
    };


}

bool Output_result(string Path, double res){
    ofstream of(Path,ios::app);
    ofstream of_res_1(KEYWORDS_1_PATH);
    ofstream of_res_2(KEYWORDS_2_PATH);

    if(of.fail() || of_res_1.fail() || of_res_2.fail()){
        cout << "error: RES_PATH null!" << endl;
        of.close();
        of_res_1.close();
        of_res_2.close();
        return false;
    }else{

        of << "本次结果: " << res << endl;
        of_res_1 << "前1000个关键词: " <<result_1 << endl;
        of_res_2 << "前1000个关键词: " <<result_2 << endl;
        of.close();
        of_res_1.close();
        of_res_2.close();
        cout << "File path: " << Path << " output " << endl;
        return 1;
    }

}

unsigned hamdist(string str1, string str2)
{
    if (str1.empty() || str2.empty()) {
        return 0;
    }
    if (str1.length() != str2.length()){
        return 0;
    }
    unsigned dist = 0;
    unsigned i = 0;
    while(i < str1.length()){
        dist += (str1[i] != str2[i]) ? 1 : 0;
        i++;
    }
    return dist;
}


bool cal_HamSim(double & res, int dis){
    res = 0.01 * (100 - dis * 100 / 128.0);
    return 1;
}


double getMold(const vector<double>& vec)
{
    int n = vec.size();
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
        sum += vec[i] * vec[i];
    return sqrt(sum);
}
double cos_distance(const vector<double>& base, const vector<double>& target)
{
    int n = base.size();
    assert(n == target.size());
    double tmp = 0.0;
    for (int i = 0; i < n; ++i)
        tmp += base[i] * target[i];
    double simility =  tmp / (getMold(base)*getMold(target));
    return simility;
}
