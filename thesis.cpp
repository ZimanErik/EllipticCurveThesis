#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <NTL/ZZ.h>

using json = nlohmann::json;
using namespace NTL;
using namespace std;


#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <nlohmann/json.hpp>

#include <exception>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <streambuf>

NTL_CLIENT

class NullStreamBuf : public std::streambuf {
    int overflow(int ch) override {
        return ch; 
    }
};
   
std::streambuf* originalCoutBuf = std::cout.rdbuf();

void turnOffCout() {
    static NullStreamBuf nullStreamBuf;  
    cout.rdbuf(&nullStreamBuf);
}

void turnOnCout() {
    cout.rdbuf(originalCoutBuf);
}
   

struct ECPoint {
    ZZ x, y;
    bool isPointAtInfinity;

    ECPoint() : isPointAtInfinity(true) {}  
    ECPoint(ZZ x, ZZ y) : x(x), y(y), isPointAtInfinity(false) {}

    void print() const {
        if (isPointAtInfinity) cout << "IPoint(∞)" << endl;
        else cout << "(" << x << ", " << y << ")" << endl;
    }

    string toString() const {
        ostringstream oss;
        if (isPointAtInfinity)
            oss << "IPoint(∞)";
        else
            oss << "(" << x << ", " << y << ")";
        return oss.str();
    }

    bool operator==(const ECPoint& other) const {
        if (isPointAtInfinity && other.isPointAtInfinity) return true;
        if (isPointAtInfinity || other.isPointAtInfinity) return false;
        return (x == other.x) && (y == other.y);
    }
};

struct ZZHasher {
    size_t operator()(const ZZ& num) const {
        stringstream ss;
        ss << num;  // Convert ZZ to string
        string str = ss.str();
        return std::hash<string>{}(str);  // Hash the string representation of ZZ
    }
};

struct ECPointHasher
{
  size_t operator()(const ECPoint& P) const
  {
    size_t rowHash = ZZHasher{}(P.x);
    size_t colHash = ZZHasher{}(P.y) << 1;
    return rowHash ^ colHash;
  }
};


struct EllipticCurve {
    ZZ a, b, p;

    EllipticCurve(ZZ a, ZZ b, ZZ p) : a(a), b(b), p(p) {}

    bool isOnCurve(const ECPoint& P) const {

        if (P.isPointAtInfinity) return true;
        ZZ_p::init(p);
        
        ZZ_p x = to_ZZ_p(P.x);
        ZZ_p y = to_ZZ_p(P.y);
        
        ZZ_p left = sqr(y);
        
        ZZ_p right = x*x*x + to_ZZ_p(a)*x + to_ZZ_p(b);      
        
        return (left == right);
    }

    void print() const {
        cout << "Elliptic curve : y^2 = x^3 + " << a << "x + "<< b << " mod "<< p << endl;
    }

    string toString() const {
        ostringstream oss;
        oss << "Elliptic curve : y^2 = x^3 + " << a << "x + "<< b << " mod "<< p ;
        return oss.str();
    }
};

class PointNotOnCurveException : public exception {
private:
    string message;

public:
    PointNotOnCurveException(const ECPoint& P, const EllipticCurve& EC){
        ostringstream mess;
        mess <<"Point ("<<P.x<<","<<P.y<<") is not on curve y^2=x^3+"<<EC.a<<"x"<<"+"<<EC.b<<"mod "<<EC.p; 
        message = mess.str();
    }

    const char* what() const throw()
    {
        return message.c_str();
    }
};

ECPoint doublePoint(const ECPoint& P,const EllipticCurve& EC){
   
    if(P.y==0) return ECPoint();
    ZZ_p::init(EC.p);
    ZZ_p p_inv = inv(to_ZZ_p(2)*to_ZZ_p(P.y));

    // in reality the second '*' bellow is division because of inverse multiplication 
    ZZ_p m = ( (to_ZZ_p(3) * to_ZZ_p(sqr(P.x)) + to_ZZ_p(EC.a) ) * (p_inv) );    
    ZZ_p n = m * m - to_ZZ_p(2) * to_ZZ_p(P.x);

    ZZ x = rep(n);
    ZZ y = rep( m*(to_ZZ_p(P.x) - n) - to_ZZ_p(P.y));
    return ECPoint(x,y);
}

ECPoint addPoints(const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC){

    try {
        if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P,EC);
        if (!EC.isOnCurve(Q)) throw PointNotOnCurveException(Q,EC);
    }
    catch (PointNotOnCurveException& e) {
        cout << "Caught an exception: " << e.what() << endl;
        return ECPoint(); 
    }


    if (P.isPointAtInfinity) return Q; 
    if (Q.isPointAtInfinity) return P; 

    if (P.x == Q.x && (P.y != Q.y || P.y == ZZ(0)))
        return ECPoint(); // (P + (-P) = Identity point)

    if (P==Q){ // P = Q -->  P + Q = P + P 
        return doublePoint(P,EC);
    }

    ZZ_p::init(EC.p);

    ZZ_p P_x = to_ZZ_p(P.x);
    ZZ_p P_y = to_ZZ_p(P.y);

    ZZ_p Q_x = to_ZZ_p(Q.x);
    ZZ_p Q_y = to_ZZ_p(Q.y);

    ZZ_p lambdaNumerator = Q_y - P_y;
    ZZ_p lambdaDenominator = inv(Q_x - P_x);

    ZZ_p lambda = lambdaNumerator * lambdaDenominator;

    ZZ_p x_p = sqr(lambda) - P_x - Q_x ;
    ZZ x = rep(x_p);
    ZZ y = rep( lambda*(P_x - x_p) - P_y);

   return ECPoint(x,y);
}

// Scalar multiplication: cP
ECPoint scalarPointMultiplication(const ZZ& c, const ECPoint& P, const EllipticCurve& EC) {
    try {
        if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P,EC);
    }
    catch (PointNotOnCurveException& e) {
        cout << "Caught an exception: " << e.what() << endl;
        return ECPoint(); 
    }

    ECPoint result; 

    // done so we have modifiable lvalue for shifting c's bits and still be able to call function with ZZ(123)
    ZZ k = ZZ(c);  

    ECPoint base = P;
    turnOffCout();
    while (k > 0) {
        if (bit(k, 0)) { // check current LSB 
            result = addPoints(result, base, EC);
            cout << result.toString() << " | ";
        }
        base = doublePoint(base,EC);
        k >>= 1;
    }
    turnOnCout();
    return result;
}

ECPoint inversePoint(const ECPoint& P, const EllipticCurve& EC) {

    try {
        if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P,EC);
    }
    catch (PointNotOnCurveException& e) {
        cout << "Caught an exception: " << e.what() << endl;
        return ECPoint();
    }

    if (P == ECPoint()) return ECPoint();
    if (P.y == ZZ(0)) return P;

    ZZ_p y = to_ZZ_p(P.y);
    y = -y;
    return ECPoint(P.x,rep(y));   
}

ZZ bruteSearchSubgroupOrder(const ECPoint& P, const EllipticCurve& EC){
    if (!EC.isOnCurve(P)) {
        throw PointNotOnCurveException(P, EC);
        return (ZZ(0));
    }

    cout << "brute search order of generated subgroup from point " << P.toString() <<  endl;
    ECPoint temp = ECPoint();
    ZZ order = ZZ(0);
    do{
        temp = addPoints(temp,P,EC);
        order++;
    } while (!(temp == ECPoint()));
    return order;
}

ZZ naiveSearchECDLP(const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC, const ZZ& n){

    if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P, EC);
    if (!EC.isOnCurve(Q)) throw PointNotOnCurveException(P, EC);

    ZZ order = ZZ(0);           
    ECPoint temp = ECPoint();        
    
    while (order <= n) {
        if (temp == Q) {
            break;
        }
        temp = addPoints(temp, P, EC); 
        order++;
    }
    return order;

}

ZZ babyStepGiantStepECDLP(const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC, const ZZ& n){

    ZZ m = SqrRoot(n) + ZZ(1);
    unordered_map<ECPoint, ZZ, ECPointHasher> babySteps;

    // Compute baby steps: P, 2P, 3P, ..., jP
    ECPoint currentBabyStep = P;
    for (ZZ j = ZZ(1); j <= m; ++j) {
        babySteps[currentBabyStep] = j;
        currentBabyStep = addPoints(currentBabyStep, P, EC); 
    }

    ECPoint minusMP = inversePoint(scalarPointMultiplication(m, P, EC), EC);
    ECPoint giantStepComponent = ECPoint();
    ECPoint currentGiantStepBase = Q; 
    for (ZZ i = ZZ(0); i <= m; ++i) {
        ECPoint curr = currentGiantStepBase; 

        if (babySteps.find(curr) != babySteps.end()){
            ZZ j = babySteps[curr];
            ZZ k = (i * m + j) % n;
            return k;              
        }

        currentGiantStepBase = addPoints(currentGiantStepBase, minusMP, EC); 
    }
    return ZZ(-1); 
}


struct PollardRhoState {
    ECPoint X;
    ZZ a, b;
};

ZZ partition(const ECPoint& P) {
    return P.x % ZZ(3); 
}

PollardRhoState f(const PollardRhoState& state, const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC, const ZZ& n) {
    PollardRhoState nextState = state;
    ZZ category = partition(state.X);

    if (category == 0) {
        nextState.X = addPoints(state.X, P, EC);
        nextState.a = (state.a + 1) % n; 
    } else if (category == 1) {
        nextState.X = doublePoint(state.X,EC);
        nextState.a = (2 * state.a) % n; 
        nextState.b = (2 * state.b) % n; 
    } else {
        nextState.X = addPoints(state.X, Q, EC);
        nextState.b = (state.b + 1) % n; 
    }
    return nextState;
}

ZZ pollardRhoECDLP(const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC, const ZZ& n) {

    PollardRhoState turtle = {P, ZZ(1), ZZ(0)}; 
    PollardRhoState rabbit = {P, ZZ(1), ZZ(0)}; 

    turtle = f(turtle, P, Q, EC, n);
    rabbit = f(f(rabbit, P, Q, EC, n), P, Q, EC, n);

    while (!(turtle.X == rabbit.X)) {
        turtle = f(turtle, P, Q, EC, n);
        rabbit = f(f(rabbit, P, Q, EC, n), P, Q, EC, n);
    }
         
    ZZ numerator = (turtle.a - rabbit.a) % n;
    ZZ denominator = (rabbit.b - turtle.b) % n;

    if (denominator == 0) {
        return ZZ(-1);
    }

    ZZ_p::init(n); 

    ZZ_p p_inv = inv(to_ZZ_p(denominator));

    if (p_inv == 0) {
        cout << "Denominator has no inverse modulo n." << endl;
        return ZZ(-1);
    }

    ZZ k = rep((to_ZZ_p(numerator) * p_inv));

    return k;
        
}

// Function to extract bits from a ZZ number and store them in a vector<bool>
vector<bool> getBitVector(const ZZ& num, long bitSize) {
    vector<bool> bitVec(bitSize, false);
    
    long numBits = NumBits(num); 
    for (long i = 0; i < numBits; i++) {
        bitVec[bitSize - 1 - i] = bit(num, i); 
    }
    
    return bitVec;
}

// Function to align two numbers into bit vectors with equal length
pair<vector<bool>, vector<bool>> alignNumbers(const ZZ& a, const ZZ& b) {
    long maxBits = max(NumBits(a), NumBits(b)); // Determine required bit length
    
    vector<bool> aBits = getBitVector(a, maxBits); 
    vector<bool> bBits = getBitVector(b, maxBits);

    return {aBits, bBits}; 
}

// computing m*P + l*Q
ECPoint shamirsTrick(const ECPoint& P, const ECPoint& Q, const ZZ& m, const ZZ& l, const EllipticCurve& EC) {
    auto [sBits, tBits] = alignNumbers(m, l);
    long bitLength = sBits.size();
    
    ECPoint result = ECPoint();
    ECPoint sumOfPandQ = addPoints(P,Q,EC);

    for (long i = 0; i < bitLength; i++) {
        result = doublePoint(result,EC);

        if (sBits[i] && tBits[i]) {
            result = addPoints(result,sumOfPandQ,EC);
        } else if (sBits[i]) {
            result = addPoints(result,P,EC);
        } else if (tBits[i]) {
            result = addPoints(result,Q,EC);
        }
    }
    
    return result;
}

ECPoint naiveLinearCombination(const ECPoint& P, const ECPoint& Q, const ZZ& m, const ZZ& l, const EllipticCurve& EC) {

    ECPoint mP = scalarPointMultiplication(m,P,EC);
    ECPoint lQ = scalarPointMultiplication(l,Q,EC);
    return addPoints(mP,lQ,EC);
}

// Function to print bit vectors
void printBitVector(const vector<bool>& bits) {
    for (bool bit : bits) {
        cout << bit;
    }
    cout << endl;
}


struct ShamirJSONTestParams {
    ZZ bit_size;
    ZZ prime;
    ZZ curve_a;
    ZZ curve_b;
    ZZ curve_order;
    ZZ P_x;
    ZZ P_y;
    ZZ Q_x;
    ZZ Q_y;
};

struct ECDLP_JSONTestParams {
    ZZ bit_size;
    ZZ prime;
    ZZ curve_a;
    ZZ curve_b;
    ZZ curve_order;
    ZZ generator_x;
    ZZ generator_y;
    ZZ generator_order;
    ZZ log2_order;
};

#include <filesystem>
namespace fs = std::filesystem;

bool loadJson(const string& preferred, const string& fallback, json& outJson) {
    string filename;

    if (fs::exists(preferred)) {
        filename = preferred;
        cout << "Loading from preferred file: " << preferred << endl;
    } else if (fs::exists(fallback)) {
        filename = fallback;
        cout << "Preferred file not found. Loading from fallback: " << fallback << endl;
    } else {
        cerr << "Neither " << preferred << " nor " << fallback << " exist." << endl;
        return false;
    }

    ifstream file(filename);
    if (!file) {
        cerr << "Failed to open JSON file: " << filename << endl;
        return false;
    }

    try {
        file >> outJson;

        if (!outJson.is_array()) {
            cerr << "Expected a JSON array at the top level in " << filename << endl;
            return false;
        }

        return true;
    } catch (const json::parse_error& e) {
        cerr << "JSON parse error in " << filename << ": " << e.what() << endl;
        return false;
    }
}

vector<ShamirJSONTestParams> loadShamirParams(const string& filename) {
    json j;
    if (!loadJson(filename, "shamircurves.json",j)) {
        return {};
    }

    vector<ShamirJSONTestParams> curves;
    try {
        for (const auto& entry : j) { 
            ShamirJSONTestParams cp;
            cp.bit_size = ZZ(NTL::INIT_VAL, entry["bit_size"].get<string>().c_str());
            cp.prime = ZZ(NTL::INIT_VAL, entry["prime"].get<string>().c_str());
            cp.curve_a = ZZ(NTL::INIT_VAL, entry["curve"]["a"].get<string>().c_str());
            cp.curve_b = ZZ(NTL::INIT_VAL, entry["curve"]["b"].get<string>().c_str());
            cp.curve_order = ZZ(NTL::INIT_VAL, entry["curve_order"].get<string>().c_str());
            cp.P_x = ZZ(NTL::INIT_VAL, entry["P"]["x"].get<string>().c_str());
            cp.P_y = ZZ(NTL::INIT_VAL, entry["P"]["y"].get<string>().c_str());
            cp.Q_x = ZZ(NTL::INIT_VAL, entry["Q"]["x"].get<string>().c_str());
            cp.Q_y = ZZ(NTL::INIT_VAL, entry["Q"]["y"].get<string>().c_str());
            curves.push_back(cp);
        }
    } catch (const json::type_error& e) {
        cerr << "Type error while reading Shamir JSON in " << filename << ": " << e.what() << endl;
        return {};
    }
    return curves;
}

vector<ECDLP_JSONTestParams> loadECDLPParams(const string& filename, int loadNaiveCurves) {
    json j;
    if (loadNaiveCurves) {
        if (!loadJson(filename, "ecdlpcurvesnaive.json", j)) {
            return {};
        }
    }
    else {
        if (!loadJson(filename, "ecdlpcurves.json", j)){ 
            return {};
        }
    }

    vector<ECDLP_JSONTestParams> curves;
    try {
        for (const auto& entry : j) { 
            ECDLP_JSONTestParams cp;
            cp.bit_size = ZZ(NTL::INIT_VAL, entry["bit_size"].get<string>().c_str());
            
            cp.prime = ZZ(NTL::INIT_VAL, entry["prime"].get<string>().c_str());
            cp.curve_a = ZZ(NTL::INIT_VAL, entry["curve"]["a"].get<string>().c_str());
            cp.curve_b = ZZ(NTL::INIT_VAL, entry["curve"]["b"].get<string>().c_str());
            cp.curve_order = ZZ(NTL::INIT_VAL, entry["curve_order"].get<string>().c_str());
            cp.generator_x = ZZ(NTL::INIT_VAL, entry["generator"]["x"].get<string>().c_str());
            cp.generator_y = ZZ(NTL::INIT_VAL, entry["generator"]["y"].get<string>().c_str());
            cp.generator_order = ZZ(NTL::INIT_VAL, entry["generator_order"].get<string>().c_str());
            cp.log2_order = ZZ(NTL::INIT_VAL, entry["log2_order"].get<string>().c_str());
            curves.push_back(cp);
        }
    } catch (const json::type_error& e) {
        cerr << "Type error while reading ECDLP JSON in " << filename << ": " << e.what() << endl;
        return {};
    }
    return curves;
}

int runECDLPAnalysis(const string& filename, const string& outputFilename, ZZ (*solveECDLP)(const ECPoint&, const ECPoint&, const EllipticCurve&, const ZZ&), const string& algorithmName) {
    vector<ECDLP_JSONTestParams> curves;
    if (solveECDLP == naiveSearchECDLP){
         curves = loadECDLPParams(filename,1);
    }
    else {
        curves = loadECDLPParams(filename,0);
    }

    if (curves.empty()) {
        cout << "Exiting because JSON file is empty." << endl;
        return 1;
    }

    clock_t start;
    clock_t end;
    vector<double> results;
    stringstream ss;
    ofstream outputFile(outputFilename, std::ios::trunc);

    for (size_t i = 0; i < curves.size(); ++i) {
        const auto& curve = curves[i];

        EllipticCurve EC(curve.curve_a, curve.curve_b, curve.prime);
        ECPoint Point(curve.generator_x, curve.generator_y);
        cout << "Curve #" << i + 1 << ":" << endl;
        cout << "Key bit Size: " << curve.bit_size << endl;
        ZZ privateKey = RandomBnd(curve.curve_order) + 1;
        cout << "Generated private key: " << privateKey << endl;
        ECPoint multipliedPoint = scalarPointMultiplication(privateKey, Point, EC);

        start = clock();
        ZZ result = solveECDLP(Point, multipliedPoint, EC, curve.curve_order);
        end = clock();

        clock_t elapsed = end - start;
        results.push_back(static_cast<double>(elapsed));

        printf("Time measured by %s method: %ld ticks.\n", algorithmName.c_str(), elapsed);
        cout << algorithmName << " result: " << result << endl;
        cout << "Result is the same as generated private key: " << ((privateKey % curve.prime == result % curve.prime) ? "yes" : "no") << endl << endl;

        ss << static_cast<int>(results[i]);
        if (i % 100 != 99) {
            ss << ",";
        } else {
            outputFile << ss.str() << "\n";
            outputFile.flush();
            cout << "results : " << ss.str() << endl;
            ss = stringstream();
        }
    }

    outputFile.close();
    return 0;
}

int runBabyStepGiantStepECDLPanalysis() {
    return runECDLPAnalysis("ecdlpcurves_extra.json", "ECDLP_BSGS.csv", babyStepGiantStepECDLP, "Baby-Step Giant-Step");
}

int runPollardRhoanalysis() {
    return runECDLPAnalysis("ecdlpcurves_extra.json", "ECDLP_PollardRho.csv", pollardRhoECDLP, "Pollard's Rho");
}

int runNaiveSearchECDLPanalysis() {
    return runECDLPAnalysis("ecdlpcurvesnaive_extra.json", "ECDLP_Naive.csv", naiveSearchECDLP, "Naive search");
}

int getIntInput(const string& prompt) {
    string input;
    int value;

    while (true) {
        cout << prompt;
        getline(cin, input);
        stringstream ss(input);

        if (ss >> value && ss.eof()) {
            return value;
        } else {
            cout << "Invalid input. Please enter a whole number (no decimals or letters).\n";
        }
    }
}

void showECDLPSubMenu() {
    bool inSubMenu = true;

    while (inSubMenu) {
        cout << "\n--- Choose ECDLP algorithm to use ---\n";
        cout << "1. BSGS\n";
        cout << "2. Pollard's Rho\n";
        cout << "3. Naive search\n";

        int choice = getIntInput("Enter your choice (1-3): ");

        switch (choice) {
            case 1: {
                cout << "starting ECDLP through BSGS\n";
                runBabyStepGiantStepECDLPanalysis();
                inSubMenu = false;
                break;
            }
            case 2: {
                cout << "starting ECDLP through Pollard's Rho\n";
                runPollardRhoanalysis();
                inSubMenu = false;
                break;
            }
            case 3: {
                cout << "starting ECDLP through Naive search\n";
                runNaiveSearchECDLPanalysis();
                inSubMenu = false;
                break;
            }
            default: {
                cout << "Invalid option. Please choose 1, 2, or 3.\n";
                break;
            }
        }
    }
}

int calculateLinearCombinationOfPoints(const string& filename, const string& outputFilename, ECPoint (*calculateLinearCombination)(const ECPoint& P, const ECPoint& Q, const ZZ& m, const ZZ& l, const EllipticCurve& EC)){
    vector<ShamirJSONTestParams> curves = loadShamirParams(filename);
    if (curves.empty()) {
        cout << "Exiting because JSON file is empty." << endl;
        return 1;
    }

    clock_t start;
    clock_t end;
    clock_t startbrute;
    stringstream ss;
    vector<clock_t> vec_elapsed;
    ofstream outputFile(outputFilename, std::ios::trunc);
    int cnt = 0;
    for (size_t i = 0; i < curves.size()/1000; i++) {
        clock_t sumofelapsed = 0;
        for (int j=0;j<1000;j++){
            const auto& curve = curves[cnt++];
            cout << "Curve #" << cnt << endl;
            EllipticCurve EC(curve.curve_a, curve.curve_b, curve.prime);
            ECPoint P = ECPoint(curve.P_x, curve.P_y);
            ECPoint Q = ECPoint(curve.Q_x, curve.Q_y);
            ZZ n = RandomBnd(curve.curve_order) + 1;
            ZZ m = RandomBnd(curve.curve_order) + 1;
            cout << "Executing Shamirs trick on points " << n << "*" << P.toString() << " and " << m << "*" << Q.toString() << endl ;
            cout << "Bit size of curve order : " << curve.bit_size << endl;
            cout << "Curve order : " << curve.curve_order << endl;
            start = clock();
            ECPoint resultPoint = calculateLinearCombination(P,Q,n,m,EC);
            end = clock();
            
            clock_t elapsed = end - start;
            
            printf("Time measured: %ld ticks.\n", elapsed);

            ss << static_cast<int>(elapsed);
            ss << ",";
            sumofelapsed += elapsed;
        
        }
        outputFile << ss.str() << "\n";
        ss = stringstream();
        vec_elapsed.push_back(sumofelapsed);
    }

    outputFile.close();

    cout << endl;
    
    return 0;
}

int calculateNaiveLinearCombination() {
    return calculateLinearCombinationOfPoints("shamircurves_extra.json", "shamir_resultsnaive.csv", naiveLinearCombination);
}

int calculateShamirTrick() {
    return calculateLinearCombinationOfPoints("shamircurves_extra.json", "shamir_results.csv", shamirsTrick);
}

void showLinearCombinationOfPointsMenu() {
    bool inSubMenu = true;

    while (inSubMenu) {
        cout << "\n--- Choose algorithm calculate linear combination of two points ---\n";
        cout << "1. Shamir's trick\n";
        cout << "2. Naive search\n";

        int choice = getIntInput("Enter your choice (1-2): ");

        switch (choice) {
            case 1: {
                cout << "starting to calculate linear combination of two points using Shamir's trick\n";
                calculateShamirTrick();
                inSubMenu = false;
                break;
            }
            case 2: {
                cout << "starting to calculate linear combination of two points using Naive search\n";
                calculateNaiveLinearCombination();

                inSubMenu = false;
                break;
            }
            default: {
                cout << "Invalid option. Please choose 1 or 2.\n";
                break;
            }
        }
    }
}


int main() {

    bool running = true;

    while (running) {
        cout << "\n=== Main Menu ===\n";
        cout << "1. run tests for calculating linear combination of two points\n";
        cout << "2. run tests for solving ECDLP instances\n";
        cout << "3. Exit\n";

        int mainChoice = getIntInput("Enter your choice (1-3): ");

        switch (mainChoice) {
            case 1: {
                showLinearCombinationOfPointsMenu();
                running = false;
                cout << "Exiting application. Goodbye!\n";
                break;
            }
            case 2: {
                showECDLPSubMenu();
                running = false;
                cout << "Exiting application. Goodbye!\n";
                break;
            }
            case 3: {
                running = false;
                cout << "Exiting application. Goodbye!\n";
                break;
            }
            default: {
                cout << "Invalid option. Please choose 1, 2, or 3.\n";
                break;
            }
        }
    }

    return 0;
}
